#
#   hleaps based on Doug Bates's approach and some of his code
#
#   Author: Mark Banghart, Arun Srinivasan
#
#   examples of calls
#       out  <- hleaps(y~x1+x2+x3+x1:x2+x2:x3, method="r2", altOut=TRUE)
#
#       out2  <- hleaps(y~x1+x2+x3+x4+x5+x6+(x1+x2+x3+x4+x5+x6)^2, data=data2,
#                    nbest=5, minSize=2, maxSize=6, timeOut=300)
#
#   Issues to be completed    
#       - Test when a column (but not a term) is NA 
#       - support models with no intercept or block dropping the intercept
#       - profile code to optimize the search loop
#       - limit search to only model sizes that are in the min to max range
#       - Apply a bounded algorithm based on SSMod values
#       - consider other indicators besides order (a 1:n ordering) like percentile
#
#
######################################################################

#
#' @export

hleaps  <- function (formula, 
                     data,
                     method = "adjr2",   # critria to be reported
                     nbest = 3,            # number of best models per size
                     minSize = NULL,       # smallest size of model reported
                     maxSize = NULL,       # largest size of model reported
                     timeOut = NULL,       # number of seconds to search 
                     #  for models
                     altOut = FALSE,       # alternate display for output
                     mem = 2, # Amount of memory available for use (GB)
                     
                     ...) {
  
  
  ##  hpmodsetup 
  ##
  ##  Does all of what hpmods does except create the matrix of models
  ##
  ##  Modified from Doug Bates hpmods function
  ##  documentation for which is below
  ##
  ##' Determine hierarchy-preserving models
  ##'
  ##' Determine the hierarchy-preserving models from a formula
  ##' @title Hierarchy-preserving models
      ##' @param formula as in \code{"\link{lm}"}
      ##' @param data as in \code{"\link{lm}"}
      ##' @param subset as in \code{"\link{lm}"}
      ##' @param weights as in \code{"\link{lm}"}
      ##' @param na.action as in \code{"\link{lm}"}
      ##' @param offset as in \code{"\link{lm}"}
      ##' @param \dots optional, additional arguments.  At present, none are used.
      ##' @return an object of class \code{"hpmods"}
      ##' @author Douglas Bates
      ##' @keywords models
      ##' @examples
      ##' set.seed(12321)
      ##' fr <- data.frame(y = rnorm(900), x1 = rnorm(900), x2 = rnorm(900),
      ##'                  x3 = gl(3,10,900), x4 = gl(10,3,900))
      ##' (hpm <- hpmods(y ~ (x1 + x2 + x3 + x4)^4, fr))
      ##' @export
      hpmodsetup <- function(formula, data, subset, weights, na.action, offset, ...) {
        cl   <- match.call()
        mf   <- match.call(expand.dots = TRUE)
        m    <- match(c("formula", "data", "subset", "weights", "na.action", 
                        "offset"), names(mf), 0L)
        mf   <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- as.name("model.frame")
        mf   <- eval(mf, sys.parent(n=2))
        atts <- attributes(terms(mf))
        inci <- crossprod(atts$factor) == atts$order # incidence matrix for terms
        mods <- array(FALSE, c(nrow(inci), 1))
        rownames(mods) <- rownames(inci)
        mods <- t(mods)
        rownames(mods) <- mods %*% 2^((seq_len(ncol(inci)))-1)
        res  <- list(call=cl, incidence=inci, models=mods,
                     frame=mf, X=model.matrix(terms(mf), mf))
        
        class(res) <- "hpmods"
        res
      }
      
      ##########################################################################
      #
      #   Used for the greedy search algorithm
      #   Calculates sum of squares for all the possible nested models  
      #       within the current set of terms.  It may recalculate the
      #       sum of squares for nested models.
      #   The terms must be provided in the sorted with the lowest ordered
      #       terms first to the highest ordered terms
      #
      critVal  <- function (thisTerms, hEnv) {
        
        terms <- thisTerms
        
        #  Create the list of columns needed for largest model
        #
        iCols   <- colAssign %in% c(0, terms)  # 0 is for intercept 
        
        # select the columns from the X matrix and calculate column effects 
        #
        ared   <- colAssign[iCols]   # column assignments for this model
        thisX  <- X[,iCols]
        QR  <- qr(thisX)
        SS.X  <- qr.qty(QR,respVec)[1:length(iCols)]
        SSCol  <- SS.X^2  
        
        #  Calculate the SS and df for the intercept if it is in the model
        #
        SSintcptCol  <- SSCol[which(ared==0)] # intercept is term 0 if it is included
        jSSMod  <- sum(SSintcptCol)   # include SS for intercept
        df  <- length(SSintcptCol)
        
        #  record selection method for collinearity reduced models
        if( QR$rank != ncol( as.matrix(thisX) ) ) {
          
          # remove dropped columns from ared list 
          #
          
          SSCol  <- SSCol[1:QR$rank]
          redColNames  <- rev(colnames(QR$qr))[1:(ncol(thisX) - QR$rank)]
          redColIdx  <- which( (colnames(X) %in% redColNames) )
          iCols[redColIdx]  <- FALSE
          ared  <- colAssign[iCols]   # column assignments for this model
          
          # Record SS for models for which collinearity dropped terms
          #
          redTerms  <- unique( colAssign[redColIdx] )
          
          # calculate search method
          nonColTerms  <- terms[ !(terms %in% redTerms) ]
          iSSMod  <- sum(SSCol[ared %in% nonColTerms]) + jSSMod
          idf  <- df + sum(ared %in% nonColTerms)
          
          if ( critCode == 0 ) {            # r2 or RSS
            iCritValue  <- SSTotUnCor - iSSMod
          }
          else if ( critCode == 1 ) {       # adjr2
            iCritValue  <- (SSTotUnCor - iSSMod) / max((dfTotal - idf),1)
            
          }
          else {                          # AIC or BIC
            iCritValue  <- log(SSTotUnCor - iSSMod) + penalty*idf
          }
          for ( i in length(redTerms):1 ) {
            iModIdx  <- sum(2^(terms-1))
            SSMod[iModIdx]  <- iCritValue          
            terms  <- terms[ (terms!=redTerms[i]) ]
          }
          
          # Record SS for models with remaining terms which are not free of droped terms
          # 
          drpTerms <- redTerms
          while ( sum( notFree[drpTerms,terms] ) > 0 ) {
            iModIdx  <- sum(2^(terms-1))
            
            # calculate search method
            #
            iSSMod  <- sum(SSCol[ared %in% terms]) + jSSMod
            idf  <- df + sum(ared %in% terms)
            if ( critCode == 0 ) {            # r2 or RSS
              iCritValue <- SSTotUnCor - iSSMod
            }
            else if (critCode == 1) {       # adjr2
              iCritValue <- (SSTotUnCor - iSSMod) / (dfTotal - idf)
              
            }
            else {                          # AIC or BIC
              iCritValue <- log(SSTotUnCor - iSSMod) + penalty*idf
              
            }
            SSMod[iModIdx]  <- iCritValue          
            terms  <- terms[ -length(terms) ]
            if ( length(terms) <= 0 ) {
              return()    
            }
          }
          
          # qr does not need to be called again.  Use remaining columns of qr
        }
        
        #  Calculate the SSModel for each model in the subset
        #
        termIds  <- terms
        nSet  <- length(terms)
        hEnv$numSets  <- hEnv$numSets + 1
        
        for ( j in 1:nSet ) {  # for each term of terms
          jTerm  <- termIds[j]
          
          # which columns to sum for this term
          jCols  <- ared %in% jTerm
          jColsSS  <- SSCol[which(jCols)]  
          
          # cumulative sum for terms so far 
          jSSMod  <- jSSMod + sum(jColsSS)
          iModIdx  <- sum(2^(terms[1:j]-1))
          SSMod[iModIdx]  <- jSSMod      
          df  <- df + sum(jCols)
          
          # calculate search method
          #
          if ( critCode == 0 ) {            # r2 or RSS
            jCritValue  <- SSTotUnCor - jSSMod
          }
          else if ( critCode == 1 ) {       # adjr2
            jCritValue  <- (SSTotUnCor - jSSMod) / max( (dfTotal - df), 1)
          }
          else {                          # AIC or BIC
            jCritValue  <- log(SSTotUnCor - jSSMod) + penalty*df
          }
          
          # save only if it a current best subset
          #
          SSModM  <- hEnv$SSModM
          termsModM  <- hEnv$termsModM
          if ( ( jCritValue < ( j.max  <- max(SSModM[,j]) ) ) &
                 !(iModIdx %in% termsModM[,j]) ) { # this model
            jRow  <- which(SSModM[,j]==j.max)[1]
            hEnv$jRow  <- jRow
            hEnv$j  <- j
            hEnv$df  <- df
            hEnv$jCritValue  <- jCritValue
            hEnv$terms  <- terms
            eval(parse(text = "SSModM[jRow,j]  <- jCritValue"), envir = hEnv)
            eval(parse(text = "termsModM[jRow,j]  <- sum(2^(terms[1:j]-1))"), envir = hEnv)
            eval(parse(text = "dfModM[jRow,j] <- df"), envir = hEnv)
          }
        }
      }
      
      ##########################################################################
      #
      #   FindBest recursively calls it self with the most significant term dropped and forced
      #   The most significant term is found by calculating the subset models of all free terms
      #
      
      findBest  <- function (thisTerms,   # ordered list of terms in model 
                             thisForced,  # list of terms forced into model
                             forced=FALSE, # Indicator for forced call
                             forcedFreeSS, # list of SS for free terms if force
                             terms,
                             notFree,
                             timeOut, 
                             mModDF,
                             hEnv) 
      {
        if ( length(thisTerms) <= 0 ) return()  # bad call if no predictors in the model
        
        if ( !forced ) {
          
          #   Find the list of terms that are free to be dropped
          #  the following line is not to be checked
          #
          
          thisNoForce  <- setdiff( thisTerms,thisForced )
          if ( length(thisTerms) > 1 ) {
            thisNotFree  <- notFree[,thisTerms]
            freeRows  <- ( rowSums(thisNotFree) == 0 )
            thisFree  <- terms[freeRows]
            free  <- thisFree[(thisFree %in% thisNoForce)]
          } else {
            free  <- thisTerms
          }
          
          if ( length(free)==0 ) return()  # no terms to drop or force
          
          #   Find the SSModel for each sub-model with a single 
          #   free term dropped
          #
          
          freeSS  <- c(0)[-1]
          maxSS  <- c(0)[-1]
          iMaxSS  <- -1
          if ( length(thisTerms)>1 ) {
            for ( i in free ) { 
              
              #   find the model id for the model with the dropped
              #   free term
              #
              iTerms  <- thisTerms[thisTerms != i]
              iModIdx  <- sum( 2^( iTerms-1 ) )
              
              #   if no SSMod for this model, calculate the SSMod for the subset which included it
              #
              
              if ( SSMod[iModIdx] < 0 ) critVal(iTerms, hEnv)
              
              #   add to the list of SSMods for the free terms
              #
              
              iMod  <- SSMod[iModIdx]
              names(iMod)  <- i
              freeSS  <- c(freeSS,iMod)
            }
          }  
        } else {       
          #  if forced, used the SSModels calculated at the
          #  above level
          freeSS  <- forcedFreeSS
        }
        
        thisTime  <- proc.time()
        if ( (thisTime - startTime)[3] > timeOut ) {  # over timeOut time?
          if ( outOfTime == FALSE ) {
            warning ("Processing time exceeded ", timeOut, 
                     " partial results have been returned")
          }
          outOfTime  <- TRUE
          return()
        }
        
        if ( length(freeSS)==0 ) return()  # no terms to drop or force
        
        #   find the best model of the set of dropped free terms
        #   which is the model which is least significant
        #          
        signfFree  <- names(which( freeSS == min(freeSS))[1] )[1] 
        
        #   If there are terms after dropping, continue dropping
        #   & forcing terms
        #
        newTerms  <- thisTerms[-which(thisTerms==signfFree)]
        
        if ( length(newTerms) >= 1 ) {   # check if there terms to drop or force
          
          #  Max criteria with this SS
          #
          
          if ( critCode == 0 ) {            # r2 or RSS
            mCritValue  <- freeSS[signfFree]
          }
          else if ( critCode == 1 ) {       # adjr2
            mCritValue  <- freeSS[signfFree] * max( (dfTotal - mModDF), 1 ) / ( dfTotal - 1 )
          }
          else {                          # AIC or BIC
            mCritValue  <- freeSS[signfFree] - penalty * mModDF
          }
          
          SSModM  <- hEnv$SSModM
          if ( mCritValue > max( SSModM[,1:length(newTerms)] )[1] ) return()
          
          # drop the term first
          findBest(newTerms, thisForced, terms=terms, notFree=notFree, timeOut=timeOut, mModDF=mModDF, hEnv = hEnv)    
          
          #    force the least significant term into the model
          newForced  <- c(thisForced, signfFree)
          newFreeSS  <- freeSS[-which(names(freeSS)==signfFree)]
          
          findBest(thisTerms, newForced, forced=TRUE,
                   forcedFreeSS=newFreeSS,  terms=terms, notFree=notFree, timeOut=timeOut, mModDF=mModDF, hEnv=hEnv)     
        }
      }
      
      ##################################################################################
      #
      #   hleaps code starts here
      #
      ##################################################################################
      
      #   initial parmeter handling
      #
      
      cl  <- match.call()
      mf  <- match.call(expand.dots = TRUE)
      m   <- match(c("formula", "data", "subset", "weights", "na.action", 
                     "offset"), names(mf), 0L)
      
      weightsOpt  <- FALSE
      weightsInp  <- NULL
      offsetOpt  <- FALSE
      offsetInp  <- NULL
      if(m[4] != 0){
        weightsOpt  <- TRUE
        weightsInp  <- list(...)$weights
      }
      
      if(m[6] != 0){
        offsetOpt  <- TRUE
        offsetInp  <- list(...)$offset
      }
      
      mf  <- mf[c(1L, m)]
      argList  <- as.list(sys.call())
      
      mf$drop.unused.levels  <- TRUE
      #   use lm to determine the number of terms
      #
      
      formula  <- as.formula(formula)  # coerce formula
      if ( missing(data) ) {
        lmOut  <- lm(formula)
      } else {
        lmOut  <- lm(formula, data=data)
      }
      
      terms  <- attr(lmOut$terms,"term.labels")
      matchTerms  <- vector()
      termList  <- terms
      numTerms  <- length(terms)
      
      #   Check that the number of terms is less than or equal to the maximum
      #   number of terms from available memory. 
      
      if( !is.numeric(mem) ) {
        stop ("mem must numeric - instead of "
                 , mem)
        return()
      }
      
      if( mem < 1 ) {
        stop ("must have at least 1 GB memory available")
        return()
      }
      
      maxTerms  <- floor(log(mem, 2)) + 26
      
      if ( numTerms > maxTerms ) {
        warning ("maximum number of terms must be less than ", maxTerms, " - instead of "
                 , numTerms)
        terms  <- terms[1:maxTerms]
        numTerms  <- maxTerms 
        newFormula  <- as.formula(paste(as.character(formula[2]),"~", 
                                        paste(terms, collapse= "+"),
                                        sep=""))
      } else {
        newFormula  <- formula
      }
      
      #   Check parameters and set defaults as needed
      #
      
      if ( is.null(minSize) ) minSize  <- 1 
      if ( is.null(maxSize) ) maxSize  <- numTerms
      if ( is.null(timeOut) ) timeOut  <- Inf  
      
      if ( weightsOpt & method == "RSS" ){
        warning("calculated RSS will not be on original scale. Calculated RSS will be on weighted scale.")
      }
      if( numTerms == 1 ) {
        stop ("number of terms must be greater than 1")
        return();
      }
      if( timeOut < 0 ){
        warning ("timeOut must be greater than 0 - instead of ",
                 timeOut, " - defaulting timeOut to Inf")
        timeOut  <- Inf;
      }
      if( !is.numeric(timeOut) ) {
        warning ("timeOut must be numeric - instead of ",
                 timeOut, " - defaulting timeOut to Inf")
        timeOut  <- Inf;
      }
      if ( minSize >= numTerms ) { 
        warning ("minSize must be less than number of terms - instead of "
                 , minSize)
        minSize  <- 1
      }
      if ( minSize > numTerms ) {
        warning ("minSize cannot be greater than the number of possible terms - instead of "
                 , minSize, " - defaulting minSize to 1")
        minSize  <- 1
      }
      if ( minSize < 1 ) { 
        warning ("minSize must be at least 1 - instead of "
                 , minSize)
        minSize  <- 1
      }
      if ( !is.numeric(minSize) ) {
        warning ("minSize must be at least 1 - instead of "
                 , minSize, " - defaulting minSize to 1")
        minSize  <- 1
      }
      if( maxSize == 0 ) {
        warning ("maxSize must be at least 1 - instead of "
                , maxSize)
        max.size  <- numTerms
      }
      if ( maxSize < minSize ) { 
        warning ("maxSize must be at least as large as minSize - instead of "
                 , maxSize)
        maxSize  <- numTerms
      }
      if ( !is.numeric(maxSize) ) {
        warning ("maxSize must be at least 1 - instead of "
                 , maxSize, " - defaulting maxSize to number of terms")
        maxSize  <- numTerms
      }
      if ( !(method %in% c("adjr2","r2","RSS","AIC","BIC")) ) {
        warning ("method must me adjr2, r2, RSS, AIC, or BIC - instead of "
                 , method, " - method defaulting to adjr2")   
        method  <- "adjr2"
      }
      if ( nbest < 1 ) {
        warning ("nbest must be greater than 0 - instead of "
                 , nbest)   
        nbest  <- 3
      }    
      if ( !is.numeric(nbest) ) {
        warning ("nbest must be numeric - instead of "
                 , nbest, " - defaulting to 3")   
        nbest  <- 3
      }
      if ( !( is.logical(altOut) ) ) {
        warning ("altOut must be logical - instead of "
                 , altOut)   
        altOut  <- FALSE
      }
      if( is.null(na.action) ){
        na.action = NULL
      }
      
      #  Create list of models and grouped subsets of models using 
      #  Doug Bates's code and get the needed data
      #
      
      newFormula  <- as.formula(newFormula)  # coerce formula
      
      # check weights length
      if(is.null(weightsInp) == FALSE) {
        if(length(weightsInp) != nrow(model.frame(newFormula))) {
          stop ("weights vector must equal length of response and covariates")
          return()
        }
        if(any(is.na(weightsInp))) {
          stop ("weights input cannot contain missing values")
          return()
        }
        if(sum(weightsInp < 0) != 0) {
          stop ("weights input cannot contain negative values")
          return()
        }
      }
      
      # check offset length
      if(is.null(offsetInp) == FALSE) {
        if(length(offsetInp) != nrow(model.frame(newFormula))) {
          stop ("offset vector must equal length of response and covariates")
          return()
        }
      }
      
      cl   <- match.call()
      mf   <- match.call(expand.dots = TRUE)
      m    <- match(c("formula", "data", "subset", "weights", "na.action", 
                      "offset"), names(mf), 0L)
      hp   <- mf[c(1L, m)]
      hp$drop.unused.levels <- TRUE
      hp[[1L]] <- as.name("hpmodsetup")
      hp   <- eval(hp)
      mods  <- hp$models
      X  <- hp$X
      colAssign <-  attr(X,"assign")
      
      if( weightsOpt ) {
        zeroWeights  <- which(weightsInp == 0)
        rows  <- 1:length(weightsInp)
        X  <- X[setdiff(rows, zeroWeights),]
        hp$frame  <- hp$frame[setdiff(rows, zeroWeights),]
        weightsInp  <- weightsInp[setdiff(rows, zeroWeights)]
      }
      
      dfTotal  <- nrow(X)
      
      incidence  <- hp$incidence
      termLabels  <- attr(attr(hp$frame,"terms"),"term.labels")
      terms  <- seq(1,length(termLabels))
      termOrder  <- attr(attr(hp$frame,"terms"),"order")
      names(termOrder)  <- terms
      respVec  <- as.numeric(model.response(hp$frame))
      
      # apply offset
      if("(offset)" %in% colnames(hp$frame)){
        respVec  <- respVec - hp$frame[,"(offset)"]
      }
      
      # apply weights
      if("(weights)" %in% colnames(hp$frame)) {
        X  <- X * sqrt(hp$frame[,"(weights)"])
        respVec  <- respVec * sqrt(hp$frame[,"(weights)"])
      }
      
      SSTotUnCor  <- sum(as.numeric(respVec)^2)
      
      #   Decomposition of full model
      #
      
      hEnv  <- new.env()
      QR  <- qr(X)
      effList  <- qr.qty(QR,respVec)
      #colAssign  <- attr(X,"assign")
      mModDF  <- QR$rank
      
      # prepare variables for best subsets method
      #
      
      if ( ( method %in% c("r2","RSS") ) ) {
        critCode <- 0
      }
      else if ( (method %in% c("adjr2"))) {
        critCode <- 1
      }
      else {   # AIC or BIC
        critCode <- 2
        if ( method =="AIC" ) {
          penalty <- 2/dfTotal
        }
        else {
          penalty <- log(dfTotal)/dfTotal
        }
      }
      
      #   table to record the Sum of Squares for each model
      #
      maxMods  <- 2^(numTerms)
      SSMod  <- rep(-1,maxMods)
      
      #   Matrices to record the best Sum of Squares models
      #
      hEnv$SSModM  <- matrix(rep(Inf,nbest*numTerms), nrow=nbest, ncol=numTerms)
      hEnv$termsModM  <- matrix(rep(0,nbest*numTerms), nrow=nbest, ncol=numTerms)
      hEnv$dfModM  <- matrix(rep(0,nbest*numTerms), nrow=nbest, ncol=numTerms) 
      
      #   search models by subsets for best models
      #
      startTime  <- proc.time()
      outOfTime  <- FALSE
      
      #   incident table used to find free terms
      #
      maxOrder  <- max(termOrder)    
      notFree  <- incidence
      rownames(notFree)  <- terms
      colnames(notFree)  <- terms
      for ( i in terms ) {
        iOrder  <- termOrder[i]
        notFree[(termOrder %in% iOrder:maxOrder),i]  <-0
      }
      
      forcedTerms  <- c("")[-1]
      neededTerms  <- terms
      hEnv$numSets  <- 0
      critVal(neededTerms, hEnv)   # calc SS for all model
      
      # run greedy search
      findBest(neededTerms, forcedTerms, terms=terms, notFree=notFree, timeOut=timeOut, mModDF=mModDF, hEnv=hEnv)
      
      output  <- format.output(minSize, maxSize, nbest, altOut, method, colAssign, dfTotal, mods,
                               numTerms, respVec, startTime, SSTotUnCor, hEnv, weightsOpt, weightsInp, offsetOpt, offsetInp)
      
      return(output)
}


###########################################################################
#
#   format the compressed which matrix of subsets
#

format.subsets  <- function (sdf) {
  numRow  <- nrow(sdf)
  numCol  <- ncol(sdf)
  numDCol  <- numCol %/% 5
  if ( ( numCol %% 5 ) > 0 ) numDCol  <- numDCol + 1
  
  sd  <- data.frame(c(1:numRow))[,-1]  # create an empty data frame
  row  <- 1
  while ( row  <= numRow ) {
    col  <- 1
    dCol  <- 1
    while ( col  <= numCol ) {
      if ( ( col + 4 )  < numCol ) {
        endCol  <- col + 4
      } else {
        endCol  <- numCol
      }
      sd[row,dCol]  <- paste(paste(ifelse(sdf[row,col:endCol]==TRUE,"T",".")), collapse= "")
      col  <- col + 5
      dCol  <- dCol + 1
    }
    row  <- row + 1
  }
  
  # Create column names
  #
  col  <- 1
  dCol  <- 1
  sdNames  <- c(1)
  while ( col <= numCol ) {
    if ( ( col + 4 )  < numCol ) {
      endCol  <- col + 4
    } else {
      endCol  <- numCol
    }
    this.name  <- c( seq( col + 1 , endCol ) )
    this.index  <- this.name[1] %/% 10
    this.name  <- this.name %% 10
    sdNames[dCol]  <- paste(col,endCol,sep="-")
    col  <- col + 5
    dCol  <- dCol + 1
  }
  names(sd)  <- sdNames
  sd
}


###########################################################################
#
#
#   format the output
#
format.output  <- function(minSize, maxSize, nbest, altOut, method, colAssign, dfTotal, mods, numTerms, respVec,
                           startTime, SSTotUnCor, hEnv, weightsOpt, weights = NULL, offsetOpt, offsets = NULL) {
  
  # Get best saved models information
  #
  dfModM  <- hEnv$dfModM
  SSModM  <- hEnv$SSModM
  termsModM  <- hEnv$termsModM
  notEmpty  <- as.vector(SSModM < Inf)      # identify the models which are non empty
  size  <- rep(1:numTerms, each = nbest)[notEmpty]
  dfMod  <- as.vector(dfModM)[notEmpty]
  critValues  <- as.vector(SSModM)[notEmpty]
  modId  <- as.vector(termsModM)[notEmpty]
  
  #   Code to convert index to terms (binary) 
  #
  decId  <- modId
  modTerms  <- matrix(ncol=numTerms, nrow=length(decId))
  for ( i in numTerms:1 ) {    
    modTerms[,i]  <- decId >= 2^( i-1 )
    decId  <- decId - ifelse( modTerms[,i],2^(i-1) , 0 )
  }
  termsMatrix  <- modTerms
  
  # code to correct the critical values
  #
  
  if ( 0 %in% colAssign ) {                     # is there an intercept in the model
    if(!weightsOpt) {
      if(!offsetOpt) {
        SSTotal  <- sum( ( as.numeric(respVec) - mean(as.numeric(respVec) ) ) ^2 )  # correct for mean
      }
      else{
        SSTotal  <- sum( ( as.numeric(respVec + offsets) - mean(as.numeric(respVec + offsets) ) ) ^2 ) 
      }
    }
    else{
      inpWeights  <- weights # duplicate weights for debug
      sqWeights  <- sqrt(weights) # root weights, this is wts
      null <- sum( (sqWeights * respVec / sum(inpWeights) ) )  # null model
      wNull <- sqWeights * null # multiply null model by weights
      mssNull <- sum( wNull^2 ) # sum weighted null
      SSTotal  <- sum(respVec^2) - mssNull # subtract weighted null and sq
      #error()
      
    }
    
    dfTotal2  <- length(respVec) - 1
  }
  else {
    SSTotal  <- SSTotUnCor                      # no corrections for mean
    dfTotal2  <- length(respVec)
  }
  if ( method == "r2" ) {
    critValueCor  <- 1 - critValues / SSTotal
  } 
  else if ( method == "adjr2" ) {
    critValueCor  <- 1 - critValues * (dfTotal2 / SSTotal)
  }
  else if ( method == "AIC" ) {
    if(weightsOpt) {
      critValueCor  <- critValues + log(1/dfTotal)  + log( prod(inpWeights^(-1)) )
    }
    else {
      critValueCor  <- critValues + log(1/dfTotal)  
    }
  }
  else if ( method == "BIC" ) {
    if(weightsOpt) {
      critValueCor  <- critValues + log(1/dfTotal)  +log( prod(inpWeights^(-1)) )
    }
    critValueCor  <- critValues  + log(1/dfTotal) 
  }
  else if ( method == "RSS" ) {
    critValueCor  <- critValues
  }
  ModCriteria <- critValueCor
  
  # order the models in the return object by number of terms and then selection method 
  # 
  
  o  <- sort(critValues)
  modOrder  <- match(critValues, o)
  uSizes  <- unique(size)
  
  for( i in 1:length(uSizes) ){
    thisSize  <- which(size == uSizes[i])
    thisSSMod  <- critValueCor[thisSize]
    thisDecSS  <- order(thisSSMod,decreasing=TRUE)
    
    orderModOrder  <- modOrder[thisSize]
    modOrder[thisSize]  <- orderModOrder[thisDecSS]
    
    orderdfMod  <- dfMod[thisSize]
    dfMod[thisSize]  <- orderdfMod[thisDecSS]
    
    orderCrit  <- critValueCor[thisSize]
    ModCriteria[thisSize]  <- orderCrit[thisDecSS]
    
    orderThis  <- termsMatrix[thisSize,]
    if( length(thisDecSS) > 1 ) {
      termsMatrix[thisSize,]  <- orderThis[thisDecSS,]
    }
  }
  
  if ( !altOut ) {  # use normal output format
    colnames(termsMatrix)  <- as.character(seq(1:ncol(termsMatrix)))
    rownames(termsMatrix)  <- as.character(size)
    label <- colnames(mods)
    
    if ( 0 %in% colAssign ) {
      label <- c("(Intercept)",label)
    }
    
    BSS  <- list(which=termsMatrix, label=label, 
                 size=size, method=ModCriteria, df = dfTotal2)   
    names(BSS)[names(BSS)=="method"]  <- method
    
  } else {   # alternate formated output
    
    #   compress the which matrix
    #
    comp.mod  <- format.subsets(termsMatrix)
    
    # format execution info
    #
    totalTime  <- proc.time()[3] - startTime[3]
    
    time.usage  <- paste("Total number of model subsets examined was ",
                         hEnv$numSets,
                         " Total run time was ", 
                         format(totalTime,trim = FALSE,digits=7), 
                         " seconds.")
    modelInfo  <- cbind(size, df=dfMod, order=modOrder, method=ModCriteria, 
                        comp.mod)   
    colnames(modelInfo)[colnames(modelInfo)=="method"]  <- method
    
    BSS  <- list(modelInfo=modelInfo, label=colnames(mods),
                 executionInfo=time.usage)
  }
  return(BSS)
}

###########################################################################
#
# approximate hleaps time to run
#
#' @export
hleapsApproxTime = function(n=1000, p=10, int=TRUE, sigI=TRUE, sigC=TRUE){
  x1  <- rnorm(10000)
  x2  <- rnorm(10000)
  x3  <- rnorm(10000)
  x4  <- rnorm(10000)
  y  <- rnorm(10000)
  timeStart = proc.time()
  hleapsOut  <- hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE)
  timeEnd = proc.time()
  
  timeElapsed = timeEnd[3] - timeStart[3]
  timeElapsed = as.numeric(timeElapsed)
  time = (7.135*(10^-7) * n + .01872 * p^2 + .0824 * p + -1.492* int + -2.252 * sigI + 1.563 * sigC  - .28781) * (timeElapsed / 7)
  time = exp(1)^time

  days = time/(24*60*60)
  daysF = floor(days)
  
  hours = (time/(60*60))%%24
  hoursF = floor(hours)
  
  minutes = (time/60)%%60
  minutesF = floor(minutes)
  
  seconds = time%%60
  secondsF = floor(seconds)
 
  message ("The approximate runtime for the hleaps algorithm is ", daysF, " days, ", hoursF, " hours, ", minutesF, " minutes, ",
           "and ", secondsF, " seconds.") 
  
  timeList = list()
  timeList$days = daysF
  timeList$hours = hoursF
  timeList$minutes = minutesF
  timeList$seconds = secondsF
  return(timeList)
}