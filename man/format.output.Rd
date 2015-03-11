\name{format.output}
\alias{format.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
format output
}
\description{
Formats hleaps output and calculates returned values based on method selection.
}
\usage{
format.output(minSize, maxSize, nbest, altOut, method, colAssign, dfTotal, mods, numTerms, respVec, startTime, SSTotUnCor, hEnv, weightsOpt, weights = NULL, offsetOpt, offsets = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minSize}{
%%     ~~Describe \code{minSize} here~~
}
  \item{maxSize}{
%%     ~~Describe \code{maxSize} here~~
}
  \item{nbest}{
%%     ~~Describe \code{nbest} here~~
}
  \item{altOut}{
%%     ~~Describe \code{altOut} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{colAssign}{
%%     ~~Describe \code{colAssign} here~~
}
  \item{dfTotal}{
%%     ~~Describe \code{dfTotal} here~~
}
  \item{mods}{
%%     ~~Describe \code{mods} here~~
}
  \item{numTerms}{
%%     ~~Describe \code{numTerms} here~~
}
  \item{respVec}{
%%     ~~Describe \code{respVec} here~~
}
  \item{startTime}{
%%     ~~Describe \code{startTime} here~~
}
  \item{SSTotUnCor}{
%%     ~~Describe \code{SSTotUnCor} here~~
}
  \item{hEnv}{
%%     ~~Describe \code{hEnv} here~~
}
  \item{weights}{
%%     ~~Describe \code{hEnv} here~~
}
  \item{weightsOpt}{
%%     ~~Describe \code{hEnv} here~~
}
  \item{offsetOpt}{
  
}
  \item{offsets}{
  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (minSize, maxSize, nbest, altOut, method, colAssign, 
    dfTotal, mods, numTerms, respVec, startTime, SSTotUnCor, 
    hEnv) 
{
    dfModM <- hEnv$dfModM
    SSModM <- hEnv$SSModM
    termsModM <- hEnv$termsModM
    notEmpty <- as.vector(SSModM < Inf)
    size <- rep(1:numTerms, each = nbest)[notEmpty]
    dfMod <- as.vector(dfModM)[notEmpty]
    critValues <- as.vector(SSModM)[notEmpty]
    modId <- as.vector(termsModM)[notEmpty]
    decId <- modId
    modTerms <- matrix(ncol = numTerms, nrow = length(decId))
    for (i in numTerms:1) {
        modTerms[, i] <- decId >= 2^(i - 1)
        decId <- decId - ifelse(modTerms[, i], 2^(i - 1), 0)
    }
    termsMatrix <- modTerms
    if (0 \%in\% colAssign) {
        SSTotal <- sum((as.numeric(respVec) - mean(as.numeric(respVec)))^2)
        dfTotal2 <- length(respVec) - 1
    }
    else {
        SSTotal <- SSTotUnCor
        dfTotal2 <- length(respVec)
    }
    if (method == "r2") {
        critValueCor <- 1 - critValues/SSTotal
    }
    else if (method == "adjr2") {
        critValueCor <- 1 - critValues * (dfTotal2/SSTotal)
    }
    else if (method == "AIC") {
        critValueCor <- critValues + log(1/dfTotal)
    }
    else if (method == "BIC") {
        critValueCor <- critValues + log(1/dfTotal)
    }
    else if (method == "RSS") {
        critValueCor <- critValues
    }
    ModCriteria <- critValueCor
    o <- sort(critValues)
    modOrder <- match(critValues, o)
    uSizes <- unique(size)
    for (i in 1:length(uSizes)) {
        thisSize <- which(size == uSizes[i])
        thisSSMod <- critValueCor[thisSize]
        thisDecSS <- order(thisSSMod, decreasing = TRUE)
        orderModOrder <- modOrder[thisSize]
        modOrder[thisSize] <- orderModOrder[thisDecSS]
        orderdfMod <- dfMod[thisSize]
        dfMod[thisSize] <- orderdfMod[thisDecSS]
        orderCrit <- critValueCor[thisSize]
        ModCriteria[thisSize] <- orderCrit[thisDecSS]
        orderThis <- termsMatrix[thisSize, ]
        if (length(thisDecSS) > 1) {
            termsMatrix[thisSize, ] <- orderThis[thisDecSS, ]
        }
    }
    if (!altOut) {
        colnames(termsMatrix) <- as.character(seq(1:ncol(termsMatrix)))
        rownames(termsMatrix) <- as.character(size)
        label <- colnames(mods)
        if (0 \%in\% colAssign) {
            label <- c("(Intercept)", label)
        }
        BSS <- list(which = termsMatrix, label = label, size = size, 
            method = ModCriteria)
        names(BSS)[names(BSS) == "method"] <- method
    }
    else {
        comp.mod <- format.subsets(termsMatrix)
        totalTime <- proc.time()[3] - startTime[3]
        time.usage <- paste("Total number of model subsets examined was ", 
            hEnv$numSets, " Total run time was ", format(totalTime, 
                trim = FALSE, digits = 7), " seconds.")
        modelInfo <- cbind(size, df = dfMod, order = modOrder, 
            method = ModCriteria, comp.mod)
        colnames(modelInfo)[colnames(modelInfo) == "method"] <- method
        BSS <- list(modelInfo = modelInfo, label = colnames(mods), 
            executionInfo = time.usage)
    }
    return(BSS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
