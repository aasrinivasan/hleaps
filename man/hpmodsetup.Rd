\name{hpmodsetup}
\alias{hpmodsetup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hpmodsetup
}
\description{
Creation of full model design matrix and incidence matrix
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
Douglas Bates
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
function (formula, data, subset, weights, na.action, offset, 
    ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    atts <- attributes(terms(mf))
    inci <- crossprod(atts$factor) == atts$order
    mods <- array(FALSE, c(nrow(inci), 1))
    rownames(mods) <- rownames(inci)
    mods <- t(mods)
    rownames(mods) <- mods \%*\% 2^((seq_len(ncol(inci))) - 1)
    res <- list(call = cl, incidence = inci, models = mods, frame = mf, 
        X = model.matrix(terms(mf), mf))
    class(res) <- "hpmods"
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
