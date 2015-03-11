\name{format.subsets}
\alias{format.subsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
format subsets
}
\description{
Formats the selected subsets for alternative output.
}
\usage{
format.subsets(sdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{
%%     ~~Describe \code{sdf} here~~
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
function (sdf) 
{
    numRow <- nrow(sdf)
    numCol <- ncol(sdf)
    numDCol <- numCol\%/\%5
    if ((numCol\%\%5) > 0) 
        numDCol <- numDCol + 1
    sd <- data.frame(c(1:numRow))[, -1]
    row <- 1
    while (row <= numRow) {
        col <- 1
        dCol <- 1
        while (col <= numCol) {
            if ((col + 4) < numCol) {
                endCol <- col + 4
            }
            else {
                endCol <- numCol
            }
            sd[row, dCol] <- paste(paste(ifelse(sdf[row, col:endCol] == 
                TRUE, "T", ".")), collapse = "")
            col <- col + 5
            dCol <- dCol + 1
        }
        row <- row + 1
    }
    col <- 1
    dCol <- 1
    sdNames <- c(1)
    while (col <= numCol) {
        if ((col + 4) < numCol) {
            endCol <- col + 4
        }
        else {
            endCol <- numCol
        }
        this.name <- c(seq(col + 1, endCol))
        this.index <- this.name[1]\%/\%10
        this.name <- this.name\%\%10
        sdNames[dCol] <- paste(col, endCol, sep = "-")
        col <- col + 5
        dCol <- dCol + 1
    }
    names(sd) <- sdNames
    sd
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
