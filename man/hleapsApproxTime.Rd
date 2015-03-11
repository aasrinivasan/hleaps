\name{hleapsApproxTime}
\alias{hleapsApproxTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate hleaps runtime
}
\description{
This method estimates the number of seconds needed to execute the hleaps algorithm for a given set of parameters.
}
\usage{
hleapsApproxTime(n, p, int, sigI, sigC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
an integer - the number of data points in the data set
}
  \item{p}{
an integer - the number of variables to be used in subset selection
}
  \item{int}{
a boolean - are interactions included in model?
}
  \item{sigI}{
a boolean - are interactions expected to be significant?
}
  \item{sigC}{
a boolean - are covariates expected to be significant?
}

}
\details{ The function prints out the estimated run time for the hleaps algorithm, attained from data collected through timing estimation runs. Multiple timing runs were preformed with varying values for n, p, int, sigI and sig C. A linear regression was conducted on the log(times) collected for each run. A log normal correction was applied to the model. The model found is:

logTime = 7.135*(10^-7) * n + .01872 * p^2 + .0824 * p + -1.492* int + -2.252 * sigI + 1.563 * sigC  - .28781
}
\value{
An object of class list.

\itemize{
\item{days} {
- an integer. The approximate number of days needed for hleaps algorithm
}
\item{hours} {
- an integer. The approximate number of hours needed for hleaps algorithm
}
\item{minutes} {
- an integer. The approximate number of minutes needed for hleaps algorithm
}
\item{seconds} {
- an integer. The approximate number of seconds needed for hleaps algorithm
}
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Arun Srinivasan <aasrinivasan@wisc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
hleapsApproxTime(1000, 18, TRUE, FALSE, FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
