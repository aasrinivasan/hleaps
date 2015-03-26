\name{hleapsApproxTime}
\alias{hleapsApproxTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate hleaps runtime
}
\description{
This method estimates the time needed to execute the hleaps algorithm
for a given model size.
}
\usage{
hleapsApproxTime(n=1000, p=10, int=TRUE, sigI=TRUE, sigC=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
an integer - the number of data observations in the data set
}
  \item{p}{
an integer - the number of variables to be used in subset selection.  This count is only for the main effects.
}
  \item{int}{
a boolean - are interactions included in model?
It TRUE, the estimate is for a model with all second order
interaction terms.
This provides an upper bound on the algorithm time.
}
  \item{sigI}{
a boolean - are interactions expected to be significant?
If TRUE, the estimate is for a model with all the
second order interaction terms being significant.
This provides an upper bound on the algorithm time.
False assumes none of the interaction terms are significant.
}
  \item{sigC}{
a boolean - are covariates expected to be significant?
If TRUE, the estimate is for a model with all the
main effects being significant.
This provides an upper bound on the algorithm time.
False assumes none of the main effects are significant.
}

}
\details{ 
The function prints out the estimated run time for the hleaps algorithm.
This function times how long it takes to run hleaps on a calibrated data set
on the currently running computer.
This calibrated time is then scaled to the model using a formula
provided below.

logTime = (7.135*(10^-7) * n + .01872 * p^2 + .0824 * p + -1.492* int + -2.252 * sigI + 1.563 * sigC  - .28781) * (timeOnUserMachine)/(timeOnEstimationRuns)

This estimate is attained from running  from data collected through timing estimation runs. 
The formula was generated from a regression on the results of
multiple timing runs of hleaps with varying values for 
n, p, int, sigI and sig C. 
The linear regression was done on the log(times) collected for each run. 
A log normal correction was applied to the model.

This function will take at least 10 seconds to run in order to ensure a stable result.
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
