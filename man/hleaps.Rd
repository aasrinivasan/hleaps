\name{hleaps}
\alias{hleaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchy supporting leaps
}
\description{
This method performs an exhaustive, recursive, search though all sets of hierarchy preserving ordinary least squares models, to identify the best subsets of provided terms. For each model the search criteria specified by the user is calculated. The calculated search criteria for each model are then used to determine the best subsets of provided terms in the model.
}
\usage{
hleaps(formula, data, method = "adjr2", nbest = 3, minSize = NULL, maxSize = NULL, timeOut = NULL, altOut = FALSE, mem = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a character scalar or \code{\link{formula}} formula object - a two-sided linear formula object describing the model, with the response on the left of a ~ operator and the terms on the right, using any of the formula symbols for the lm function.
}
  \item{data}{
an object of class data.frame - an optional data frame containing the variables named in formula. By default the variables are taken from the environment from which hleaps is called.
}
  \item{method}{
a character scalar - criteria to be used for best model selection. The valid set of model selection criteria is \dQuote{RSS}, \dQuote{r2}, \dQuote{adjr2}, \dQuote{AIC}, and \dQuote{BIC}.
}
  \item{nbest}{
an integer scalar - number of best models per size
}
  \item{minSize}{
an integer scalar - smallest size of model reported, where size is measured in number of model terms. NULL defaults to size of 1.
}
  \item{maxSize}{
integer scalar - largest size of model reported, where size is measured in number of model terms. NULL defaults to the number of terms specified in the formula parameter.}
\item{timeOut}{
a numeric scalar - number of seconds to search for models. NULL defaults to Inf.
}
  \item{altOut}{
a logical scalar - TRUE indicates a request for the alternate returned object.
}
  \item{mem}{
a numeric scalar - Should identify the amount of memory available for the the hleap function.
}
  \item{\dots}{
other potential arguments. These are arguments which would be passed to lm such as weights, subsets, etc.
}
}
\details{The search algorithm uses a binary tree approach. The tree branches based on models which contain or do not contain a particular model term. The search is done first on the branch with the term dropped from the model. The selection of which term to drop is done via the method provided. 

Thus the algorithm attempts to find the best models, as determined by method, first. The best subsets are grouped by number of terms. The degrees of freedom for models with the same number of terms may vary.

If the algorithm is stopped by user run time limit, the algorithm may not have identified the best subsets for each number of terms requested.
}
\value{
An object of class list.

For the altOut = FALSE
\itemize{
\item{which} {
- a logical matrix. Each row can be used to select the columns of x in the respective model. The logical matrix includes
}
\item{size} {
- an integer vector. The number of variables, including intercept if any, in the model
}
\item{\dQuote{RSS}, \dQuote{r2}, \dQuote{adjr2}, \dQuote{AIC}, or \dQuote{BIC}} {
- a numeric vector lists the value of the chosen model selection statistic for each model
}
\item{label} {
a character vector - names for the columns, terms, of x
}
}
For  altOut = TRUE
\itemize{
\item{modelInfo} {
- a data.frame each row represents an individual model. First column indicates the number of terms in the select model. Second column represents the number of degrees of freedom in model. Third column states the individual ranking of the model. Fourth column is the value of the selected search criteria for model. Remaining columns can be used to select the terms in model, based on ordering label.
}
\item{label} {
- a character vector a vector containing term names in matching order to those found in modelInfo
}
\item{executionInfo} {
- a character vector detailing the number of model subsets examined and the amount of time the algorithm ran.
}
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark Banghart <mbanghart@ssc.wisc.edu>, Douglas Bates, Arun Srinivasan <aasrinivasan@wisc.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

set.seed(21564)
x1  <- rnorm(100)
x2  <- rnorm(100)
x3  <- rnorm(100)
x4  <- rnorm(100)
y  <- rnorm(100)

hleapsOut  <- hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
