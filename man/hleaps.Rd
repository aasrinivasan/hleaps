\name{hleaps}
\alias{hleaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchy supporting leaps
}
\description{
This method performs a search of all sets of hierarchy preserving models. 
The search criteria for best fit can be R^2, adjusted R^2, AIC, and BIC.
the best models for each size of models are reported.
}
\usage{
hleaps(formula, data, method = "adjr2", nbest = 3, 
        minSize = NULL, maxSize = NULL, timeOut = NULL, 
        altOut = FALSE, mem = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a fomula <this should be alink to the formula page> object or a character scalar which can be coereced to formula object.
}
  \item{data}{
an object of class data.frame - an optional data frame containing the variables named in formula. By default the variables are taken from the environment from which hleaps is called.
}
  \item{method}{
a character scalar - criteria to be used for best model selection. The valid set of model selection criteria is \dQuote{RSS}, \dQuote{r2}, \dQuote{adjr2}, \dQuote{AIC}, and \dQuote{BIC}.
}
  \item{nbest}{
an integer scalar - number of best models per size model where size is the
number of model terms.
}
  \item{minSize}{
an integer scalar - smallest size of model to be reported. 
NULL defaults to size of 1.
}
  \item{maxSize}{
integer scalar - largest size of model to be reported. 
NULL defaults to the number of terms specified in the formula parameter.
}
\item{timeOut}{
a numeric scalar - number of seconds to after which the search for models is terminated. 
NULL defaults to Inf.
}
  \item{altOut}{
a logical scalar - TRUE indicates a request for the alternate returned object.
}
  \item{mem}{
a numeric scalar - identifies the amount of memory, measured in gig, available for the the hleap function.
}
  \item{\dots}{
other potential arguments. These are arguments which would be passed to lm such as weights, subsets, etc.
}
}
\details{
The search algorithm uses a binary tree approach. 
The tree branches based on models which contain or do not contain a particular model term. 
The algorithm starts with all provided terms in the model.
The dropped term branch is searched first. 
The term with the smallest reduction in the
search method criteria is dropped. 
Thus the algorithm attempts to search the best models first, 
as determined by the search method critera. 

The search algorithm may not need to build all possible
models.
The algorithm determines if all further reductions of the
current model will not produce a best subset model.
If so, the remainder of the branch will be not be searched.
When the algorithm fits a model, the search criteria
for that model and the set of model nested within
this model by dropping the right most terms dropped
are also calculated.

The best subsets are grouped by number of terms. 
The degrees of freedom for models with the same number of terms may vary.

If the algorithm is stopped by user provided search time limit, 
the algorithm may not have identified the best subsets for each number of terms requested.

The amount of memory needed to run the algorithm can be large.
If the algorithm needs more memory than what the mem
parameter states is available, an error message is produced 
and the function stops.



offsets are assumed to be in the column space of the X matrix by lm and summary.lm.
The use of offsets which remove a column of the X matrix in algorithm which 
searches all models requires care.
hleaps will consider models which will result in the offsets is not in the column
space of the X matrix.
The calculation of R^2 or adjusted R^2 from hleaps will not match the calculation
from summary.lm in this case.
If one has an application of offset which is not in the column space of X,
hleaps will correctly prefrom the search based on R^2 or adjusted R^2.
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
\item{df} {
- an integer vector. The number of degrees of freedom for the model
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
\author{
Mark Banghart <mbanghart@ssc.wisc.edu>, Douglas Bates, Arun Srinivasan <ara.arun@gmail.com>
}

\examples{

set.seed(21564)
x1  <- rnorm(100)
x2  <- rnorm(100)
x3  <- rnorm(100)
x4  <- rnorm(100)
y  <- rnorm(100)

# Example method with alternative output
hleapsOutAlt  <- hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE)

# Example method with standard leaps output
hleapsOutStandard  <- hleaps(y~(x1+x2+x3+x4)^2, altOut = FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
