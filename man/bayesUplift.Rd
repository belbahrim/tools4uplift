\name{bayesUplift}
\alias{bayesUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Direct Uplift Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bayesUplift(data, outcome, treatment, predictors, folds=10, prior=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame containing the response, the treatment and the variables to include in the model.
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{treatment}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{predictors}{
vector of names or indexes representing the explanatory variables to include in the model.
}
  \item{folds}{
number of folds - default is 10. Although folds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is \code{folds=3}.
}
  \item{prior}{
indication on how you want to evaluate the propensity - default is \code{TRUE}. This means only control group observations will be used. If \code{FALSE}, all observations will be used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
An object of class \code{bayesUplift}, which is a list with the following components:
  \item{data}{an augmented data frame containing fitted uplift.}
  \item{outcome}{a name of the binary response.}
  \item{treatment}{a name of the binary treatment assignment variable.}
  \item{predictors}{a vector of names or indexes representing the explanatory variables.}
  \item{cv.fit}{an object of class \code{cv.glmnet} for the propensity model.}
  \item{cv.fit0}{an object of class \code{cv.glmnet} for the first conditional model.}
  \item{cv.fit1}{an object of class \code{cv.glmnet} for the second conditional model.}
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

library(tools4uplift)
#do not output the functions plots
scene <- generateScene(plotit=FALSE)

#generate 1000 observations with 20 predictors
data <- generateUplift(scene, 1000, 20)

#fit a Bayes Uplift Model
fit1 <- bayesUplift(data, outcome="y", treatment="treat", predictors=colnames(data[,3:22]), folds=10, prior=FALSE)

#save fitted uplift in a data frame
data_fit1 <- fit1[[1]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
