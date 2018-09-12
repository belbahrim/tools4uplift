\name{bestFeatures}
\alias{bestFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature Selection for the Interaction Uplift Model Estimator
}
\description{
Penalized logistic regression in order to select the features that maximize the Qini coefficient. 
}
\usage{
bestFeatures(data, outcome, treat, predictors, nb_lambda=100, nb_group=10, value=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame containing the response, the treatment and the predictors.
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{predictors}{
a vector of names representing the predictors to consider in the model.
}
  \item{nb_lambda}{
the number of lambda values - Default is 100.
}
  \item{nb_group}{
the number of groups for computing the Qini coefficient - Default is 10.
}
  \item{value}{
if TRUE, the values of the best lambda and qini coefficient will be printed - Default is FALSE.
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
library(tools4uplift)
scene <- generateScene(plotit=FALSE)
data <- generateUplift(scene, 1000, 10)

features <- bestFeatures(data, "y", "treat", colnames(data[,3:12]), nb_lambda = 100, nb_group = 10, value = FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
