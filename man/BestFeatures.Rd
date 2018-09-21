\name{BestFeatures}
\alias{BestFeatures}
\title{
Feature Selection for the Interaction Uplift Model Estimator
}
\description{
Penalized logistic regression (LASSO) in order to select the features that maximize the Qini coefficient. 
}
\usage{
BestFeatures(data, treat, outcome, predictors, nb.lambda = 100, nb.group = 10, value = FALSE)
}
\arguments{
  \item{data}{
a data frame containing the treatment, the outcome and the predictors.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{predictors}{
a vector of names representing the predictors to consider in the model.
}
  \item{nb.lambda}{
the number of lambda values - Default is 100.
}
  \item{nb.group}{
the number of groups for computing the Qini coefficient - Default is 10.
}
  \item{value}{
if TRUE, the values of the best lambda and Qini coefficient will be printed - Default is FALSE.
}
}
\details{
The regularization parameter is chosen based on the interaction uplift model that maximizes the Qini coefficient. Using the LASSO penalty, some predictors have coefficients set to zero.
}
\value{
a vector of names representing the selected best features from the penalized logistic regression.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

features <- BestFeatures(data = train, treat = "treat", outcome = "y", 
                         predictors = colnames(data[,3:12]))
features

}
\keyword{LASSO}