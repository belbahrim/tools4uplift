\name{qLHS}
\alias{qLHS}
\title{
Qini-based uplift regression
}
\description{
A Qini-based LHS (Latin Hypercube Sampling) uplift model.
}
\usage{
qLHS(data, treat, outcome, predictors, 
     lhs_points = 50, lhs_range = 1, 
     adjusted = TRUE, rank.precision = 2, equal.intervals = FALSE, 
     nb.group = 10, validation = TRUE, p = 0.3)
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
  \item{lhs_points}{
number of LHS points to sample for each regularization constant.
}
  \item{lhs_range}{
a multiplicative scalar that controls the variance of the LHS search - Default is 1, the LHS procedure samples points uniformly with variance equal to the variance of the maximum likelihood estimator.
}
  \item{adjusted}{
if TRUE, the adjusted Qini coefficient is used instead of the Qini coefficient.
}
  \item{rank.precision}{
precision for the ranking quantiles to compute the Qini coefficient. Must be 1 or 2. If 1, the ranking quantiles will be rounded to the first decimal. If 2, to the second decimal.
}
  \item{equal.intervals}{
flag for using equal intervals (with equal number of observations) or the true ranking quantiles which result in an unequal number of observations in each group to compute the Qini coefficient.
}
  \item{nb.group}{
the number of groups for computing the Qini coefficient if equal.intervals is TRUE - Default is 10.
}
  \item{validation}{
if TRUE, the best LHS model is selected based on cross-validation - Default is TRUE.
}
  \item{p}{
if validation is TRUE, the desired proportion for the validation set. p is a value between 0 and 1 expressed as a decimal, it is set to be proportional to the number of observations per group - Default is 0.3.
}
}
\details{
The regularization parameter is chosen based on the interaction uplift model that maximizes the Qini coefficient of the LHS search.
}
\value{
the models with LHS coefficients of class \code{InterUplift}.
}
\references{
Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2020) Qini-based Uplift Regression, 
<https://arxiv.org/pdf/1911.12474.pdf>
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
data("SimUplift")

upliftLHS <- qLHS(data = SimUplift, treat = "treat", outcome = "y", 
                  predictors = colnames(SimUplift[,3:7]), lhs_points = 5,
                  lhs_range = 1, adjusted = TRUE, equal.intervals = TRUE, 
                  nb.group = 5, validation = FALSE)

}
\keyword{lhs}
\keyword{LASSO}
\keyword{Qini}
\keyword{logistic}
