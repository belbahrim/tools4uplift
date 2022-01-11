\name{InterUplift}
\alias{InterUplift}
\alias{InterUplift.default}
\alias{InterUplift.formula}
\title{
Interaction estimator
}
\description{
Fit the interaction uplift model estimator.
}
\usage{

\method{InterUplift}{formula}(formula, treat, data, ...)

\method{InterUplift}{default}(data, treat, outcome, predictors, input = "all", ...)

}
\arguments{
  \item{data, formula}{a data frame containing the treatment, the outcome and the predictors or a formula describing the model to be fitted.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{predictors}{
a vector of names representing the explanatory variables to include in the model.
}
  \item{input}{
an option for \code{predictors} argument. If \code{"all"} (default), the model assumes that the model has to create the interaction of all varibles with \code{treat}. If \code{"best"}, the model assumes that the \code{predictors} vector is the output of the \code{BestFeatures} function.
}
  \item{\dots}{additional arguments (other than \code{formula}, \code{family}, and \code{data}) to be passed to \code{glm} function for the interaction model.
}
}
\value{
an interaction model
}
\references{
Lo, V., S., Y. (2002) The true lift model: a novel data mining approach to response modeling in database marketing. ACM SIGKDD Explorations Newsletter, Vol. 4(2), 78-86.

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
data("SimUplift")

fit <- InterUplift(SimUplift, "treat", "y", colnames(SimUplift[, 3:12]))

}
\keyword{uplift}
