\name{DualUplift}
\alias{DualUplift}
\alias{DualUplift.default}
\alias{DualUplift.formula}
\alias{print.DualUplift}
\alias{summary.DualUplift}
\title{
Two-model estimator
}
\description{
Fit the two-model uplift model estimator.
}
\usage{

\method{DualUplift}{formula}(formula, treat, data, ...)

\method{DualUplift}{default}(data, treat, outcome, predictors, ...)

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
  \item{\dots}{additional arguments (other than \code{formula}, \code{family}, and \code{data}) to be passed to \code{glm} function for each sub-model.
}
}
\value{

\item{model0}{Fitted model for control group}
\item{model1}{Fitted model for treatment group}

}
\references{
Hansotia, B., J., and Rukstales B. (2001) Direct marketing for multichannel retailers: Issues, challenges and solutions. Journal of Database Marketing and Customer Strategy Management, Vol. 9(3), 259-266.

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
data("SimUplift")

fit <- DualUplift(SimUplift, "treat", "y", predictors = colnames(SimUplift[, 3:12]))

print(fit)
summary(fit)
}
\keyword{uplift}
