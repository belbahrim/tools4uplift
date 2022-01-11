\name{predict.InterUplift}
\alias{predict.InterUplift}
\title{
Predictions from an interaction estimator
}
\description{
Predictions from the interaction uplift model estimator with associated model performance.
}
\usage{

\method{predict}{InterUplift}(object, newdata, treat, ...)

}
\arguments{
  \item{object}{an object of class \code{InterUplift}, as that created by the function \code{InterUplift}.
}
  \item{newdata}{a data frame containing the treatment, the outcome and the predictors of observations at which predictions are required.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{\dots}{additional arguments to be passed to \code{predict.glm} function for the interaction model.
}
}
\value{
a vector of predicted uplift
}
\references{
Lo, V., S., Y. (2002) The true lift model: a novel data mining approach to response modeling in database marketing. ACM SIGKDD Explorations Newsletter, Vol. 4(2), 78-86.

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\seealso{
\code{InterUplift}
}
\examples{

library(tools4uplift)
data("SimUplift")

fit <- InterUplift(SimUplift, "treat", "y", colnames(SimUplift[, 3:12]))

pred <- predict(fit, SimUplift, "treat")

}
\keyword{prediction}
\keyword{uplift}
