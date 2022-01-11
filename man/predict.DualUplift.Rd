\name{predict.DualUplift}
\alias{predict.DualUplift}
\title{
Predictions from a two-model estimator
}
\description{
Predictions from the two-model uplift model estimator with associated model performance.
}
\usage{

\method{predict}{DualUplift}(object, newdata, ...)
}
\arguments{
  \item{object}{an object of class \code{DualUplift}, as that created by the function \code{DualUplift}.
}
  \item{newdata}{a data frame containing the treatment, the outcome and the predictors of observations at which predictions are required.
}
  \item{\dots}{additional arguments to be passed to \code{predict.glm} function for each sub-model.
}
}
\value{
a vector of predicted uplift
}
\references{
Hansotia, B., J., and Rukstales B. (2001) Direct marketing for multichannel retailers: Issues, challenges and solutions. Journal of Database Marketing and Customer Strategy Management, Vol. 9(3), 259-266.

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\seealso{
\code{DualUplift}
}
\examples{

library(tools4uplift)
data("SimUplift")

fit <- DualUplift(SimUplift, "treat", "y", predictors = colnames(SimUplift[, 3:12]))

pred <- predict(fit, SimUplift)

}
\keyword{prediction}
\keyword{uplift}
