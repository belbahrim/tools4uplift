\name{predict.BinUplift}
\alias{predict.BinUplift}
\title{
Prediction from univariate quantization
}
\description{
Predictions from the univariate quantization method, i.e. this function transforms a continuous variable into a categorical one.
}
\usage{

\method{predict}{BinUplift}(object, newdata, ...)
}
\arguments{
  \item{object}{an object of class \code{BinUplift}, as that created by the function \code{BinUplift}.
}
  \item{newdata}{the variable that was quantized in \code{object}.
}
  \item{\dots}{additional arguments to be passed to \code{cut} function.
}
}
\value{
a quantized variable
}
\references{
Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\seealso{
\code{BinUplift}
}
\examples{

library(tools4uplift)
data("SimUplift")

binX1 <- BinUplift(data = SimUplift, treat = "treat", outcome = "y", x = "X1", 
                  n.split = 100, alpha = 0.01, n.min = 30)

quantizedX1 <- predict(binX1, SimUplift$X1)

}
\keyword{quantization}
\keyword{uplift}
