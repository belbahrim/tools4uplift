\name{QiniArea}
\alias{QiniArea}
\alias{QiniArea.default}
\alias{QiniArea.PerformanceUplift}
\title{
Qini coefficient
}
\description{
Computes the area under the Qini curve.
}
\usage{
\method{QiniArea}{PerformanceUplift}(x, adjusted=FALSE, ...)
}
\arguments{
  \item{x}{
a table that must be the output of \code{PerformanceUplift} function.
}
  \item{adjusted}{
if TRUE, returns the Qini coefficient adjusted by the Kendall's uplift rank correlation.
}
  \item{...}{Generic S3 Method argument.
}
}
\value{
the Qini or the adjusted Qini coefficient
}
\references{
Radcliffe, N. (2007). Using control groups to target on predicted lift: Building and assessing uplift models. Direct Marketing Analytics Journal, An Annual Publication from the Direct Marketing Association Analytics Council, pages 14-21.

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2020) Qini-based Uplift Regression, 
<https://arxiv.org/pdf/1911.12474.pdf>

Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\seealso{
\code{PerformanceUplift}
}
\examples{

library(tools4uplift)
data("SimUplift")

model <- BinUplift2d(SimUplift, "X1", "X2", "treat", "y")

#performance of the heat map uplift estimation on the training dataset
perf <- PerformanceUplift(data = model, treat = "treat", 
                  outcome = "y", prediction = "Uplift_X1_X2", 
                  equal.intervals = TRUE, nb.group = 5)
QiniArea(perf)

}
\keyword{Qini}
\keyword{performance}
