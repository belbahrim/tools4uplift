\name{barplot.PerformanceUplift}
\alias{barplot.PerformanceUplift}
\title{
Uplift barplot
}
\description{
Barplot of observed uplift with respect to predicted uplift sorted from the highest to the lowest.
}
\usage{

\method{barplot}{PerformanceUplift}(height, ...)
}
\arguments{
  \item{height}{
a table that must be the output of \code{PerformanceUplift} function.
}
  \item{...}{
additional barplot arguments.
}
}
\value{
a barplot and the associated Kendall's uplift rank correlation
}
\references{
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

barplot(perf)

}
\keyword{performance}
