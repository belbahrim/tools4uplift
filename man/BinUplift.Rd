\name{BinUplift}
\alias{BinUplift}
\title{
Univariate quantization
}
\description{
Univariate optimal partitionning for Uplift Models. The algorithm quantizes a single variable into bins with significantly different observed uplift.
}
\usage{
BinUplift(data, treat, outcome, x, n.split = 10, alpha = 0.05, n.min = 30)
}
\arguments{
  \item{data}{
a data frame containing the treatment, the outcome and the predictor to quantize.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{x}{
name of the explanatory variable to quantize.
}
  \item{n.split}{
number of splits to test at each node. For continuous explanatory variables only (must be > 0). If n.split = 10, the test will be executed at each decile of the variable.
}
  \item{alpha}{
significance level of the statistical test (must be between 0 and 1).
}
  \item{n.min}{
minimum number of observations per child node.
}
}
\value{

\item{out.tree}{Descriptive statistics for the different nodes of the tree}
}
\references{
Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\seealso{
\code{predict.BinUplift}
}
\examples{

library(tools4uplift)
data("SimUplift")

binX1 <- BinUplift(data = SimUplift, treat = "treat", outcome = "y", x = "X1", 
                  n.split = 100, alpha = 0.01, n.min = 30)

}
\keyword{tree}
\keyword{quantization}
