\name{SplitUplift}
\alias{SplitUplift}
\title{
Split data with respect to uplift distribution
}
\description{
Split a dataset into training and validation subsets with respect to the uplift sample distribution.
}
\usage{
SplitUplift(data, p, group)
}
\arguments{
  \item{data}{
a data frame of interest that contains at least the response and the treatment variables.
}
  \item{p}{
The desired sample size. p is a value between 0 and 1 expressed as a decimal, it is set to be proportional to the number of observations per group.
}
  \item{group}{
Your grouping variables. Generally, for uplift modelling, this should be a vector of treatment and response variables names, e.g. c("treat", "y").
}
}
\value{
\item{train}{a training data frame of p percent}
\item{valid}{a validation data frame of 1-p percent}
}
\references{
Belbahri, M., Murua, A., Gandouet, O., and Partovi Nia, V. (2021) Uplift Regression : The R Package tools4uplift,
<https://arxiv.org/pdf/1901.10867.pdf>
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
data("SimUplift")

split <- SplitUplift(SimUplift, 0.8, c("treat", "y"))
train <- split[[1]]
valid <- split[[2]]

}
\keyword{sampling}