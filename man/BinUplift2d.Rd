\name{BinUplift2d}
\alias{BinUplift2d}
\title{
Bivariate quantization
}
\description{
A non-parametric heat map representing the observed uplift in rectangles that explore a bivariate dimension space. The function also returns the individual uplift based on the heatmap. 
}
\usage{
BinUplift2d(data, var1, var2, treat, outcome, valid = NULL, 
            n.split = 3, n.min = 30, plotit = FALSE, nb.col = 20)
}
\arguments{
  \item{data}{
a data frame containing uplift models variables. 
}
  \item{var1}{
x-axis variable name. Represents the first dimension of interest.
}
  \item{var2}{
y-axis variable name. Represents the second dimension of interest.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{valid}{
a validation data frame containing uplift models variables. 
}
  \item{n.split}{
the number of intervals to consider per explanatory variable. Must be an integer > 1. 
}
  \item{n.min}{
minimum number of observations per group (treatment and control) within each rectangle. Must be an integer > 0.
}
  \item{plotit}{
if TRUE, a heatmap of observed uplift per rectangle is plotted.
}
  \item{nb.col}{
number of colors for the heatmap. Default is 20. Must be an integer and should greater than \code{n.split} for better visualization.
}
}
\value{
returns an augmented dataset with \code{Uplift_var1_var2} variable representing a predicted uplift for each observation based on the rectangle it belongs to. The function also plots a heat map of observed uplifts.
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

heatmap <- BinUplift2d(SimUplift, "X1", "X2", "treat", "y")

}
\keyword{heatmap}
\keyword{uplift}
