\name{UpliftPerCat}
\alias{UpliftPerCat}
\title{
Uplift barplot for categorical variables
}
\description{
Computes the observed uplift per category and creates a barplot.
}
\usage{
UpliftPerCat(data, treat, outcome, x, ...)
}
\arguments{
  \item{data}{
a data frame containing the treatment, the outcome and the variable of interest.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{x}{
name of the explanatory variable of interest.
}
  \item{...}{extra parameters for the barplot.
}
}
\value{
returns a barplot representing the uplift per category.
}
\author{
Mouloud Belbahri
}
\examples{

library(tools4uplift)
data("SimUplift")

binX1 <- BinUplift(data = SimUplift, treat = "treat", outcome = "y", x = "X1", 
                  n.split = 100, alpha = 0.01, n.min = 30)

SimUplift$quantizedX1 <- predict(binX1, SimUplift$X1)
UpliftPerCat(data = SimUplift, treat = "treat", outcome = "y", 
            x = "quantizedX1", xlab='Quantized X1', ylab='Uplift',
            ylim=c(-1,1))

}
\keyword{barplot}
