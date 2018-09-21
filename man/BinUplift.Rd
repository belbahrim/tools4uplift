\name{BinUplift}
\alias{BinUplift}
\title{
Optimal Partitionning for Uplift Models
}
\description{
Univariate optimal partitionning for Uplift Models. The algorithm categorizes a single variable into bins with significantly different uplift values.
}
\usage{
BinUplift(data, treat, outcome, x, n.split = 10, alpha = 0.05, 
                       n.min = 30, ylim = NULL, ylab = "Uplift",
                       title = "Binning Results", color = NULL)
}
\arguments{
  \item{data}{
a data frame containing the treatment, the outcome and the predictor to categorize.
}
  \item{treat}{
name of a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
name of a binary response (numeric) vector (coded as 0/1).
}
  \item{x}{
name of the explanatory variable to categorize.
}
  \item{n.split}{
number of splits to test at each node. For continuous explanatory variables only (must be > 0).
}
  \item{alpha}{
significance level of the statistical test (must be between 0 and 1).
}
  \item{n.min}{
minimum number of observations per child node.
}
  \item{ylim}{
a range for the y axis.
}
  \item{ylab}{
a title for the y axis.
}
  \item{title}{
an overall title for the plot.
}
  \item{color}{
a color for the plot. If ommitted, the color will be set by default to a custom light blue. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mouloud Belbahri
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~ BinUpliftEnhanced
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

binX1 <- BinUplift(data = train, treat = "treat", outcome = "y", x = "X1", 
                  alpha = 0.10, n.min = 3, title = "Binning for X1")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
