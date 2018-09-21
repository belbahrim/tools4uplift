\name{BinUpliftEnhanced}
\alias{BinUpliftEnhanced}
\title{
Optimal Partitionning for Uplift Models Augmented Dataset 
}
\description{
Univariate optimal partitionning for Uplift Models. The algorithm categorizes several variables into bins and creates an augmented dataset with the binned variables.
}
\usage{
BinUpliftEnhanced(data, treat, outcome, var.list, n.split = 10, alpha = 0.05, 
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
  \item{var.list}{
a vector of names representing the explanatory variables to categorize.
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\value{
%%  ~Describe the value returned
an augmented data frame with categorized variables. If a variable is enhanced, the function returns automatically a barplot.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\author{
Olivier Gandouet
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~ BinUpliftEnhanced
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

train.enhanced <- BinUpliftEnhanced(data = train, treat = "treat", outcome = "y", 
                                    var.list = colnames(train[,3:12]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line