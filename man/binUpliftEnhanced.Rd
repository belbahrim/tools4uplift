\name{binUpliftEnhanced}
\alias{binUpliftEnhanced}
\docType{data}
\title{
Augmented Data With Binned Variables
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
binUpliftEnhanced(data, 
                  var_list, 
                  response_var_name, 
                  group_var_name, 
                  n.split=10, alpha=0.05, 
                  n.min=30, ylim=NULL, 
                  ylab="Uplift", 
                  title="Binning Results",
                  ...)
}
\arguments{
  \item{data}{
a data frame containing the response, the treatment and the variable to bin.
}
  \item{var_list}{
a list representing the explanatory variables to bin.
}
  \item{response_var_name}{
a binary response (numeric) vector (coded as 0/1).
}
  \item{group_var_name}{
a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{n.split}{
number of splits to test at each node. For continuous varX only (must be > 0).
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
  \item{...}{
Additional arguments passed to plot.
}
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\value{
%%  ~Describe the value returned
an augmented data frame with categorized variables. If a variable is enhanced, the function returns automatically a barplot.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

library(tools4uplift)

#generate 1000 observations with 10 predictors
scene <- generateScene(plotit=FALSE)
data <- generateUplift(scene, 1000, 10)

#augment dataset with categorized variables
enhanced_data <- binUpliftEnhanced(data, 
                  var_list=colnames(data[,3:12]), 
                  response_var_name="y", 
                  group_var_name="treat", 
                  n.split=10, 
                  alpha=0.05, 
                  n.min=30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
