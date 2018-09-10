\name{binUplift}
\alias{binUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal partitionning for Uplift Models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binUplift(df, 
          response, 
          group, 
          varX, 
          n.split = 10, 
          alpha=0.05, 
          n.min=30, 
          ylim = NULL, 
          ylab = "Uplift", 
          title="Binning Results", 
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
a data frame containing the response, the treatment and the variable to bin.
}
  \item{response}{
a binary response (numeric) vector (coded as 0/1).
}
  \item{group}{
a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{varX}{
a vector representing the explanatory variable to bin.
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(tools4uplift)
scene <- generateScene(plotit=FALSE)
train <- generateUplift(scene, 1000, 10)

binX1 <- binUplift(df=train, response="y", group="treat", varX="X1", alpha=0.10, n.min=3, title="Binning for X1")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
