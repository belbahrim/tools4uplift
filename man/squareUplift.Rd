\name{squareUplift}
\alias{squareUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-parametric heat map for Uplift
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
squareUplift(train, valid, var1, var2, treat, outcome, step, 
         n.min=1, 
         plotit = TRUE, nbGroup = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{
a training data frame containing uplift models variables. 
}
  \item{valid}{
a validation data frame containing uplift models variables. 
}
  \item{var1}{
x-axis variable name. Represents the first dimension of interest.
}
  \item{var2}{
y-axis variable name. Represents the second dimension of interest.
}
  \item{treat}{
a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
a binary response (numeric) vector (coded as 0/1).
}
  \item{step}{
the size of a side of the squares.
}
  \item{n.min}{
minimum number of observations per square.
}
  \item{plotit}{
if TRUE, the heatmap will be plotted for both training and validation datasets.
}
  \item{nbGroup}{
number of groups of equal observations in which to partition the data set to create a cluster new variable.
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
library("tools4uplift")
set.seed(12345)
train <- simUplift(n=100, p=2, c0=1, c1=5, sigma0=0.1, sigma1=0.1)
valid <- simUplift(n=100, p=2, c0=1, c1=5, sigma0=0.1, sigma1=0.1)

square <- squareUplift(train, valid, "X1", "X2", "treat", "y", step=0.3, nbGroup=5, plotit=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
