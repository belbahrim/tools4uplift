\name{qiniTable}
\alias{qiniTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the performance of an Uplift model
}
\description{
Table of performance of an uplift model. This table is used in order to vizualise the performance of an uplift model and to compute the qini coefficient.
}
\usage{
qiniTable(data, nbGroup, treatment, outcome, prediction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame containing the response, the treatment and predicted uplift.
}
  \item{nbGroup}{
number of groups of equal observations in which to partition the data set to show results.
}
  \item{treatment}{
a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
a binary response (numeric) vector (coded as 0/1).
}
  \item{prediction}{
a predicted uplift (numeric) vector to sort the observations from highest to lowest uplift.
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
data <- generateUplift(scene, 1000, 10)

split <- splitUplift(data, 0.5, c("treat", "y"))
train <- split[[1]]
valid <- split[[2]]

square <- squareUplift(train, valid, "X1", "X2", "treat", "y", 
                      step=0.3, nbGroup=5, plotit=FALSE)

#performance of the heat map uplift estimation on the training dataset
perf <- qiniTable(data=square[[1]], nbGroup=4, treatment="treat", 
                  outcome="y", prediction="uplift")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
