\name{QiniTable}
\alias{QiniTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the performance of an Uplift model
}
\description{
Table of performance of an uplift model. This table is used in order to vizualise the performance of an uplift model and to compute the qini coefficient.
}
\usage{
QiniTable(data, treat, outcome, prediction, nb.group = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data frame containing the response, the treatment and predicted uplift.
}
  \item{treat}{
a binary (numeric) vector representing the treatment assignment (coded as 0/1).
}
  \item{outcome}{
a binary response (numeric) vector (coded as 0/1).
}
  \item{prediction}{
a predicted uplift (numeric) vector to sort the observations from highest to lowest uplift.
}
  \item{nb.group}{
number of groups of equal observations in which to partition the data set to show results.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

split <- SplitUplift(train, 0.5, c("treat", "y"))
train <- split[[1]]
valid <- split[[2]]

square <- SquareUplift(train, valid, "X1", "X2", "treat", "y", 
                      step = 0.3, nb.group = 5, plotit = FALSE)

#performance of the heat map uplift estimation on the training dataset
perf <- QiniTable(data = square[[1]], treat = "treat", 
                  outcome = "y", prediction = "uplift", nb.group = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
