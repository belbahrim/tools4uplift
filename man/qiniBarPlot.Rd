\name{qiniBarPlot}
\alias{qiniBarPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the Uplift barplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qiniBarPlot(x, title = "Model Performance: Uplift by Group")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a table that must be the output of qiniTable function.
}
  \item{title}{
an overall title for the plot.
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

square <- squareUplift(train, valid, "X1", "X2", "treat", "y", step=0.3, nbGroup=5, plotit=FALSE)

#performance of the heat map uplift estimation on the training dataset
perf <- qiniTable(data=square[[1]], nbGroup=4, treatment="treat", outcome="y", prediction="uplift")

qiniBarPlot(perf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
