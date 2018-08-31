\name{qiniArea}
\alias{qiniArea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Qini coefficient
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qiniArea(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a table that must be the output of qiniTable function.
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

square <- squareUplift(train, valid, "X1", "X2", "treat", "y", step=0.3, nbGroup=5, plotit=FALSE)

#performance of the heat map uplift estimation on the training dataset
perf <- qiniTable(data=square[[1]], nbGroup=4, treatment="treat", outcome="y", prediction="uplift")

qiniArea(perf)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
