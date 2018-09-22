\name{InterUplift}
\alias{InterUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The interaction Uplift Model Estimator
}
\description{
Fit the interaction uplift model estimator.
}
\usage{
InterUplift(data, treat, outcome, predictors, input = c("p", "bf"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
}
  \item{outcome}{
%%     ~~Describe \code{outcome} here~~
}
  \item{predictors}{
%%     ~~Describe \code{predictors} here~~
}
  \item{input}{
%%     ~~Describe \code{input} here~~
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
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

fit <- InterUplift(train, "treat", "y", data[,3:12], "p")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
