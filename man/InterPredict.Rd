\name{interPredict}
\alias{interPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions from the interaction Uplift Model Estimator
}
\description{
Predictions from the interaction uplift model estimator with associated model performance.
}
\usage{
interPredict(data, treat, outcome, model, nb_group = 10, plotit = FALSE)
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
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{nb_group}{
%%     ~~Describe \code{nb_group} here~~
}
  \item{plotit}{
%%     ~~Describe \code{plotit} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
