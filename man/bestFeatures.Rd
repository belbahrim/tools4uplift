\name{bestFeatures}
\alias{bestFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature Selection for Uplift Penalized Logistic Regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bestFeatures(data, outcome, treat, predictors, nb_lambda=100, nb_group=10, value=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{training} here~~
}
  \item{outcome}{
%%     ~~Describe \code{formula} here~~
}
  \item{treat}{
%%     ~~Describe \code{treat} here~~
}
  \item{predictors}{
%%     ~~Describe \code{treat} here~~
}
  \item{nb_lambda}{
%%     ~~Describe \code{nb_lambda} here~~
}
  \item{nb_group}{
%%     ~~Describe \code{nb_group} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
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

features <- bestFeatures(data, "y", "treat", colnames(data[,3:12]) nb_lambda = 100, nb_group = 10, value = FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
