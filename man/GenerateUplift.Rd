\name{GenerateUplift}
\alias{GenerateUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulations for Personalized Treatment Effects
}
\description{
Numerical simulation for treatment effect heterogeneity estimation.
}
\usage{
GenerateUplift(object, nbobs, dimension)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a scene of class as that created by the function \code{generateScene}.
}
  \item{nbobs}{
number of observations.
}
  \item{dimension}{
number of predictors.
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
Olivier Gandouet
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)
train <- GenerateUplift(scene, nbobs = 1000, dimension = 10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
