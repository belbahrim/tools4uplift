\name{GenerateScene}
\alias{GenerateScene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scenarios for Personalized Treatment Effects
}
\description{
Creation of scenarios for treatment effect heterogeneity estimation.
}
\usage{
GenerateScene(space_dimension = 2, 
              Principal_effect = c(0.1, 1), 
              Principal_std = c(1, 4), 
              Control_effect_positive = c(0.2, 0.3), 
              Control_std_positive = c(1.5, 2), 
              Control_effect_negative = c(0.3, 0.4), 
              Control_std_negative = c(1.5, 2), 
              Nb_point = c(20, 10, 20), 
              plotit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{space_dimension}{
%%     ~~Describe \code{space_dimension} here~~
}
  \item{Principal_effect}{
%%     ~~Describe \code{Principal_effect} here~~
}
  \item{Principal_std}{
%%     ~~Describe \code{Principal_std} here~~
}
  \item{Control_effect_positive}{
%%     ~~Describe \code{Control_effect_positive} here~~
}
  \item{Control_std_positive}{
%%     ~~Describe \code{Control_std_positive} here~~
}
  \item{Control_effect_negative}{
%%     ~~Describe \code{Control_effect_negative} here~~
}
  \item{Control_std_negative}{
%%     ~~Describe \code{Control_std_negative} here~~
}
  \item{Nb_point}{
%%     ~~Describe \code{Nb_point} here~~
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
Olivier Gandouet
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(tools4uplift)
scene <- GenerateScene(plotit = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
