\name{simUplift}
\alias{simUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates synthetic data from a randomized experiment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simUplift(n, p, c0, c1, sigma0, sigma1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{c0}{
%%     ~~Describe \code{c0} here~~
}
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{sigma0}{
%%     ~~Describe \code{sigma0} here~~
}
  \item{sigma1}{
%%     ~~Describe \code{sigma1} here~~
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
Mouloud Belbahri <mouloud.belbahri@gmail.com>
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
simUplift(n=10, p=2, c0=1, c1=5, sigma0=0.1, sigma1=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
