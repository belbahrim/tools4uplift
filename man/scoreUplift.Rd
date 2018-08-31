\name{scoreUplift}
\alias{scoreUplift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Uplift Model Performance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreUplift(data, treat, outcome, model, nb_group = 10, plotit = FALSE)
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, treat, outcome, model, nb_group = 10, plotit = FALSE) 
{
    data1 <- data
    data1[[treat]] <- 1
    pr.y1_ct1 <- predict(model, newdata = data1, type = "response")
    data0 <- data
    data0[[treat]] <- 0
    pr.y1_ct0 <- predict(model, newdata = data0, type = "response")
    data$uplift_prediction <- pr.y1_ct1 - pr.y1_ct0
    perf <- qiniTable(data, nb_group, treat, outcome, "uplift_prediction")
    if (plotit == TRUE) {
        qiniCurve(perf, title = "Qini Curve for proposed method model")
        qiniBarPlot(perf, title = "Uplift Bar Plot for proposed method model")
    }
    qini <- qiniArea(perf)
    return(qini)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
