\name{bayesPredict}
\alias{bayesPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions from a Fitted Bayesian Direct Uplift Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bayesPredict(object, newdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class as that created by the function \code{bayesUplift}.
}
  \item{newdata}{
a data frame containing the values at which predictions are required.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
an augmented data frame with predicted uplift for each observation.
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
#do not output the functions plots
scene <- generateScene(plotit=FALSE)

#generate 1000 observations with 20 predictors
train <- generateUplift(scene, 1000, 20)
valid <- generateUplift(scene, 1000, 20)

#fit a Bayes Uplift Model on training data frame
fit <- bayesUplift(train, outcome="y", treatment="treat", predictors=colnames(data[,3:22]), folds=10, prior=FALSE)

#predict on a validation data frame
valid_pred <- bayesPredict(fit, newdata=valid)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{uplift}% use one of  RShowDoc("KEYWORDS")
\keyword{bayes}% __ONLY ONE__ keyword per line
