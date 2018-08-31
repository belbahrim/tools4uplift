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
bestFeatures(training, formula, treat, nb_lambda = 100, nb_group = 10, value = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training}{
%%     ~~Describe \code{training} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{treat}{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (training, formula, treat, nb_lambda = 100, nb_group = 10, 
    value = FALSE) 
{
    path <- lassoPath(training, formula, nb_lambda)
    path <- path[path[, "dimension"] > 0, ]
    lambda_qini <- c()
    for (k in 1:nrow(path)) {
        features <- path[k, -c(1, 2)]
        features <- features[features != 0]
        lambda_formula <- paste(paste(colnames(model.frame(formula, 
            training))[1], "~ "), paste(names(features), collapse = " + "))
        lambda_model <- glm(lambda_formula, family = binomial(link = "logit"), 
            training)
        lambda_qini[k] <- scoreUplift(training, treat, colnames(model.frame(formula, 
            training))[1], lambda_model, nb_group)
    }
    best_model <- cbind(path[, c(1, 2)], lambda_qini)
    best_model <- best_model[which.max(best_model[, "lambda_qini"]), 
        ]
    best_features <- path[path[, "lambda"] == best_model["lambda"], 
        -c(1, 2)]
    best_features <- names(best_features[best_features != 0])
    if (value == TRUE) {
        print(best_model)
    }
    return(best_features)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
