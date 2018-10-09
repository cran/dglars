\name{coef.cvdglars}
\alias{coef.cvdglars}

\title{Extract the Coefficients Estimated by \code{cvdglars}}

\description{\code{coef.cvdglars} is used to extract the coefficients estimated by \eqn{k}-fold cross-validation deviance.}

\usage{
\method{coef}{cvdglars}(object, ...)
}

\arguments{
\item{object}{fitted \code{cvdglars} object}
\item{\dots}{additional argument used to ensure the compatibility with the generic method function \dQuote{\code{\link{coef}}}.}
}

\value{
\code{coef.cvdglars} returns a named list with components \code{beta}, i.e., the estimate of the coefficient vector, and \code{phi} the estimate of the dispersion parameter.
}

\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}

\seealso{
\code{\link{cvdglars}} function.
}

\examples{
###########################
# Logistic regression model
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- cvdglars.fit(X, y, family = binomial)
coef(fit)
}

\keyword{models}
\keyword{regression}
