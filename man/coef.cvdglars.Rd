\name{coef.cvdglars}
\alias{coef.cvdglars}
\title{Extract the coefficients estimated by \code{cvdglars}}
\description{\code{coef.cvdglars} is used to extract the coefficients estimated by the \eqn{k}-fold cross-validation deviance.}
\usage{
\method{coef}{cvdglars}(object, ...)
}
\arguments{
  \item{object}{fitted \code{cvdglars} object}
  \item{\dots}{additional arguments}
}
\value{
\code{coef.cvdglars} returns a vector with the coefficients estimated by the \eqn{k}-fold cross-validation deviance.
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
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- cvdglars.fit(X, y, family = "binomial")
coef(fit)
}
\keyword{models}
\keyword{regression}
