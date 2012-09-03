\name{coef.dglars}
\alias{coef.dglars}
\title{
Extract the path of the coefficients
}
\description{
\code{coef.dglars} is used to extract the path of the coefficients computed using the dgLARS method.
}
\usage{
\method{coef}{dglars}(object, ...)
}
\arguments{
  \item{object}{fitted \code{dglars} object}
  \item{\dots}{additional arguments}
}
\value{
\code{coef.dglars} returns a matrix with the solution path computed using the dgLARS method.}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} function.
}
\examples{
######################
# Logistic regression model

set.seed(123)

n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(X, y, family = "binomial")
coef(fit)
}
\keyword{models}
\keyword{regression}
