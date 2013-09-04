\name{print.cvdglars}
\alias{print.cvdglars}
\title{Print a cvdglars object}
\description{
Print information about the points of the dgLARS solution curve identified by the \eqn{k}-fold cross-validation deviance.
}
\usage{
\method{print}{cvdglars}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{dglars} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
The call that produced the object \code{x} is printed, followed by the estimate of the coefficients of the variables included in the active set. Such
estimates are obtained using the whole data set while the optimal value of the tuning parameter is estimated by \eqn{k}-fold cross-validation 
deviance. The last part of the print gives information about the algorithm and the method used to compute the solution curve. The code about 
the convergence of the used algorithm is also showed.
}
\value{The vector of the estimates is silently returned.}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{cvdglars}} function.
}
\examples{
######################
# Logistic regression model

set.seed(123)

n <- 100
p <- 10
X <- matrix(rnorm(n*p) ,n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- cvdglars.fit(X, y, family = "binomial")
fit
}
\keyword{models}
\keyword{regression}
