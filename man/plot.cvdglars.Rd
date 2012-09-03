\name{plot.cvdglars}
\alias{plot.cvdglars}
\title{Plot from a cvdglars object}
\description{
Plots the cross-validation curve, and upper and lower standard deviation curves, as a function of the used \eqn{\gamma} values.
}
\usage{
\method{plot}{cvdglars}(x, ...)
}
\arguments{
  \item{x}{fitted \code{cvdglars} object.}
  \item{\dots}{additional graphical parameters to plot.}
}
\details{
A plot for a \code{cvdglars} object is produced. The plot shows the curve of the cross-validation deviance and the upper and lower standard deviation curves.
The optimal value of the \eqn{\gamma} parameter, which corresponds to the minimum of the cross-validation deviance, is identified by a vertical dashed red line.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}\seealso{
\code{\link{cvdglars}} function.
}
\examples{
######################
# Logistic regression model
#
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- cvdglars(X, y, family = "binomial")
plot(fit)
}
\keyword{models}
\keyword{regression}
