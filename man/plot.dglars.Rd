\name{plot.dglars}
\alias{plot.dglars}
\title{Plot from a dglars object}
\description{
Produces different plots to study the sequence of models identified by the dgLARS method.
}
\usage{
\method{plot}{dglars}(x, ...)
}
\arguments{
  \item{x}{fitted \code{dglars} object.}
  \item{\dots}{additional graphical parameters.}
}
\details{
\code{plot.dglars} method produces different plots to study the sequence of models identified by the dgLARS method. The first two plots
show the Bayesian Information Criterion (BIC) and the Akaike Information Criterion (AIC) as functions of the used \eqn{\gamma} values. The number
of non zero coefficients is used as measure of the complexity of a model. The value of the tuning parameter corresponding to the minimum of the 
BIC/AIC curve is labeled by a vertical dashed red line, while the \eqn{\gamma} values at which corresponds a change in the active set are labeled by 
vertical dashed gray lines. Next plot shows the coefficient profile plot while, if the control parameter \code{algorithm} is equal to \code{"pc"}, the last 
plot shows the Rao's score test statistics as functions of \eqn{\gamma}. Also in these two plots, the optimal \eqn{\gamma} and the \eqn{\gamma} values 
at which corresponds a change in the active set are labeled by vertical dashed lines, red and gray respectively.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} function.
}
\examples{
######################
# Logistic regression model

n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit_pc <- dglars(X, y, family = "binomial")
par(mfrow = c(2, 2))
plot(fit_pc)
fit_ccd <- dglars(X, y, family = "binomial", control = list(algorithm = "ccd"))
par(mfrow = c(1, 3))
plot(fit_ccd)
}
\keyword{models}
\keyword{regression}
