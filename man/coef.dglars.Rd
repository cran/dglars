\name{coef.dglars}
\alias{coef.dglars}
\title{Extract the dgLARS Coefficient Path}
\description{\code{coef.dglars} is used to extract the coefficient path 
computed by dgLARS method.}
\usage{
\method{coef}{dglars}(object, type = c("pearson", "deviance", "mle"), 
\dots)
}
\arguments{
\item{object}{fitted \code{dglars} object.}
\item{type}{a description of the estimator used for the dispersion parameter.}
\item{\dots}{this is the argument by means of to pass
the argument \code{g} to \code{predict.dglars} and \code{phihat}.}
}
\details{
\code{coef.dglars} is a wrapper function calling \dQuote{\code{\link{predict.dglars}}}
and \dQuote{\code{\link{phihat}}}. By default, this function returns the sequence of
the penalized coefficients and the sequence of the penalized estimate of the dispersion
parameter \eqn{\phi}{phi}; by using the argument \dQuote{\code{\dots}}, the user can obtain
the prediction of the estimates for any value of the tuning parameter (see the example below).
}
\value{
\code{coef.dglars} returns a named list with component:
\item{beta}{the sequence of the penalized estimates of the regression coefficients;}
\item{phi}{the penalized estimates of the dispersion parameter;}
\item{g}{the vector of the values of the tuning parameter.}
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{predict.dglars}} and \code{\link{phihat}}.
}
\examples{
###########################
# Logistic regression model
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
b <- 1:2
eta <- b[1] + X[, 1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(y ~ X, family = binomial)
coef(fit)
coef(fit, g = seq(4, 0.5, length = 10))

###########################
# Gamma family
n <- 100
p <- 10
X <- matrix(abs(rnorm(n * p)), n, p)
b <- 1:2
eta <- b[1] + X[, 1] * b[2]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))
coef(fit, type = "pearson")
coef(fit, type = "deviance")
coef(fit, type = "mle")
}
\keyword{models}
\keyword{regression}
