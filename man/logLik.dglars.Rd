\name{logLik.dglars}
\alias{logLik}
\alias{logLik.dglars}
\alias{print.loglik_dglars}

\title{Extract Log-Likelihood}

\description{logLik method for an object with class \sQuote{dglars}.}

\usage{
\method{logLik}{dglars}(object, phi = c("pearson", "deviance", "mle", "grcv"),
	g = NULL, \dots)
}

\arguments{
\item{object}{any fitted \code{dglars} object from which the log-likelihood values can be extracted.}
\item{phi}{a description of the estimator used to estimate the dispersion pamater (see below for more details).}
\item{g}{vector of values of the tuning parameter.}
\item{\dots}{further arguments passed to the function \code{\link{grcv}}.}
}

\details{
\code{logLik.dglars} returns the sequence of the log-likelihood values of the models fitted by \dQuote{\code{dglars}}. For the binomial and Poisson family, the dispersion parameter is assumed known and equal to one while for the other families the dispersion parameter is estimated using the method specified by the argument \dQuote{\code{phi}} (see \code{\link{phihat}} for more details). The optional argument \code{g} is used to specified the values of the tuning parameter; if not specified (default), the log-likelihood values are computed for the sequence of models storage in \code{object} otherwise \code{\link{predict.dglars}} is used to compute the estimate of the parameters needed to evaluate the log-likelihood function (see the example below).
}

\value{
\code{logLik.dglars} returns an object of class \dQuote{\code{loglik_dglars}}. This is a named list containing the following components:
\item{loglik}{the log-likelihood of the sequence of models fitted by dglars method.}
\item{df}{the numbers of non-zero estimates corresponding to the used \eqn{\gamma}{g}-values.}
\item{object}{the fitted \code{dglars} object.}
\item{g}{the sequence of \eqn{\gamma}{g}-values.}
\item{phi}{a description of the estimator used to estimate the dispersion pamater.}
\item{phih}{the sequence of estimated dispersion parameter.}
}

\seealso{
\code{\link{dglars}}, \code{\link{phihat}}, \code{\link{predict.dglars}}, \code{\link{AIC.dglars}} and \code{\link{BIC.dglars}}.
}

\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}

\examples{
#################################
# y ~ Poisson

library(dglars)
set.seed(123)
n <- 100
p <- 5
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + X[, 1]  + X[, 2]
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)
out <- dglars(y ~ X, poisson)
logLik(out)
logLik(out, g = seq(2, 0.5, by = -0.1))

############################
# y ~ Gamma

n <- 100
p <- 5
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[, 1L]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
out <- dglars(y ~ X, Gamma("log"))

# generalized Pearson statistic
logLik(out, phi = "pearson")

# deviance estimator
logLik(out, phi = "deviance")

# mle estimator
logLik(out, phi = "mle")

# grcv estimator
logLik(out, phi = "grcv")
}
\keyword{models}
