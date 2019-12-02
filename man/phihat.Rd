\name{phihat}
\alias{phihat}

\title{Estimate the Dispersion Parameter}

\description{
\code{phihat} returns the estimates of the dispersion parameter.
}

\usage{
phihat(object, type = c("pearson", "deviance", "mle", "grcv"), g = NULL, \dots)
}

\arguments{
\item{object}{fitted \code{dglars} object.}
\item{type}{a description of the used estimator.}
\item{g}{vector of values of the tuning parameter.}
\item{\dots}{further arguments passed to the function \code{\link{grcv}}.}
}

\details{
\code{phihat} implements four different estimators of the dispersion parameter, i.e, the generalized Pearson statistic (\code{type = "pearson"}), the deviance estimator (\code{type = "deviance"}), the maximum likelihood estimator (\code{type = "mle"}) and general refitted cross-Validation estimator (\code{type = "grcv"}) proposed in Pazira et al. (2018). For regression models with Gamma family, the maximum likelihood estimator of the dispersion parameter is computed using the approximation proposed in Cordeiro et al. (1997).
}

\value{\code{phihat} returns a vector with the estimates of the dispersion parameter.}

\references{
Cordeiro G. M. and McCullagh P. (1991) <doi:10.2307/2345592>
Bias Correction in Generalized Linear Models,
\emph{Journal of the Royal Statistical Society. Series B.}, Vol 53(3), 629--643.

Jorgensen B. (1997, ISBN:0412997188)
\emph{The Theory of Dispersion Models},
Chapman \& Hall, Great Britain.

Pazira H., Augugliaro L. and Wit E.C. (2018) <doi:10.1007/s11222-017-9761-7>
Extended differential-geometric LARS for high-dimensional GLMs with general dispersion parameter,
\emph{Statistics and Computing}, Vol 28(4), 753-774.
}

\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}

\seealso{
\code{\link{grcv}}, \code{\link{coef.dglars}}, \code{\link{logLik.dglars}}, \code{\link{AIC.dglars}} and \code{\link{BIC.dglars}}.
}

\examples{
############################
# y ~ Gamma

library("dglars")
set.seed(321)
n <- 100
p <- 50
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[, 1L]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))
g <- seq(range(fit$g)[1L], range(fit$g)[2L], length = 10)

# generalized Pearson statistic
phihat(fit, type = "pearson")
phihat(fit, type = "pearson", g = g)

# deviance estimator
phihat(fit, type = "deviance")
phihat(fit, type = "deviance", g = g)

# mle
phihat(fit, type = "mle")
phihat(fit, type = "mle", g = g)

# grcv
phihat(fit, type = "grcv")
phihat(fit, type = "grcv", g = g)
}
\keyword{models}
\keyword{regression}
