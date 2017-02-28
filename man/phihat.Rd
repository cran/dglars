\name{phihat}
\alias{phihat}
\title{
Estimate the Dispersion Parameter
}
\description{
\code{phihat} returns the estimates of the dispersion parameter.
}
\usage{
phihat(object, type = c("pearson", "deviance", "mle"), g = NULL)
}
\arguments{
\item{object}{fitted \code{dglars} object.}
\item{type}{a description of the used estimator.}
\item{g}{vector of values of the tuning parameter.}
}
\details{
\code{phihat} implements three different estimators of the dispersion parameter, i.e,
the generalized Pearson statistic (\code{type = "pearson"}), the deviance estimator
(\code{type = "deviance"}) and the maximum likelihood estimator (\code{type = "mle"}).
For regression models with Gamma family, the maximum likelihood estimator of the
dispersion parameter is computed using the approximation proposed in Cordeiro et al.
(1997).
}
\value{
\code{phihat} returns a vector with the estimates of the dispersion parameter.}
\references{

Gauss M. Cordeiro and Peter McCullagh (1991)
\emph{Bias Correction in Generalized Linear Models}, 
\emph{Journal of the Royal Statistical Society. Series B.}, Vol 53(3), 629--643.

Bent Jorgensen (1997)
\emph{The Theory of Dispersion Models},
Chapman \& Hall, Great Britain.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{coef.dglars}}, \code{\link{logLik.dglars}}, \code{\link{AIC.dglars}} and \code{\link{BIC.dglars}}.
}
\examples{
############################
# y ~ Gamma
set.seed(321)
n <- 100
p <- 100
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[,1]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))
g <- seq(range(fit$g)[1], range(fit$g)[2], length = 10)
# generalized Pearson statistic
phihat(fit, type = "pearson")
phihat(fit, type = "pearson", g = g)
# deviance estimator
phihat(fit, type = "deviance")
phihat(fit, type = "deviance", g = g)
# mle
phihat(fit, type = "mle")
phihat(fit, type = "mle", g = g)
}
\keyword{models}
\keyword{regression}
