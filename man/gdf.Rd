\name{gdf}
\alias{gdf}
\alias{print.gof_dglars}
\title{
Estimate the Generalized Degrees-of-Freedom
}
\description{
\code{gdf} returns to estimate of the generalized degrees-of-freedom.
}
\usage{
gdf(object)
}
\arguments{
\item{object}{fitted \code{dglars} object.}
}
\details{
For a general nonlinear modelling procedure, a rigorous definition of 
degrees-of-freedom is obtained using the covariance penalty theory 
(Efron, 2004). This theory was used in Augugliaro et al. (2013) to
define a measure of model complexity for the dgLARS method, called
\dQuote{generalized degrees-of-freedom}. The \code{gdf} function
implements the estimator proposed in Augugliaro et al. (2013).
}
\value{
\code{gdf} returns a vector of length \code{np} with the generalized degrees-of-freedom.}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2014)
\emph{dglars: An R Package to Estimate Sparse Generalized Linear Models}, 
\emph{Journal of Statistical Software}, Vol 59(8), 1-40. \url{http://www.jstatsoft.org/v59/i08/}.

Augugliaro L., Mineo A.M. and Wit E.C. (2013)
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, 
\emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.

Efron B. (2004)
\emph{The estimation of prediction error: covariance penalties and cross-validation}, 
\emph{Journal of the American Statistical Association}, Vol. 99(467), 619-632.
}
\author{Luigi Augugliaro and Hassan Pazira\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}}, \code{\link{AIC.dglars}}, \code{\link{BIC.dglars}} and \code{\link{summary.dglars}}.
}
\examples{
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(y ~ X, binomial)
gdf(fit)
}
\keyword{models}
\keyword{regression}
