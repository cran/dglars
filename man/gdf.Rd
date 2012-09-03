\name{gdf}
\alias{gdf}
\title{
Compute the generalized degrees-of-freedom.
}
\description{
\code{gdf} is used to compute the first order estimator of the generalized degrees-of-freedom for logistic regression model.
}
\usage{
gdf(X, y, beta = NULL, control = list())
}
\arguments{
  \item{X}{Design matrix of dimension \eqn{n\times p}.}
  \item{y}{Response vector.}
  \item{beta}{The \eqn{(p+1)\times\code{np}} matrix corresponding to the dgLARS solution curve.}
  \item{control}{The list of control parameters passed to \code{dglars} function.}
}
\details{
For a general nonlinear modelling procedure, a more rigorous definition of degrees-of-freedom is obtained using the covariance penalty theory (Efron, 2004). 
When we work with a logistic regression model defined in a low dimensional setting, the \code{gdf} function can be used to compute the first order estimator 
proposed in Augugliaro et al. (2013). How to define an estimator of the generalized degrees-of-freedom in a high-dimensional setting is still an open 
question. Simulation studies seem to show that for a Poisson regression model the number of nonzero coefficients can be considered a satisfying approximation 
to the generalized degrees-of-freedom, for this reason the first order estimator is not implemented for this model.
}
\value{
\code{gdf} returns a vector of length \code{np} with the generalized degrees-of-freedom.}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2013)
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, \emph{Journal of the Royal Statistical Society. Series B.} Vol 75(3), 471-498.\cr\cr
Efron B. (2004)
\emph{The estimation of prediction error: covariance penalties and cross-validation}, \emph{Journal of the American Statistical Association}, Vol. 99(467), 619-632.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} function.
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
fit <- dglars(X, y, family = "binomial")
beta <- coef(fit)
gdf(X, y, beta)
}
\keyword{models}
\keyword{regression}
