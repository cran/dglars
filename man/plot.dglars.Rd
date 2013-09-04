\name{plot.dglars}
\alias{plot.dglars}
\title{Plot from a dglars object}
\description{
Produces different plots to study the sequence of models identified by the dgLARS method.
}
\usage{
\method{plot}{dglars}(x, k = c("BIC","AIC"), complexity = c("df","gdf"), 
	g.gof = NULL, ...)
}
\arguments{
  \item{x}{fitted \code{dglars} object.}
  \item{k}{character/numeric argument used to specify the 'weight' of the complexity part in the measure of goodness-of-fit used to select the best model (see below for more details). 
  Default is \code{k = "BIC"};}
  \item{complexity}{measure used to define the complexity of a model. If \code{complexity = "df"} then the complexity is defined as the number of nonzero coefficients. For logistic 
  regression model it is possible to use the generalized degrees of freedom to measure the complexity of a model (\code{complexity = "gdf"}). Default is \code{complexity = "df"};}
  \item{g.gof}{a user specified value of the tuning parameter;}
  \item{\dots}{additional graphical parameters.}
}
\details{
\code{plot.dglars} method produces different plots to study the sequence of models identified by the dgLARS method. 

The first plot show the goodness-of-fit (Gof) criterion  \deqn{Dev + k \times complexity,}{Dev + k  complexity,} as functions of the used \eqn{\gamma} values. In the GoF formula the \eqn{Dev} component is the residual deviance,  \eqn{complexity} is the term used to measure the complexity of the fitted model while \eqn{k} is the term used to 'weight' the complexity part in the GoF formula. By default, \code{plot.dlgars} uses the BIC criterion as measure of goodness-of-fit, i.e \code{k = "BIC"} and \code{complexity = "df"}, but the AIC criterion can be easily computed setting \code{k = "AIC"} and \code{complexity = "df"}. See also \code{summary.dglars} for more details about the arguments \code{k} and \code{complexity}. The value of the tuning parameter corresponding to the minimum of the GoF curve is labeled by a vertical dashed red line, while the \eqn{\gamma} values at which corresponds a change in the active set are labeled by  vertical dashed gray lines. Next plot shows the coefficient profile plot while, if the control parameter \code{algorithm} is equal to \code{"pc"}, the last plot shows the Rao's score test statistics as functions of \eqn{\gamma}. Also in these two plots, the optimal \eqn{\gamma} and the \eqn{\gamma} values 
at which corresponds a change in the active set are labeled by vertical dashed lines, red and gray respectively.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} and \code{\link{summary.dglars}} function.
}
\examples{
###########################
# Logistic regression model

set.seed(123)

n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit_pc <- dglars.fit(X, y, family = "binomial")
par(mfrow = c(2, 3))
plot(fit_pc, k = "BIC", complexity = "df")
plot(fit_pc, k = "AIC", complexity = "df")

fit_ccd <- dglars.fit(X, y, family = "binomial", 
control = list(algorithm = "ccd"))
par(mfrow = c(2, 2))
plot(fit_ccd, k = "BIC", complexity = "df")
plot(fit_ccd, k = "AIC", complexity = "df")
}
\keyword{models}
\keyword{regression}
