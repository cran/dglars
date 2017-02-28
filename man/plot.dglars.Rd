\name{plot.dglars}
\alias{plot.dglars}
\title{Plot from a dglars Object}
\description{
Produces plots to study the sequence of models identified by dgLARS method.
}
\usage{
\method{plot}{dglars}(x, type = c("both", "AIC", "BIC"), ...)
}
\arguments{
\item{x}{fitted \code{dglars} object.}

\item{type}{ a description of the measure of goodness-of-fit used to compare
the sequence of models fitted by \code{dglars} or \code{dglars.fit}. See
below for a more details.}
  	
\item{\dots}{further arguments passed to the functions \code{AIC.dglars} or 
  	\code{BIC.dglars}.}
}
\details{
\code{plot.dglars} method produces different plots to study the sequence of 
models fitted by dgLARS method. 

First plot gives information about the goodness-of-fit of the sequence of models 
fitted by dgLARS method. The user can plot the sequence of AIC (\code{type = "AIC"})
or BIC values (\code{type = "BIC"}). By default, AIC and BIC values are shown on the 
same plot (\code{type = "both"}). More general measures of goodness-of-fit can be 
specified by using the argument \dQuote{\dots} to pass futher arguments to function 
\code{\link{AIC.dglars}} (see the examples below). The value of the tuning parameter 
corresponding to the minimum of the used measure of goodness-of-fit is indentified by
a vertical dashed red line, while the \eqn{\gamma}{g} values at which corresponds a 
change in the active set are labeled by vertical dashed gray lines. Second plot shows 
the coefficient profile plot; if the predictor-corrector algorithm is used to fit the
model, the third plot shows the Rao's score test statistics as function of \eqn{\gamma}{g}.
}
\author{Luigi Augugliaro and Hassan Pazira\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}}, \code{\link{summary.dglars}} and \code{\link{AIC.dglars}}.
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
fit <- dglars.fit(X, y, family = binomial)
plot(fit) 
plot(fit, type = "AIC")
plot(fit, type = "BIC")
plot(fit, type = "AIC", k = 5)
plot(fit, type = "AIC", complexity = "gdf")
}
\keyword{models}
\keyword{regression}
