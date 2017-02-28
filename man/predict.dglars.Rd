\name{predict.dglars}
\alias{predict.dglars}
\title{
Predict Method for dgLARS Fits.
}
\description{
\code{predict.dglars} is used to obtain general predictions from a \code{dglars} object.
}
\usage{
\method{predict}{dglars}(object, xnew, ynew, g = NULL, 
type = c("coefficients", "nnonzero", "predictors", "eta", 
"mu", "probability", "class", "deviance"), \dots)
}
\arguments{
\item{object}{fitted \code{dglars} object.}
\item{xnew}{matrix of new values of the predictors at which predictions are to be made. 
This argument is not used for \code{type} \dQuote{\code{coefficients}}, \dQuote{\code{nnonzero}}
and \dQuote{\code{predictors}}.}
\item{ynew}{vector of new values of the responce variable. This argument is used 
only when \code{type} ie equal to \dQuote{\code{deviance}}.}
\item{g}{value(s) of the tuning parameter \eqn{\gamma}{g} at which the predictions are required. 
By default, the predictions are made using the sequence of \eqn{\gamma}{g} values storaged in \code{dglars}.}
\item{type}{type of prediction required; see below for more details.}
\item{\dots}{additional argument used to ensure the compatibility with 
the generic method function \dQuote{\code{\link{predict}}}.}
}
\value{
The object returned by \code{predict.dglars} depends on \code{type} argument:

\item{\code{coefficients}:}{a named list with components \dQuote{\code{beta}}, i.e., the matrix corresponding 
to the dgLARS solution curve, and \dQuote{\code{phi}}, i.e., the sequence of Pearson estimates of the 
dispersion parameter;}

\item{\code{nnonzero}:}{the number of nonzero estimates;}

\item{\code{predictors}:}{a named list; each component is a vector containing the indices of the variables 
that are in the active set;}

\item{\code{eta}:}{a matrix with the linear preditors. If \code{xnew} is not specified then the linear predictors
are computed using \code{object$X};}

\item{\code{mu}:}{a matrix with the fitted expeted values, obtained by transforming the linear predictor by the
inverse of the link function. For models with \sQuote{binomial} family, canonical link function (\sQuote{\code{logit}})
and response a vector with elements 0 (failure) or 1 (success), \code{type = "mu"} and \code{type = "probability"} give
the same result. If \code{xnew} is not spefified then these values are computed using \code{object$X};}

\item{\code{probability}:}{available only for \sQuote{\code{binomial}} family. In this case \code{predict.dglars}
returns a matrix with the fitted probabilities; furthermore, if the model is specified by using the canonical 
link function (\sQuote{\code{logit}}) and response a vector with entries 0/1, \code{type = "mu"} and \code{type = "probability"} give the same result. If \code{xnew} is not specified then
these values are computed using \code{object$X};}

\item{\code{class}:}{available only for \sQuote{\code{binomial}} family. In this case \code{predict.dglars}
returns a matrix with the fitted class. If \code{xnew} is not specified then these values are computed using
\code{object$X} otherwise \code{xnew} is used to compute the fitted probabilities;}

\item{\code{deviance}:}{a vector with the scaled residual deviances.}
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}
\seealso{
\code{\link{dglars}} and \code{\link{coef.dglars}}.
}
\examples{
######################
# Logistic regression model
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
Xnew <- matrix(rnorm(n * p), n, p)
b <- 1:2
eta <- b[1] + X[, 1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars.fit(X, y, binomial)
coef(fit)
predict(fit, type = "coefficients")
g <- seq(3, 1, by = -0.1)
coef(fit, g = g)
predict(fit, type = "coefficients", g = g)
predict(fit, type = "nnonzero")
predict(fit, type = "nnonzero", g = g)
predict(fit, type = "predictors")
predict(fit, type = "predictors", g = g)
predict(fit, type = "eta", g = g)
predict(fit, type = "eta", g = g, xnew = Xnew)
predict(fit, type = "mu", g = g)
predict(fit, type = "mu", g = g, xnew = Xnew)
predict(fit, type = "probability", g = g)
predict(fit, type = "probability", g = g, xnew = Xnew)
predict(fit, type = "class", g = g)
predict(fit, type = "class", g = g, xnew = Xnew)

}
\keyword{models}
\keyword{regression}
