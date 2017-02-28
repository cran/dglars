\name{summary.dglars}
\alias{summary.dglars}
\title{Summaryzing dgLARS Fits}
\description{
Summary method for an object with class 'dglars'.
}
\usage{
\method{summary}{dglars}(object, type = c("AIC", "BIC"), 
digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{object}{fitted \code{dglars} object.}
\item{type}{a description of the used measure of goodness-of-fit,
i.e., the Akaike Information Criterion (\code{type = "AIC"}) or the 
Bayesian Information Criterion (\code{type = "BIC"}). See below for 
the description of how to define more general measures of 
goodness-of-fit.}
\item{digits}{significant digits in printout.}
\item{\dots}{additional arguments passed to \dQuote{\code{AIC.dglars}} 
or \dQuote{\code{BIC.dglars}}. See below for more details.}
}
\details{
\code{summary.dglars} gives information about the sequence of models 
estimtaed by dgLARS method.

To select the best fitted model, summary method uses a measure of
goodness-of-fit (GoF) defined as follows: 
\deqn{-2 \mbox{log-likelihood} + k \mbox{comp},}{-2 * log-likelihood + k * comp,}
where \dQuote{comp} represents the term used to measure the complexity
of the fitted models, and \eqn{k} is the 'weight' of the complexity 
in the previous formula. This quantity is computed using the functions
\dQuote{\code{AIC.dglars}} or \dQuote{\code{BIC.dglars}}.

By default, \code{summary.dglars} function computes the AIC criterion, 
but the user can use \dQuote{\code{dots}} to pass to the function 
\code{AIC.dglars} the additional arguments needed to compute a more
general measure of goodness-of-fit, i.e, \dQuote{\code{g}}, 
\dQuote{\code{phi}}, \dQuote{\code{k}} and \dQuote{\code{complexity}}
(see \dQuote{\code{\link{AIC.dglars}}} for the description of these 
arguments). Below we give some examples on how to use these additional
arguments. 

The output of the summary method is divided in two sections. 

The first section shows the call that produced the object \code{x} 
followed by a \code{data.frame}. The column named \dQuote{\code{Sequence}}
gives information on how is changed the active set along the path. 
The column \dQuote{\code{g}} shows the sequence of the \eqn{\gamma}{g}-values
used to compute the solution curve, while the column \dQuote{\code{\%Dev}}
shows the the fraction of explained deviance. The remaining columns 
show the complexity measure, the used measure of goodness-of-fit 
and the corresponding ranking of the fitted models. 

The second section shows the details of the selected
model, i.e. family and link function used to specify the generalized 
linear model, the penalized estimate of the coefficient vector,
the value of the tuning parameter, the null and residual deviance, and
finaly the value of the used measure of goodness-of-fit. Information 
about the method and the algorithm used to compute the solution curve 
is also provided.
}
\value{
\code{summary.dglars} function silently returns a named list with components:
\item{table}{a \code{data.frame} with the information about the sequence 
of model fitted by \code{dglars} function;}
\item{formula.gof}{if the model is specified by the formula in \code{dglars}, 
then the formula of the selected model is reported;}
\item{b.gof}{the estimates of the coefficients of the selected model;}
\item{phi.gof}{the estimate of the dispersion parameter of the selected model;}
\item{nulldev}{the null residual deviance;}
\item{resdev.gof}{the residual deviance of the selected model;}
\item{type}{a description of the measure of goodness-of-fit used to select the
model;}
\item{k}{the \sQuote{weight} used to compute the measure of goodness-of-fit;}
\item{complexity}{a description of the method used to measure the complexity
of the fitted models;}
\item{phi}{a description of the method used the estimate the dispersion
parameter.}
}
\author{Luigi Augugliaro and Hassan Pazira\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}}, \code{\link{AIC.dglars}}, \code{\link{BIC.dglars}} and \code{\link{gdf}} functions.
}
\examples{
###########################
# Logistic regression model
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[, 1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(y ~ X, family = binomial)

summary(fit, type = "AIC")
summary(fit, type = "AIC", k = 0.1)
summary(fit, type = "AIC", complexity = "gdf")
summary(fit, type = "BIC", complexity = "df")
summary(fit, type = "BIC", complexity = "gdf")

############################
# y ~ Gamma
n <- 100
p <- 10
X <- matrix(abs(rnorm(n * p)), n, p)
eta <- 1 + 2 * X[, 1]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))

summary(fit, phi = "pearson")
summary(fit, phi = "deviance")
summary(fit, phi = "mle")
}
\keyword{models}
\keyword{regression}
