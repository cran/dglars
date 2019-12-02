\name{AIC.dglars}
\alias{AIC}
\alias{AIC.dglars}
\alias{BIC}
\alias{BIC.dglars}

\title{Akaike's An Information Criterion}

\description{
\code{AIC.dglars} is used to compute the Akaike's \sQuote{An Information Criterion} for the sequence of models estimated by \dQuote{\code{dglars}}.
}

\usage{
\method{AIC}{dglars}(object, phi = c("pearson", "deviance", "mle", "grcv"),
k = 2, complexity = c("df", "gdf"), g = NULL, ...)

\method{BIC}{dglars}(object, \dots)
}

\arguments{
\item{object}{a fitted \code{dglars} object.}
\item{phi}{a description of the estimator of the dispersion parameter (see below for more details).}
\item{k}{non negative value used to weight the complexity of the fitted dglars model (see below for more details).}
\item{complexity}{argument used to specify the method to measure the complexity of a fitted dglars model, i.e. the number of non-zero estimates (\code{complexity = "df"}) of the generalized degrees-of-freedom (\code{complexity = "gdf"}); see below for more details.}
\item{g}{vector of values of the tuning parameter.}
\item{\dots}{further arguments passed to the function \code{link{grcv}}.}
}

\details{
The values returned by \code{AIC.dglars} are computed according to the following formula of a generic measure of Goodness-of-Fit (GoF): \deqn{-2 \mbox{log-likelihood} + k \mbox{comp},}{-2 * log-likelihood + k * comp,} where \dQuote{comp} represents the term used to measure the complexity of the fitted model, and \eqn{k} is the \sQuote{weight} of the complexity in the previous formula.

For binomial and Poisson family, the log-likelihood function is evaluated assuming that the dispersione parameter is known and equal to one while for the remaining families the dispersion parameter is estimated by the method specified by \code{phi} (see \code{\link{phihat}} for more details).
  
According to the results given in Augugliaro et. al. (2013), the complexity of a model fitted by dglars method can be measured by the classical notion of \sQuote{Degrees-of-Freedom} (\code{complexity = "df"}), i.e., the number of non-zero estimated, or by the notion of \sQuote{Generalized Degrees-of-Freedom} (\code{complexity = "gdf"}).

By the previous formula, it is easy to see that the standard AIC-values are obtained setting \code{k = 2} and \code{complexity = "df"} (default values for the function \code{AIC.dglars}) while the so-called BIC-values (Schwarz's Bayesian criterion) are obtained setting \code{k = log(n)}, where \eqn{n}{n} denotes the sample size, and \code{complexity = "df"} (default values for the function \code{BIC.dglars}).

The optional argument \code{g} is used to specify the values of the tuning parameter; if not specified (default), the values of the measure of goodness-of-fit are computed for the sequence of models storage in \code{object} otherwise \code{\link{predict.dglars}} is used to compute the estimate of the parameters needed to evaluate the log-likelihood function (see the example below).
}

\value{
\code{AIC.dglars} and \code{BIC.dglars} return a named list with class \dQuote{\code{gof_dglars}} and components:
\item{val}{the sequence of AIC/BIC-values;}
\item{g}{the sequence of \eqn{\gamma}{g}-values;}
\item{loglik}{the sequence of log-likelihood values used to compute the AIC or BIC;}
\item{k}{the non negative value used to weight the complexity of the fitted dglars model;}
\item{comp}{the measures of model complexity used to compute the measure of goodness-of-fit. It is equal to \code{npar} when code{complexity = "df"};}
\item{npar}{the seqeunce of the number of non-zero estimates}
\item{phi}{a description of the estimator used to estimate the dispersion pamater;}
\item{phih}{the vector of penalized estimate of the dispersion parameter used to evaluate the log-likelihood function;}
\item{complexity}{character specifying the method to measure the complexity of a fitted dglars model;}
\item{object}{the fitted \code{dglars} object;}
\item{type}{character specifying the type of used measure-of-goodness of fit, i.e., AIC, BIC or GoF.}

In order to summarize the information about the AIC-valuse, a \code{print} method is available for an object with class \dQuote{\code{gof_dglars}}.
}

\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2013) <doi:10.1111/rssb.12000>
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, 
\emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.

Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986, ISBN:978-90-277-2253-9)
\emph{Akaike Information Criterion Statistics}.
KTK Scientific Publishers, 1986.
}

\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}

\seealso{
\code{\link{logLik.dglars}}, \code{\link{predict.dglars}}, \code{\link{dglars}} and \code{\link{summary.dglars}}.
}

\examples{
#################################
# y ~ Pois

library("dglars")
set.seed(123)
n <- 100
p <- 5
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + X[, 1]  + X[, 2]
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)
out <- dglars(y ~ X, poisson)
out
AIC(out)
AIC(out, g = seq(2, 1, by = -0.1))
AIC(out, complexity = "gdf")
AIC(out, k = log(n)) #BIC-values
BIC(out)

#################################
# y ~ Gamma

n <- 100
p <- 50
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[, 1L]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
out <- dglars(y ~ X, Gamma("log"))

AIC(out, phi = "pearson")
AIC(out, phi = "deviance")
AIC(out, phi = "mle")
AIC(out, phi = "grcv")
}
\keyword{models}
