\name{print.dglars}
\alias{print.dglars}
\title{
Printing a dgLARS Object
}
\description{
Print information about the sequence of models estimated by dgLARS method.
}
\usage{
\method{print}{dglars}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{x}{fitted \code{dglars} object}
\item{digits}{significant digits in printout}
\item{\dots}{additional print arguments}
}
\details{
The call that produced the object \code{x} is printed, followed by a 
five-column \code{data.frame} with columns \dQuote{\code{Sequence}}, 
\dQuote{\code{g}}, \dQuote{\code{Dev}}, \dQuote{\code{\%Dev}} and 
\dQuote{\code{n. non zero}}. The column named \dQuote{\code{Sequence}} 
gives information on how is changed the active set along the path. The 
column \dQuote{\code{g}} shows the sequence of \eqn{\gamma}{g} 
values used to compute the solution curve, while the columns \dQuote{\code{Dev}} 
and \dQuote{\code{\%Dev}} show the corresponding deviance and the fraction 
of explained deviance, respectively. Finally the \dQuote{\code{n. non zero}}
column shows the number of nonzero coefficients. The last part gives information 
about the algorithm and the method used to compute the solution curve. 
The code about the convergence of the used algorithm is also showed.
}
\value{The \code{data.frame} above is silently returned.}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} function.
}
\examples{
#############################
# y ~ Binomial
set.seed(123)
n <- 100
p <- 100
X <- matrix(rnorm(n * p), n, p)
eta <- 1 + 2 * X[,1]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(y ~ X, family = binomial)
fit

# adaptive dglars method
b_wght <- coef(fit)$beta[, 20]
fit <- dglars(y ~ X, family = binomial, b_wght = b_wght) 
fit 

# the first three coefficients are not penalized
fit <- dglars(y ~ X, family = binomial, unpenalized = 1:3) 
fit 

# 'probit' link function
fit <- dglars(y ~ X, family = binomial("probit"))
fit
}
\keyword{models}
\keyword{regression}
