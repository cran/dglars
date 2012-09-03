\name{print.dglars}
\alias{print.dglars}
\title{
Print a dglars object
}
\description{
Print information about the sequence of models identified by the dgLARS method.
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
The call that produced the object \code{x} is printed, followed by a five-column \code{data.frame} with columns \code{Sequence}, \code{g}, \code{Dev}, \code{\%Dev} 
and \code{df}. The column named \code{Sequence} gives us information on how is changed the active set along the path. The column named \code{g} shows the 
sequence of the \eqn{\gamma} values used to compute the solution curve, while the columns \code{Dev} and \code{\%Dev} show the corresponding deviance and the 
fraction of explained deviance, respectively. Finally the \code{df} column shows the number of nonzero coefficients. The last part gives information about the algorithm 
and the method used to compute the solution curve. The code about the convergence of the used algorithm is also showed.
}
\value{The \code{data.frame} above is silently returned.}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} function.
}
\examples{
######################
# Logistic regression model
#
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit <- dglars(X, y, family = "binomial")
fit
}
\keyword{models}
\keyword{regression}
