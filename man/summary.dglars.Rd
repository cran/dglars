\name{summary.dglars}
\alias{summary.dglars}
\title{Print the summary for a dglars object}
\description{
Summary method for class 'dglars'.
}
\usage{
\method{summary}{dglars}(object, complexity = c("df", "gdf"),
	digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{fitted \code{dglars} object}
  \item{complexity}{measure used to define the complexity of a model. If \code{complexity = "df"} then the complexity is defined as the number of nonzero coefficients. For logistic 
  regression model it is possible to use the generalized degrees of freedom (\code{complexity = "gdf"}) to measure the complexity of a model. Default is \code{complexity = "df"};}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
\code{summary.dglars} gives us information about the sequence of models identified by dgLARS method. The output is divided in three sections. The first section shows the
call that produced the object \code{x} followed by a \code{data.frame}. The column named \code{Sequence} gives us information on how is changed the active set along the path. 
The column named \code{g} shows the sequence of the \eqn{\gamma} values used to compute the solution curve, while the column \code{Dev} shows the corresponding 
deviance. The column named \code{Complexity} shows the type of measure used to define the complexity of a model. The columns named \code{AIC} and \code{Rank.AIC} 
show the Akaike Information Criterion (AIC) and the ranking of the models obtained using these measures of goodness-of-fit. Similarly, the columns named \code{BIC} and 
\code{Rank.BIC} give information on the sequence of models based on the Bayesian Information Criterion (BIC). The best models identified by these two criteria are labeled 
by an arrow. Finally, the second and the third section show the estimated coefficients of the best models identified by the AIC and the BIC criteria. Like for the \code{print.dglars} 
method, information about the method, the algorithm and the corresponding convergence are also provided.
}
\value{
A list with components \code{table}, \code{b.aic} and \code{b.bic} is silently returned. The \code{table} component is the \code{data.frame} previously described while the
components \code{b.aic} and \code{b.bic} are the two vectors of the estimated coefficients corresponding to the best models identified using the AIC and the BIC criterion.
Alternatively, the user can use the internal function \code{make_summary_table} to obtain this list (see example below).
}
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
summary(fit, complexity = "df")
tbl <- make_summary_table(fit, complexity = "df")
str(tbl)
summary(fit, complexity = "gdf")
tbl <- make_summary_table(fit, complexity = "gdf")
str(tbl)

######################
# Poisson regression model
#
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)
fit <- dglars(X, y, family = "poisson")
summary(fit, complexity = "df")
tbl <- make_summary_table(fit, complexity = "df")
str(tbl)
summary(fit, complexity = "gdf")
}
\keyword{models}
\keyword{regression}
