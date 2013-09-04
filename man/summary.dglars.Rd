\name{summary.dglars}
\alias{summary.dglars}
\title{Print the summary for a dglars object}
\description{
Summary method for class 'dglars'.
}
\usage{
\method{summary}{dglars}(object, k = c("BIC", "AIC"), complexity = c("df", "gdf"), 
digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object}{fitted \code{dglars} object}
  \item{k}{character/numeric argument used to specify the 'weight' of the complexity part in the measure of goodness-of-fit used to select the best model (see below for more details). 
  Default is \code{k = "BIC"};}
  \item{complexity}{measure used to define the complexity of a model. If \code{complexity = "df"} then the complexity is defined as the number of nonzero coefficients. For logistic 
  regression model it is possible to use the generalized degrees of freedom to measure the complexity of a model (\code{complexity = "gdf"}). Default is \code{complexity = "df"};}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
\code{summary.dglars} gives us information about the sequence of models identified by dgLARS method. 

To select the best model the summary method uses a measure of goodness-of-fit (GoF) defined as follows: \deqn{Dev + k \times complexity,}{Dev + k  complexity,} 
where \eqn{Dev} is the residual deviance, \eqn{complexity} is the term used to measure the complexity of the fitted model and \eqn{k} is the term used to 'weight'  
the complexity part in the GoF formula.  By default the summary method computed the BIC criterion to select the best model, which means that the argument \code{complexity} is 
equal \code{"df"} (the number of nonzero coefficients) and the argument \code{k} is equal to \code{"BIC"}. The AIC criterion can be easily computed setting \code{k = "AIC"} and 
\code{complexity = "df"}. The user can define other measures of goodness-of-fit, specifying \code{k} as any non-negative integer and, for example, using the generalized degrees-of-
freedom to measure the complexity of the fitted model (see the \code{gdf} function for more details). 

The output of the summary method is divided in two sections. The first section shows the call that produced the object \code{x} followed by a \code{data.frame}. 
The column named \code{Sequence} gives us information on how is changed the active set along the path. The column named \code{g} shows the sequence of the \eqn{\gamma} 
values used to compute the solution curve, while the column \code{Dev} shows the corresponding residual deviance. The remaining columns show the complexity measure, the 
used measure of goodness-of-fit and the corresponding ranking of the models. Finally, the second section show the estimated coefficients of the best models identified by the used 
GoF criterion. Like for the \code{print.dglars} method, information about the method, the algorithm and the corresponding convergence are also provided.
}
\value{
A list with components \code{table} and \code{b.gof} is silently returned. The \code{table} component is the \code{data.frame} previously described while the
component \code{b.gof} is the vector of the estimated coefficients corresponding to the best models identified using the GoF criterion.
Alternatively, the user can use the internal function \code{make_summary_table} to obtain this list (see the example below).
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{dglars}} and \code{\link{gdf}} functions.
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
fit <- dglars.fit(X, y, family = "binomial")

## the BIC criterion is used to select the best model
summary(fit, k= "BIC",complexity = "df")
## the AIC criterion is used to select the best model
summary(fit, k= "AIC",complexity = "df")

tbl <- make_summary_table(fit, k=log(n), complexity = "df")
str(tbl)
summary(fit, complexity = "gdf")
tbl <- make_summary_table(fit, k=log(n), complexity = "gdf")
str(tbl)

##########################
# Poisson regression model

n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)
fit <- dglars.fit(X, y, family = "poisson")

## the BIC criterion is used to select the best model
summary(fit, k= "BIC",complexity = "df")
## the AIC criterion is used to select the best model
summary(fit, k= "AIC",complexity = "df")

tbl <- make_summary_table(fit, k=log(n), complexity = "df")
str(tbl)
summary(fit, k=log(n), complexity = "gdf")
}
\keyword{models}
\keyword{regression}
