\name{grcv}
\alias{grcv}

\title{
General Refitted Cross-Validation Estimator
}

\description{
\code{grcv} computes the estimate of the dispersion parameter using the general refitted cross-validation method.
}

\usage{
grcv(object, type = c("BIC", "AIC"), nit = 10L, trace = FALSE,
     control = list(), \dots)
}

\arguments{
\item{object}{fitted \code{dglars} object.}
\item{type}{the measure of goodness-of-fit used in Step 2 to select the two set of variables (see section \bold{Description} for more details). Default is \code{type = BIC}.}
\item{control}{a list of control parameters passed to the function \code{\link{dglars}}.}
\item{nit}{integer specifying the number of times that the general refitted cross-validation method is repeated (see section \bold{Description} for more details). Default is \code{nit = 10L}.}
\item{trace}{flag used to print out information about the algorithm. Default is \code{trace = FALSE}.}
\item{\dots}{further arguments passed to the functions \code{\link{AIC.dglars}} or \code{\link{BIC.dglars}}.}
}

\details{
The general refitted cross-validation (grcv) estimator (Pazira et al., 2018) is an estimator of the dispersion parameter of the exponential family based on the following four stage procedure:
\tabular{cl}{
Step \tab Description \cr
1. \tab randomly split the data set \eqn{D = (y, X)}{D = (y, X)} into two even datasets, denoted by \eqn{D_1}{D_1} and \eqn{D_2}{D_2}.\cr
2. \tab fit dglars model to the dataset \eqn{D_1}{D_1} to select a set of variables \eqn{A_1}{A_1}.\cr
   \tab fit dglars model to the dataset \eqn{D_2}{D_2} to select a set of variables \eqn{A_2}{A_2}.\cr
3. \tab fit the glm model to the dataset \eqn{D_1}{D_1} using the variables that are in \eqn{A_2}{A_2}; then estimate the \cr
   \tab disporsion parameter using the Pearson method. Denote by \eqn{\hat{\phi}_1(A_2)}{hat{phi}_1(A_2)} the resulting estimate.\cr
   \tab fit the glm model to the dataset \eqn{D_2}{D_2} using the variables that are in \eqn{A_1}{A_1}; then estimate the \cr
   \tab disporsion parameter using the Pearson method. Denote by \eqn{\hat{\phi}_2(A_1)}{hat{phi}_2(A_1)} the resulting estimate.\cr
4. \tab estimate \eqn{\phi}{phi} using the following estimator: \eqn{\hat{\phi}_{grcv} = (\hat{\phi}_1(A_2) + \hat{\phi}_2(A_1)) / 2}{\hat{phi}_{grcv} = (hat{phi}_1(A_2) + hat{phi}_2(A_1)) / 2}.
}
In order to reduce the random variabilty due to the splitting of the dataset (Step 1), the previous procedure is repeated \sQuote{\code{nit}}-times; the median of the resulting estimates is used as final estimate of the dispersion parameter. In Step 3, the two sets of variables are selected using the \code{\link{AIC.dglars}} and \code{\link{BIC.dglars}}; in this step, the Pearson method is used to obtain a first estimate of the dispersion parameter. Furthermore, if the function \code{\link[stats]{glm}} does not converge the function \code{\link{dglars}} is used to compute the maximum likelihood estimates.
}

\value{\code{grcv} returns the estimate of the dispersion parameter.}

\references{
Pazira H., Augugliaro L. and Wit E.C. (2018) <doi:10.1007/s11222-017-9761-7>
\emph{Extended differential-geometric LARS for high-dimensional GLMs with general dispersion parameter},
\emph{Statistics and Computing}, Vol 28(4), 753-774.
}

\author{Luigi Augugliaro and Hassan Pazira\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}

\seealso{
\code{\link{phihat}}, \code{\link{AIC.dglars}} and \code{\link{BIC.dglars}}.
}

\examples{
############################
# y ~ Gamma
set.seed(321)
n <- 100
p <- 50
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[,1]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))

phi
grcv(fit, type = "AIC")
grcv(fit, type = "BIC")
}
\keyword{models}
\keyword{regression}
