\name{dglars}
\alias{dglars}
\title{dgLARS solution curve for GLM}
\description{
\code{dglars} function is used to estimate the solution curve implicitly defined by the dgLARS method for logistic and Poisson regression model.
}
\usage{
dglars (X, y, family = c("binomial", "poisson"), control = list())
}
\arguments{
  \item{X}{Design matrix of dimension \eqn{n\times p}.}
  \item{y}{Response vector.}
  \item{family}{A description of the error distribution used in the model (see below for more details).}
  \item{control}{A list of control parameters. See 'Details'.}
 }
\details{
\code{dglars} function implements the differential geometric generalization of the least angle regression method (Efron et al., 2004) proposed
in Augugliaro et al. (2013). Actual version of the package can be used to estimate the solution curve for a logistic regression model 
(\code{family = "binomial"}) and for a Poisson regression model (\code{family = "poisson"}).
\cr\cr
The dgLARS solution curve can be estimated using two different algorithms, i.e. the predictor-corrector method and the cyclic coordinate descent method (see below for 
more details about the control parameter \code{algorithm}). The first algorithm is based on two steps. In the first step, called predictor step, an approximation of the point 
that lies on the solution curve is computed. If the control parameter \code{dg_max} is equal to zero, in this step it is also computed an approximation of the optimal 
step size using a generalization of the method proposed in Efron et al. (2004). The optimal step size is defined as the reduction of the tuning parameter, denoted by
\eqn{d\gamma}, such that at \eqn{\gamma-d\gamma} there is a change in the active set. In the second step, called corrector step, a Newton-Raphson algorithm is used to
correct the approximation to the solution point  computed in the previous step. The main problem of this algorithm is that the number of arithmetic operations required to 
compute the approximation of the point that lies on the solution curve scales as the cube of the variables, this means that such algorithm is cumbersome in a high dimensional 
setting. To overcome this problem, the second algorithm compute the dgLARS solution curve using an adaptive version of the cyclic coordinate descent method proposed 
in Friedman et al. (2010).
\cr\cr
The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{algorithm}}{a string to specify the algorithm used to fit the dgLARS solution curve. If \code{algorithm = "pc"} (default)
	the predictor-corrector method is used, while the cyclic coordinate descent method is used if \code{algorithm = "ccd"};}
	\item{\code{method}}{a string to specify the method used to define the dgLARS solution curve. If \code{method = "dgLASSO"} (default)
	the algorithm computes the solution curve defined by the differential geometric generalization of the LASSO estimator; otherwise, if \code{method = "dgLAR"}, the
	differential geometric generalization of the least angle regression method is used;}
	\item{\code{nv}}{control parameter for the \code{pc} algorithm. An integer value belonging to the interval \eqn{[1;min(n,p)]} used to specify the maximum number 
	of variables included in the final model. Default is \code{nv = min(n-1,p)};}
	\item{\code{np}}{control parameter for the \code{pc/ccd} algorithm. A non negative integer used to define the maximum number of points of the solution curve. For the
	predictor-corrector algorithm \code{np} is set to \eqn{50 \cdot min(n-1,p)} (default) while for the cyclic coordinate descent method is set to 100 (default), i.e. the number 
	of values of the tuning parameter \eqn{\gamma};}
	\item{\code{g0}}{control parameter for the \code{pc/ccd} algorithm. Set the smallest value for the tuning parameter \eqn{\gamma}. Default is \code{g0 = ifelse(p<n, 1.0e-04, 0.05)};}
	\item{\code{dg_max}}{control parameter for the \code{pc} algorithm. A non negative value used to specify the maximum length of the step size. Setting \code{dg_max = 0} 
	(default) the predictor-corrector algorithm computes an approximation of the optimal step size (see Augugliaro et al. (accepted) for more details);}
	\item{\code{nNR}}{control criterion parameter for the \code{pc} algorithm. A non negative integer used to specify the maximum number of iterations of the Newton-Raphson algorithm 
	used in the corrector step. Default is \code{nNR = 50};}
	\item{\code{NReps}}{control parameter for the \code{pc} algorithm. A non negative value used to define the convergence of the Newton-Raphson algorithm. Default is 
	\code{NReps = 1.0e-06};}
	\item{\code{ncrct}}{control parameter for the \code{pc} algorithm. When one of the following conditions is satisfied
	\describe{
		\item{\code{i.}}{the Newton-Raphson algorithm does not converge}
		\item{\code{ii.}}{exists a non active variable such that, at the solution point, the absolute value of the corresponding Rao's score test statistics is greater than 
		\eqn{\gamma + }\code{eps}}	
	}
	then the step size (\eqn{d\gamma}) is reduced by \eqn{d\gamma = cf \cdot d\gamma} and the corrector step is repeated. \code{ncrct} is a non negative integer used to specify 
	the maximum number of trials of the corrector step. Default is \code{ncrct = 50};
	}
	\item{\code{cf}}{control parameter for the \code{pc} algorithm. The contractor factor is a real value belonging to the interval \eqn{[0,1]} used to reduce the step size as previously
	described. Default is \code{cf = 0.5};}
	\item{\code{nccd}}{control parameter for the \code{ccd} algorithm. A non negative integer used to specify the maximum number of steps of the cyclic coordinate descent algorithm.
	Default is \code{1.0e+05}.}
	\item{\code{eps}}{control parameter for the \code{pc/ccd} algorithm. The meaning of this parameter is related to the algorithm used to estimate the dgLARS solution curve, namely
	\describe{
		\item{\code{i.}}{when \code{algorithm = "pc"}, \code{eps} is used
		\describe{
			\item{\code{a.}}{to identify a variable that will be included in the active set, i.e. when the absolute value of the corresponding Rao's score test statistic belongs to 
			\eqn{[\gamma-\code{eps},\gamma+\code{eps}]};}
			\item{\code{b.}}{as previously described, to establish if the corrector step must be repeated;}
			\item{\code{c.}}{to define the convergence of the algorithm, i.e. the actual value of the tuning parameter belongs to the interval \eqn{[\code{g0-eps},\code{g0+eps}];}}
		}
		}
		\item{\code{ii.}}{when \code{algorithm = "ccd"}, \code{eps} is used to define the convergence of a single solution point, i.e. each inner 
		coordinate-descent loop continues until the maximum change in the Rao's score test statistic, after any coefficient update, is less than \code{eps}.}
	}
	Default is \code{eps = 1.0e-05.}
	}
	}
}
\value{
\code{dglars} returns an object with S3 class \code{"dglars"}, i.e. a list containing the following components:
	\item{call}{the call that produced this object;}
	\item{family}{a description of the error distribution used in the model;}
	\item{np}{the number of points of the dgLARS solution curve;}
	\item{beta}{the \eqn{(p+1)\times\code{np}} matrix corresponding to the dgLARS solution curve;}
	\item{ru}{the matrix of the Rao's score test statistics of the variables included in the final model;}
	\item{dev}{the \code{np} dimensional vector of the deviance corresponding to the values of the tuning parameter \eqn{\gamma};}
	\item{df}{the sequence of number of nonzero coefficients for each value of the tuning parameter \eqn{\gamma};}
	\item{g}{the sequence of \eqn{\gamma} values used to compute the solution curve;}
	\item{X}{the used design matrix;}
	\item{y}{the used response vector;}
	\item{action}{a \code{np} dimensional vector of characters used to show how is changed the active set for each value of the tuning parameter \eqn{\gamma};}
	\item{conv}{an integer value used to encode the warnings and the errors related to the algorithm used to compute the solution curve. The values returned are:
		        \describe{
	         			\item{\code{0}}{convergence of the algorithm has been achieved,}
	         			\item{\code{1}}{problems related with the predictor-corrector method: error in predictor step,}
	         			\item{\code{2}}{problems related with the predictor-corrector method: error in corrector step,}
	         			\item{\code{3}}{maximum number of iterations has been reached,}
	         			\item{\code{4}}{error in dynamic allocation memory;}
	         		}
		}
	\item{control}{the list of control parameters used to compute the dgLARS solution curve.}
}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2013)
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, \emph{Journal of the Royal Statistical Society. Series B.} Vol 75(3), 471-498.\cr\cr
Augugliaro L., Mineo A.M. and Wit E.C. (2012)
\emph{Differential geometric LARS via cyclic coordinate descent method}, in \emph{Proceeding of COMPSTAT 2012}, pp. 67-79. Limassol, Cyprus.
\cr\cr
Efron B., Hastie T., Johnstone I. and Tibshirani R. (2004)
\emph{Least Angle Regression}, \emph{The Annals of Statistics} Vol. 32(2), 407-499.
\cr\cr
Friedman J., Hastie T. and Tibshirani R. (2010) 
\emph{Regularization Paths for Generalized Linear Models via Coordinate Descent}, \emph{Journal of Statistical Software}, Vol. 33(1), 1-22.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{coef.dglars}}, \code{\link{plot.dglars}}, \code{\link{print.dglars}} and \code{\link{summary.dglars}} methods.
}
\examples{
#############################
# Logistic regression model #

set.seed(123)

# low dimensional setting
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
system.time(fit <- dglars(X, y, family = "binomial"))
system.time(
fit <- dglars(X, y, family = "binomial", control = list(algorithm = "ccd"))
)

# high dimensional setting
n <- 100
p <- 1000
X <- matrix(rnorm(n*p), n, p)
b <- 1:2
eta <- b[1] + X[,1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
system.time(fit <- dglars(X, y, family = "binomial"))
system.time(
fit <- dglars(X, y, family = "binomial", control = list(algorithm = "ccd"))
)
}
\keyword{models}
\keyword{regression}
