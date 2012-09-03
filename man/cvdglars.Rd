\name{cvdglars}
\alias{cvdglars}
\title{Cross-validation deviance for dgLARS}
\description{
Uses the \eqn{k}-fold cross-validation deviance to estimate the solution point of the dgLARS solution curve.
}
\usage{
cvdglars(X, y, family = c("binomial", "poisson"), control = list())
}
\arguments{
  \item{X}{Design matrix of dimension \eqn{n\times p}.}
  \item{y}{Response vector.}
  \item{family}{A description of the error distribution used in the model.}
  \item{control}{A list of control parameters. See 'Details'.}
 }
\details{
\code{cvdglars} function runs \code{dglars} \code{nfold}+1 times. The deviance is stored, and the average and its standard deviation over the folds are computed.
\cr\cr
The \code{control} argument is a list that can supply any of the following components:
\describe{
	\item{\code{algorithm}}{a string to specify the algorithm used to fit the dgLARS solution curve. If \code{algorithm = "pc"} (default)
	the predictor-corrector method is used while the cyclic coordinate descent method is used if \code{algorithm = "ccd"};}
	\item{\code{method}}{a string to specify the method used to define the dgLARS solution curve. If \code{method = "dgLASSO"} (default)
	the algorithm computes the solution curve defined by the differential geometric generalization of the LASSO estimator; otherwise, if \code{method = "dgLAR"}, the
	differential geometric generalization of the least angle regression method is computed;}
	\item{\code{nfold}}{a non negative integer used to specify the number of folds. Although \code{nfolds} can be as large as the sample size (leave-one-out CV), it 
	is not recommended for large datasets. Default is \code{nfold = 10};}
	\item{\code{foldid}}{a \eqn{n}-dimensional vector of integers, between 1 and \eqn{n}, used to define the folds for the cross-validation. By default \code{foldid} is 
	randomly generated;}
	\item{\code{ng}}{number of values of the tuning parameter used to compute the cross-validation deviance. Default is \code{ng = 100};}
	\item{\code{nv}}{control parameter for the \code{pc} algorithm. An integer value belonging to the interval \eqn{[1;min(n,p)]} (default is \code{nv = min(n-1,p)}) used to 
	specify the maximum number of variables included in the final model;}
	\item{\code{np}}{control parameter for the \code{pc/ccd} algorithm. A non negative integer used to define the maximum number of points of the solution curve. For the
	predictor-corrector algorithm \code{np} is set to \eqn{50 \cdot min(n-1,p)} (default), while for the cyclic coordinate descent method is set to 100 (default), i.e. the number 
	of values of the tuning parameter \eqn{\gamma};}
	\item{\code{g0}}{control parameter for the \code{pc/ccd} algorithm. Set the smallest value for the tuning parameter \eqn{\gamma}. Default is \code{g0 = ifelse(p<n, 1.0e-04, 0.05)};}
	\item{\code{dg_max}}{control parameter for the \code{pc} algorithm. A non negative value used to specify the maximum length of the step size. Setting \code{dg_max = 0} 
	(default) the predictor-corrector algorithm uses the optimal step size (see Augugliaro et al. (accepted) for more details) to approximate the value of the tuning parameter corresponding to the 
	inclusion/exclusion of a variable from the model;}
	\item{\code{nNR}}{control parameter for the \code{pc} algorithm. A non negative integer used to specify the maximum number of iterations of the Newton-Raphson algorithm 
	used in the corrector step. Default is \code{nNR = 50};}
	\item{\code{NReps}}{control parameter for the \code{pc} algorithm. A non negative value used to define the convergence criterion of the Newton-Raphson algorithm. Default is 
	\code{NReps = 1.0e-06};}
	\item{\code{ncrct}}{control parameter for the \code{pc} algorithm. When one of the following conditions is satisfied
	\describe{
		\item{\code{i.}}{the Newton-Raphson algorithm does not converge}
		\item{\code{ii.}}{exist a non active variable such that, at the solution point, the absolute value of the corresponding Rao's score test statistics is greater than 
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
\code{cvdglars} returns an object with S3 class \code{"cvdglars"}, i.e. a list containing the following components:
	\item{call}{the call that produced this object;}
	\item{family}{a description of the error distribution used in the model;}
	\item{beta}{the vector of the coefficients estimated by cross-validation;}
	\item{dev_m}{a vector of length \code{ng} used to store the mean cross-validation deviance;}
	\item{dev_v}{a vector of length \code{ng} used to store the variance of the mean cross-validation deviance;}
	\item{g_hat}{the value of the tuning parameter corresponding to the minimum of the cross-validation deviance;}
	\item{X}{the used design matrix;}
	\item{y}{the used response vector;}
	\item{conv}{an integer value used to encode the warnings and the errors related to the algorithm used to dgLARS solution curve. The values returned are:
		        \describe{
         				\item{\code{0}}{convergence of the algorithm has been achieved,}
         				\item{\code{1}}{problems related with the predictor-corrector method: error in predictor step,}
         				\item{\code{2}}{problems related with the predictor-corrector method: error in corrector step,}
         				\item{\code{3}}{maximum number of iterations has been reached,}
         				\item{\code{4}}{error in dynamic allocation memory;}
         			}
		}
	\item{control}{the list of control parameters used to compute the cross-validation deviance.}
}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2013)
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, \emph{Journal of the Royal Statistical Society. Series B.} Vol 75(3), 471-498.\cr\cr
Augugliaro L., Mineo A.M. and Wit E.C. (2012)
\emph{Differential geometric LARS via cyclic coordinate descent method}, in \emph{Proceeding of COMPSTAT 2012}, pp. 67-79. Limassol, Cyprus.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{coef.cvdglars}}, \code{\link{print.cvdglars}}, \code{\link{plot.cvdglars}} methods
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
fit_cv <- cvdglars(X, y, family = "binomial")
fit <- dglars(X, y, family = "binomial", control = list(g0=fit_cv$g_hat))
fit_cv
fit$beta[,fit$np]
}
\keyword{models}
\keyword{regression}
