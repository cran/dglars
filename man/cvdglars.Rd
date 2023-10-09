\name{cvdglars}
\alias{cvdglars}
\alias{cvdglars.fit}

\title{Cross-Validation Method for dgLARS}

\description{
Uses the \eqn{k}-fold cross-validation deviance to estimate the solution point of the dgLARS solution curve.
}

\usage{
cvdglars(formula, family = gaussian, g, unpenalized, 
b_wght, data, subset, contrasts = NULL, control = list())

cvdglars.fit(X, y, family = gaussian, g, unpenalized,
b_wght, control = list())
}

\arguments{
\item{formula}{an object of class \dQuote{\code{\link{formula}}}:
a symbolic description of the model to be fitted. When the
\code{binomial} family is used, the responce can be a vector
with entries 0/1 (failure/success) or, alternatively, a
matrix where the first column is the number of \dQuote{successes}
and the second column is the number of \dQuote{failures}.}

\item{family}{a description of the error distribution and link
function used to specify the model. This can be a character string
naming a family function or the result of a call to a family function
(see \code{\link{family}} for details). By default the gaussian family
with identity link function is used.}

\item{g}{argument available only for \code{ccd} algorithm. When the \code{ccd}
algorithm is used to fit the dgLARS model, this argument can be used to specify
the values of the tuning parameter.}
 
\item{unpenalized}{a vector used to specify the unpenalized estimators;
\code{unpenalized} can be a vector of integers or characters specifying
the names of the predictors with unpenalized estimators.}

\item{b_wght}{a vector, with length equal to the number of columns of
the matrix \code{X}, used to compute the weights used in the
adaptive dgLARS method. \code{b_wght} is used to specify the
initial estimates of the parameter vector.}
  
\item{data}{an optional data frame, list or environment (or object
coercible by \sQuote{as.data.frame} to a data frame) containing the
variables in the model. If not found in \sQuote{data}, the
variables are taken from \sQuote{environment(formula)}.}
    	
\item{subset}{an optional vector specifying a subset of observations
to be used in the fitting process.}
    	
\item{contrasts}{an optional list. See the \sQuote{contrasts.arg} of
\sQuote{model.matrix.default}.}
    	
\item{control}{a list of control parameters. See \sQuote{Details}.}
    
\item{X}{design matrix of dimension \eqn{n\times p}.}
	
\item{y}{response vector. When the \code{binomial} family is used,
this argument can be a vector with entries 0 (failure) or 1
(success). Alternatively, the response can be a matrix where
the first column is the number of \dQuote{successes} and the second
column is the number of \dQuote{failures}.}
}

\details{
\code{cvdglars} function runs \code{dglars} \code{nfold}+1 times. 
The deviance is stored, and the average and its standard deviation 
over the folds are computed.

\code{cvdglars.fit} is the workhorse function: it is more efficient 
when the design matrix have already been calculated. For this reason 
we suggest to use this function when the dgLARS method is applied in 
a high-dimensional setting, i.e. when \code{p>n}.

The \code{control} argument is a list that can supply any of the following components:
\describe{
\item{\code{algorithm}:}{a string specifying the algorithm used to
compute the solution curve. The predictor-corrector algorithm is used
when \code{algorithm = ''pc''} (default), while the cyclic coordinate d
escent method is used setting \code{algorithm = ''ccd''};}
	
\item{\code{method}:}{a string by means of to specify the kind of solution curve.
If \code{method = ''dgLASSO''} (default) the algorithm computes the solution
curve defined by the differential geometric generalization of the LASSO
estimator; otherwise, if \code{method = ''dgLARS''}, the differential geometric
generalization of the least angle regression method is used;}
	
\item{\code{nfold}:}{a non negative integer used to specify the number of folds. 
Although \code{nfolds} can be as large as the sample size (leave-one-out CV), it 
is not recommended for large datasets. Default is \code{nfold = 10};}

\item{\code{foldid}}{a \eqn{n}-dimensional vector of integers, between 1 and \eqn{n}, 
used to define the folds for the cross-validation. By default \code{foldid} is 
randomly generated;}

\item{\code{ng}:}{number of values of the tuning parameter used to compute the 
cross-validation deviance. Default is \code{ng = 100};}

\item{\code{nv}:}{control parameter for the \code{pc} algorithm. An integer value 
belonging to the interval \eqn{[1;min(n,p)]} (default is \code{nv = min(n-1,p)}) 
used to specify the maximum number of variables included in the final model;}

\item{\code{np}:}{control parameter for the \code{pc/ccd} algorithm. A non negative 
integer used to define the maximum number of points of the solution curve. For the
predictor-corrector algorithm \code{np} is set to \eqn{50 \cdot min(n-1,p)} (default), 
while for the cyclic coordinate descent method is set to 100 (default), i.e. the number 
of values of the tuning parameter \eqn{\gamma}{g};}

\item{\code{g0}:}{control parameter for the \code{pc/ccd} algorithm. Set the smallest 
value for the tuning parameter \eqn{\gamma}{g}. Default is \code{g0 = ifelse(p<n, 1.0e-06, 0.05)};}

\item{\code{dg_max}:}{control parameter for the \code{pc} algorithm. A non negative value 
used to specify the maximum length of the step size. Setting \code{dg_max = 0} (default) 
the predictor-corrector algorithm uses the optimal step size (see Augugliaro et al. (2013) 
for more details) to approximate the value of the tuning parameter corresponding to the 
inclusion/exclusion of a variable from the model;}

\item{\code{nNR}:}{control parameter for the \code{pc} algorithm. A non negative integer 
used to specify the maximum number of iterations of the Newton-Raphson algorithm 
used in the corrector step. Default is \code{nNR = 200};}

\item{\code{NReps}:}{control parameter for the \code{pc} algorithm. A non negative 
value used to define the convergence criterion of the Newton-Raphson algorithm. 
Default is \code{NReps = 1.0e-06};}

\item{\code{ncrct}:}{control parameter for the \code{pc} algorithm. When the Newton-Raphson 
algorithm does not converge, the step size (\eqn{d\gamma}{dg}) is reduced by 
\eqn{d\gamma = cf \cdot d\gamma}{dg = cf * dg} and the corrector step is repeated. \code{ncrct} 
is a non negative integer used to specify the maximum number of trials for the corrector step. 
Default is \code{ncrct = 50};}


\item{\code{cf}:}{control parameter for the \code{pc} algorithm. The contractor factor 
is a real value belonging to the interval \eqn{[0,1]} used to reduce the step size 
as previously described. Default is \code{cf = 0.5};}
	
\item{\code{nccd}:}{control parameter for the \code{ccd} algorithm. A non negative integer 
used to specify the maximum number for steps of the cyclic coordinate descent algorithm.
Default is \code{1.0e+05}.}
	
\item{\code{eps}}{control parameter for the \code{pc/ccd} algorithm. The meaning of 
this parameter is related to the algorithm used to estimate the solution curve:
\describe{			
	\item{\code{i.}}{if \code{algorithm = ''pc''} then \code{eps} is used
	\describe{
		\item{\code{a.}}{to identify a variable that will be included in the active 
						set (absolute value of the corresponding Rao's score test 
						statistic belongs to 
						\eqn{[\gamma - \code{eps}, \gamma + \code{eps}]}{[g - eps, g + eps]});}
					
		\item{\code{b.}}{to establish if the corrector step must be repeated;}
				
		\item{\code{c.}}{to define the convergence of the algorithm, i.e., the 
						actual value of the tuning parameter belongs to the interval 
						\eqn{[\code{g0 - eps},\code{g0 + eps}]}{g0 - eps, g0 + eps};}
	}}
		
	\item{\code{ii.}}{if \code{algorithm = ''ccd''} then \code{eps} is used to define the 
	convergence for a single solution point, i.e., each inner coordinate-descent loop 
	continues until the maximum change in the Rao's score test statistic, after any 
	coefficient update, is less than \code{eps}.}}
		
Default is \code{eps = 1.0e-05.}}}
}

\value{
\code{cvdglars} returns an object with S3 class \dQuote{\code{cvdglars}}, i.e. a list 
containing the following components:

\item{call}{the call that produced this object;}

\item{formula_cv}{if the model is fitted by \code{cvdglars}, the used formula is returned;}

\item{family}{a description of the error distribution used in the model;}
	
\item{var_cv}{a character vector with the name of variables selected by cross-validation;}

\item{beta}{the vector of the coefficients estimated by cross-validation;}

\item{phi}{the cross-validation estimate of the disperion parameter;}

\item{dev_m}{a vector of length \code{ng} used to store the mean cross-validation 
deviance;}

\item{dev_v}{a vector of length \code{ng} used to store the variance of the mean 
cross-validation deviance;}

\item{g}{the value of the tuning parameter corresponding to the minimum of 
the cross-validation deviance;}

\item{g0}{the smallest value for the tuning parameter;}

\item{g_max}{the value of the tuning parameter corresponding to the starting 
point of the dgLARS solution curve;}

\item{X}{the used design matrix;}

\item{y}{the used response vector;}

\item{w}{the vector of weights used to compute the adaptive dglars method;}

\item{conv}{an integer value used to encode the warnings and the errors 
related to the algorithm used to fit the dgLARS solution curve. The values returned are:
\describe{
	\item{\code{0}}{convergence of the algorithm has been achieved,}
	\item{\code{1}}{problems related with the predictor-corrector method: error in predictor step,}
	\item{\code{2}}{problems related with the predictor-corrector method: error in corrector step,}
	\item{\code{3}}{maximum number of iterations has been reached,}
	\item{\code{4}}{error in dynamic allocation memory;}
}}
		
\item{control}{the list of control parameters used to compute the cross-validation deviance.}
}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2014) <doi:10.18637/jss.v059.i08>
\emph{dglars: An R Package to Estimate Sparse Generalized Linear Models}, \emph{Journal of Statistical Software}, Vol 59(8), 1-40. \url{https://www.jstatsoft.org/v59/i08/}.

Augugliaro L., Mineo A.M. and Wit E.C. (2013) <doi:10.1111/rssb.12000>
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, \emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{coef.cvdglars}}, \code{\link{print.cvdglars}}, \code{\link{plot.cvdglars}} methods
}
\examples{
###########################
# Logistic regression model
# y ~ Binomial
set.seed(123)
n <- 100
p <- 100
X <- matrix(rnorm(n * p), n, p)
b <- 1:2
eta <- b[1] + X[, 1] * b[2]
mu <- binomial()$linkinv(eta)
y <- rbinom(n, 1, mu)
fit_cv <- cvdglars.fit(X, y, family = binomial)
fit_cv
}
\keyword{models}
\keyword{regression}
