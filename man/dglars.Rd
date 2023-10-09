\name{dglars}
\alias{dglars}
\alias{dglars.fit}
\title{dgLARS Solution Curve for GLM}
\description{
\code{dglars} function is used to estimate the solution curve defined by dgLARS method.
}
\usage{
dglars(formula, family = gaussian, g, unpenalized, 
b_wght, data, subset, contrasts = NULL, control = list())

dglars.fit(X, y, family = gaussian, g, unpenalized, 
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
    
\item{g}{argument available only for \code{ccd} algorithm. When the model
is fitted by using the \code{ccd} algorithm, this argument can be used to
specify the values of the tuning parameter.}

\item{unpenalized}{a vector used to specify the unpenalized estimators;
\code{unpenalized} can be a vector of integers or characters specifying
the names of the predictors with unpenalized estimators (see example
below for more details).}
    	
\item{b_wght}{a \eqn{p+1}{p+1}-dimensional vector used to compute
the weights in the adaptive dgLARS method. \code{b_wght} is used
to specify the initial estimates of the parameter vector.}
    		
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
\code{dglars} function implements the differential geometric generalization 
of the least angle regression method (Efron et al., 2004) proposed in 
Augugliaro et al. (2013) and Pazira et al. (2017).

As in \dQuote{\code{\link{glm}}}, the user can specify family and link function using
the argument \code{family}. When the \code{binomial} family is used, the responce
can be a vector with entries 0/1 (failure/success) or, alternatively, a matrix where
the first column is the number of \dQuote{successes} and the second column is the number
of \dQuote{failures}. Starting with the version 2.0.0, the model can be specified combining
family and link functions as describted in the following table:

\tabular{ll}{
\bold{Family}   \tab \bold{Link} \cr
\code{gaussian} \tab \sQuote{\code{identity}}, \sQuote{\code{log}} and \sQuote{\code{inverse}} \cr
\code{binomial} \tab \sQuote{\code{logit}}, \sQuote{\code{probit}}, \sQuote{\code{cauchit}}, \sQuote{\code{log}} and \sQuote{\code{cloglog}} \cr
\code{poisson}  \tab \sQuote{\code{log}}, \sQuote{\code{identity}}, and \sQuote{\code{sqrt}} \cr
\code{Gamma}    \tab \sQuote{\code{inverse}}, \sQuote{\code{identity}} and \sQuote{\code{log}} \cr
\code{inverse.gaussian} \tab \sQuote{\code{1/mu^2}}, \sQuote{\code{inverse}}, \sQuote{\code{identity}}, and \sQuote{\code{log}}
}

The \code{R} code for binomial, Gamma and inverse gaussian families is due to
Hassan Pazira while the fortran version is due to Luigi Augugliaro.

\code{dglars.fit} is a workhorse function: it is more efficient when the design
matrix does not require manipulations. For this reason we suggest to use this function
when the dgLARS method is applied in a high-dimensional setting, i.e., when \code{p>n}.

When gaussian, gamma or inverse.gaussian is used to model the error distribution, \code{dglars}
returns the vector of the estimates of the dispersion parameter \eqn{\phi}{phi}; by
default, the generalized Pearson statistic is used as estimator but the user can use
the function \code{phihat} to specify other estimators (see \code{\link{phihat}} for
more details).

The dgLARS solution curve can be estimated using two different algorithms, i.e. the 
predictor-corrector method and the cyclic coordinate descent method (see below for 
more details about the argument \code{algorithm}). The first algorithm is
based on two steps. In the first step, called predictor step, an approximation of 
the point that lies on the solution curve is computed. If the control parameter 
\code{dg_max} is equal to zero, in this step it is also computed an approximation 
of the optimal step size using a generalization of the method proposed in Efron 
et al. (2004). The optimal step size is defined as the reduction of the tuning parameter, 
denoted by \eqn{d\gamma}{dg}, such that at \eqn{\gamma-d\gamma}{g - dg} there is a change in the
active set. In the second step, called corrector step, a Newton-Raphson algorithm is used to
correct the approximation previously computed. The main problem of this algorithm is that the
number of arithmetic operations required to compute the approximation scales as the cube
of the variables, this means that such algorithm is cumbersome in a high dimensional setting. 
To overcome this problem, the second algorithm compute the dgLARS solution curve using an 
adaptive version of the cyclic coordinate descent method proposed in Friedman et al. (2010).

The argument \code{control} is a list that can supply any of the following components:

\describe{
	
\item{\code{algorithm}:}{a string specifying the algorithm used to
compute the solution curve. The predictor-corrector algorithm is used
when \code{algorithm = ''pc''} (default), while the cyclic coordinate
descent method is used setting \code{algorithm = ''ccd''};}
		
\item{\code{method}:}{a string by means of to specify the kind of solution curve.
If \code{method = ''dgLASSO''} (default) the algorithm computes the solution
curve defined by the differential geometric generalization of the LASSO
estimator; otherwise (\code{method = ''dgLARS''}) the differential geometric
generalization of the least angle regression method is used;}
		
\item{\code{nv}:}{control parameter for the \code{pc} algorithm. An integer value
between 1 and \eqn{\min(n,p)}{min(n, p)} used to specify the maximum number of
variables in the final model. Default is \code{nv = min(n - 1, p)};}
		
\item{\code{np}:}{control parameter for the \code{pc/ccd} algorithm. A non negative 
integer used to define the maximum number of solution points. For the predictor-corrector 
algorithm \code{np} is set to \eqn{50\times\min(n - 1, p)}{50 * min{n-1, p}} (default); for 
the cyclic coordinate descent method, if \code{g} is not specified, this argument is
set equal to 100 (default);}
		
\item{\code{g0}:}{control parameter for the \code{pc/ccd} algorithm. This parameter is 
used to set the smallest value for the tuning parameter \eqn{\gamma}{g}. Default is 
\code{g0 = ifelse(p < n, 1.0e-04, 0.05)}; this argument is not required when \code{g} is
used with the cyclic coordinate descent algorithm;}
		
\item{\code{dg_max}:}{control parameter for the \code{pc} algorithm. A non negative 
value used to specify the largest value for the step size. Setting \code{dg_max = 0} 
(default) the predictor-corrector algorithm computes an approximation of the optimal 
step size (see Augugliaro et al. (2013) for more details);}
		
\item{\code{nNR}:}{control criterion parameter for the \code{pc} algorithm. A non 
negative integer used to specify the maximum number of iterations of the Newton-Raphson 
algorithm. Default is \code{nNR = 50};}
		
\item{\code{NReps}:}{control parameter for the \code{pc} algorithm. A non negative 
value used to define the convergence of the Newton-Raphson algorithm. Default is 
\code{NReps = 1.0e-06};}
	
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
	
\code{dglars} returns an object with S3 class \dQuote{\code{dglars}}, i.e., a list containing the 
following components:

\item{call}{the call that produced this object;}
	
\item{formula}{if the model is fitted by \code{dglars}, the used formula is returned;}
	
\item{family}{a description of the error distribution used in the model;}

\item{unpenalized}{the vector used to specify the unpenalized estimators;}

\item{np}{the number of points of the dgLARS solution curve;}
	
\item{beta}{the \eqn{(p + 1)\times\code{np}}{(p + 1) * np} matrix corresponding to the 
dgLARS solution curve;}
		
\item{phi}{the \code{np} dimensional vector of the Pearson estimates of the disperion 
parameter;}
		
\item{ru}{the matrix of the Rao's score test statistics of the variables included in 
the final model. This component is reported only if the predictor-corrector algorithm 
is used to fit the model;}
		
\item{dev}{the \code{np} dimensional vector of the deviance corresponding to the 
values of the tuning parameter \eqn{\gamma}{g};}
		
\item{nnonzero}{the sequence of number of nonzero coefficients for each value of the 
tuning parameter \eqn{\gamma}{g};}
		
\item{g}{the sequence of \eqn{\gamma}{g} values used to compute the solution curve;}
	
\item{X}{the used design matrix;}
	
\item{y}{the used response vector;}
	
\item{w}{the vector of weights used to compute the adaptive dglars method;}
	
\item{action}{a \code{np} dimensional vector of characters used to show how is 
changed the active set for each value of the tuning parameter \eqn{\gamma}{g};}
		
\item{conv}{an integer value used to encode the warnings and the errors related 
to the algorithm used to fit the model. The values returned are:
	\describe{
		\item{\code{0}}{convergence of the algorithm has been achieved;}
	         			
		\item{\code{1}}{problems related with the predictor-corrector 
						method: error in predictor step;}
	         				
		\item{\code{2}}{problems related with the predictor-corrector 
						method: error in corrector step;}
	         				
		\item{\code{3}}{maximum number of iterations has been reached;}
	         			
		\item{\code{4}}{error in dynamic allocation memory;}
	         			
		\item{\code{5}}{fitted expected value is out of range;}
	         			
		\item{\code{6}}{does not exist dgLARS estimator;}
	         			
		\item{\code{7}}{maximum number of solution points ('code{np}') reached.}
	}
}

\item{control}{the list of control parameters used to compute the dgLARS solution curve.}
}
\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2016) <doi:10.1093/biomet/asw023>
\emph{A differential-geometric approach to generalized linear models with grouped predictors}, 
\emph{Biometrika}, Vol 103(3), 563-577.
	
Augugliaro L., Mineo A.M. and Wit E.C. (2014) <doi:10.18637/jss.v059.i08>
\emph{dglars: An R Package to Estimate Sparse Generalized Linear Models}, \emph{Journal of Statistical Software}, Vol 59(8), 1-40. \url{https://www.jstatsoft.org/v59/i08/}.

Augugliaro L., Mineo A.M. and Wit E.C. (2013) <doi:10.1111/rssb.12000>
\emph{dgLARS: a differential geometric approach to sparse generalized linear models}, \emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.

Efron B., Hastie T., Johnstone I. and Tibshirani R. (2004) <doi:10.1214/009053604000000067>
\emph{Least Angle Regression}, \emph{The Annals of Statistics}, Vol. 32(2), 407-499.

Friedman J., Hastie T. and Tibshirani R. (2010) <doi:10.18637/jss.v033.i01>
\emph{Regularization Paths for Generalized Linear Models via Coordinate Descent}, \emph{Journal of Statistical Software}, Vol. 33(1), 1-22.

Pazira H., Augugliaro L. and Wit E.C. (2018) <doi:10.1007/s11222-017-9761-7>
\emph{Extended di erential geometric LARS for high-dimensional GLMs with general dispersion parameter},
\emph{Statistics and Computing}, Vol. 28(4), 753-774.
}
\author{Luigi Augugliaro and Hassan Pazira\cr
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{coef.dglars}}, \code{\link{phihat}}, \code{\link{plot.dglars}}, \code{\link{print.dglars}} and \code{\link{summary.dglars}} methods.
}
\examples{
set.seed(123)

#############################
# y ~ Binomial
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

############################
# y ~ Poisson 
n <- 100
p <- 100
X <- matrix(rnorm(n * p), n, p)
eta <- 2 + 2 * X[,1]
mu <- poisson()$linkinv(eta)
y <- rpois(n, mu)
fit <- dglars(y ~ X, family = poisson)
fit

############################
# y ~ Gamma
n <- 100
p <- 100
X <- matrix(abs(rnorm(n*p)),n,p)
eta <- 1 + 2 * X[,1]
mu <- drop(Gamma()$linkinv(eta))
shape <- 0.5
phi <- 1 / shape
y <- rgamma(n, scale = mu / shape, shape = shape)
fit <- dglars(y ~ X, Gamma("log"))
fit
}
\keyword{models}
\keyword{regression}
