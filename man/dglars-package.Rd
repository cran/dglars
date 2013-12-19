\name{dglars-package}
\alias{dglars-package}
\docType{package}
\title{
Differential geometric LARS (dgLARS) method for some generalized linear models
}
\description{
  This package implements the differential geometric extension of the least angle 
  regression method (Efron et al., 2004) proposed in Augugliaro et al. (2013).
  Using this version of the package it is possible to compute the dgLARS solution
  curve for the logistic and Poisson regression models. The solution curve can
  be computed using two different algorithms, i.e., a predictor-corrector algorithm and a 
  cyclic coordinate descent algorithm.
}
\details{
\tabular{ll}{
Package: \tab dglars\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2013-12-19\cr
License: \tab GPL (>=2)
}
}
\author{
Luigi Augugliaro\cr
Maintainer: Luigi Augugliaro  <luigi.augugliaro@unipa.it>
}
\references{Augugliaro L., Mineo A.M. and Wit E.C. (2013)
\emph{dgLARS: a differential geometric approach to sparse generalized linear 
  models}, \emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.
 
Augugliaro L., Mineo A.M. and Wit E.C. (2012)
\emph{Differential geometric LARS via cyclic coordinate descent method}, in \emph{Proceeding of COMPSTAT 2012}, pp. 67-79. Limassol, Cyprus.

Efron B., Hastie T., Johnstone I. and Tibshirani R. (2004)
\emph{Least Angle Regression}, \emph{The Annals of Statistics}, Vol. 32(2), 407-499.
}
\keyword{models}
\keyword{regression}
\keyword{package}
