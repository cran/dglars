\name{dglars-package}
\alias{dglars-package}
\docType{package}

\title{Differential Geometric Least Angle Regression}

\description{
Differential geometric least angle regression method for fitting sparse generalized linear models. In this version of the package, the user can fit models specifying Gaussian, Poisson, Binomial, Gamma and Inverse Gaussian family. Furthermore, several link functions can be used to model the relationship between the conditional expected value of the response variable and the linear predictor. The solution curve can be computed using an efficient predictor-corrector or a cyclic coordinate descent algorithm, as described in the paper linked to via the URL below.
}

\details{
\tabular{ll}{
Package: \tab dglars\cr
Type: \tab Package\cr
Version: \tab 2.1.3\cr
Date: \tab 2019-10-09\cr
License: \tab GPL (>=2)
}
}

\author{
Luigi Augugliaro\cr
Maintainer: Luigi Augugliaro  <luigi.augugliaro@unipa.it>
}

\references{
Augugliaro L., Mineo A.M. and Wit E.C. (2016) <DOI:10.1093/biomet/asw023>
\emph{A differential-geometric approach to generalized linear models with grouped predictors}, Vol 103(3), 563-577.

Augugliaro L., Mineo A.M. and Wit E.C. (2014) <DOI:10.18637/jss.v059.i08>
\emph{dglars: An R Package to Estimate Sparse Generalized Linear Models}, \emph{Journal of Statistical Software}, Vol 59(8), 1-40. \url{http://www.jstatsoft.org/v59/i08/}.

Augugliaro L., Mineo A.M. and Wit E.C. (2013) <DOI:10.1111/rssb.12000>
\emph{dgLARS: a differential geometric approach to sparse generalized linear models},
\emph{Journal of the Royal Statistical Society. Series B.}, Vol 75(3), 471-498.
 
Augugliaro L., Mineo A.M. and Wit E.C. (2012)
\emph{Differential geometric LARS via cyclic coordinate descent method}, in \emph{Proceeding of COMPSTAT 2012}, pp. 67-79. Limassol, Cyprus.

Efron B., Hastie T., Johnstone I. and Tibshirani R. (2004) <DOI:10.1214/009053604000000067>
\emph{Least Angle Regression}, \emph{The Annals of Statistics}, Vol. 32(2), 407-499.

Pazira H., Augugliaro L. and Wit E.C. (2018) <DOI: 10.1007/s11222-017-9761-7>
\emph{Extended differential-geometric LARS for high-dimensional GLMs with general dispersion parameter},
\emph{Statistics and Computing}, Vol 28(4), 753-774.
}

\keyword{models}
\keyword{regression}
\keyword{package}
