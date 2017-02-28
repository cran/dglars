d2mu_de2_mk <- function(link){
	switch(link,
		log = function(eta, mi) mi * pmax(exp(eta), .Machine$double.eps),
		inverse = function(eta, mi) 2 / eta^3,
		sqrt = function(eta, mi) rep.int(2, length(eta)),
		cloglog = function(eta, mi) mi * (1 - exp(eta)) * exp(eta - exp(eta)),
		probit = function(eta, mi) - mi * eta * pmax(dnorm(eta), .Machine$double.eps),
		cauchy = function(eta, mi) -2 * mi * eta * pmax(dcauchy(eta) / (1 + eta^2), .Machine$double.eps)
	)
}
