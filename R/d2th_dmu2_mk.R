d2th_dmu2_mk <- function(fml){
	switch(fml,
		binomial = function(mu, mi) 1 / (mi - mu)^2 - 1 / mu^2,
		poisson = function(mu, mi) - 1 / mu^2,
		Gamma = function(mu, mi) - 2 / mu^3,
		inverse.gaussian = function(mu, mi) - 3 / mu^4
	)
}
