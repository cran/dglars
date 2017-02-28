logLik.dglars <- function(object, phi = c("pearson", "deviance", "mle"), ...){
    phi <- match.arg(phi)
    dots <- list(...)
    g <- dots$g
    family_used <- object$family$family
    y <- object$y
    if(family_used == "binomial"){
        if(is.matrix(y)){
            mi <- y[, 1] + y[, 2]
            y <- y[, 1]
        } else mi <- rep(1, length(y))
        if(is.factor(y))
            y <- as.numeric(y) - 1
        prob <- predict(object, g = g, type = "probability")
    }
    else
        muh <- predict(object, g = g, type = "mu")
    npar <- predict(object, g = g, type = "nnonzero")
    if(!family_used %in% c("binomial", "poisson")){
    	n <- dim(object$X)[1]
        df <- n - npar
        dev <- predict(object, g = g, type = "deviance")
        phih <- switch(phi,
                    pearson = predict(object, g = g, type = "coefficients")$phi,
                    deviance = dev / df,
                    mle = if(family_used != "Gamma") dev / n
                            else 2 * dev / (n * (1 + sqrt(1 + 2/3 * dev / n)))
                    )
    }
    if(family_used == "gaussian"){
    	val <- -n * 0.5 * (1 + log(2 * pi * phih))
    	if(phi == "pearson") val <- val + npar / 2
    }
    if(family_used == "binomial")
        val <- apply(prob, 2, function(prb) sum(dbinom(y, mi, prb, log = TRUE)))
    if(family_used == "poisson")
        val <- apply(muh, 2, function(mu) sum(dpois(y, mu, log = TRUE)))
    if(family_used == "Gamma"){
        par <- rbind(phih, muh)
        val <- apply(par, 2, function(coef) sum(dgamma(y, 1 / coef[1], scale = coef[-1] * coef[1], log = TRUE)))
    }
    if(family_used == "inverse.gaussian"){
        val <- -0.5 * (n * log(2 * pi * phih) + 3 * sum(log(y)))
        val <- switch(phi,
                pearson = val - dev / (2 * phih),
                deviance = val - df/ 2,
                mle = val - n / 2
               )
    }
    if(is.null(g)) g <- object$g
    if(!family_used %in% c("binomial", "poisson")) npar <- npar + 1
    out <- list(loglik = val, df = npar, object = object, g = g, phi = phi)
    class(out) <- "loglik_dglars"
    out
}
