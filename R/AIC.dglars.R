AIC.dglars <- function(object, phi = c("pearson", "deviance", "mle", "grcv"), k = 2, complexity = c("df", "gdf"), g = NULL, ...){
    type <- ifelse(k == 2, "AIC", "GoF")
	phi <- match.arg(phi)
	complexity <- match.arg(complexity)
    out_loglik <- logLik(object, phi = phi, g = g, ...)
    loglik <- out_loglik$loglik
    comp <- if(complexity == "df") out_loglik$df
    		else{
    			out_gdf <- gdf(object)
    			if(!object$family$family %in% c("binomial", "poisson")) out_gdf <- out_gdf + 1
                out_gdf
    		}
    out <- list(val = - 2 * loglik + k * comp, g = out_loglik$g, loglik = loglik, k = k,
                comp = comp, npar = out_loglik$df, phi = phi, phih = out_loglik$phih,
                complexity = complexity, object = object, type = type)
    class(out) <- "gof_dglars"
    out
}
