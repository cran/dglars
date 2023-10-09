phihat <- function(object, type = c("pearson", "deviance", "mle", "grcv"), g = NULL, ...){
    if(!inherits(object, what = "dglars"))
        stop("this function works only with objects with class 'dglars'")
    family_used <- object$family$family
    type <- match.arg(type)
    if(family_used %in% c("binomial", "poisson"))
        phih <- rep.int(1, ifelse(is.null(g), length(object$g), length(g)))
    else{
        if(type != "grcv") {
            n <- dim(object$X)[1]
            npar <- predict(object, g = g, type = "nnonzero")
            df <- n - npar
            dev <- predict(object, g = g, type = "deviance")
            phih <- switch(type,
                            pearson = predict(object, g = g, type = "coefficients")$phi,
                            deviance = dev / df,
                            mle = if(family_used != "Gamma") dev / n
                                else 2 * dev / (n * (1 + sqrt(1 + 2/3 * dev / n))))
        } else {
            np <- ifelse(is.null(g), length(object$g), length(g))
            phih <- grcv(object, ...)
            phih <- rep(phih, length.out = np)
        }
    }
    phih
}
