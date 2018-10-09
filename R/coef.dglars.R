coef.dglars <- function(object, type = c("pearson", "deviance", "mle", "grcv"), g = NULL, ...){
    type <- match.arg(type)
    out <- predict(object, g = g, type = "coefficients")
    if(type != "pearson") out$phi <- phihat(object, type = type, g = g, ...)
    out$g <- if(is.null(g)) object$g
                else g
    out
}
