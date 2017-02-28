coef.dglars <- function(object, type = c("pearson", "deviance", "mle"), ...){
    type <- match.arg(type)
    dots <- list(...)
    out <- predict(object, g = dots$g, type = "coefficients")
    if(type != "pearson") out$phi <- phihat(object, type = type, g = dots$g)
    out$g <- if(is.null(dots$g)) object$g
                else dots$g
    out
}
