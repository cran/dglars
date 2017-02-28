BIC.dglars <- function(object, ...){
    arg <- names(list(...))
    if("k" %in% arg) out_aic <- AIC(object, ...)
    else{
        logn <- log(dim(object$X)[1])
        out_aic <- AIC(object, k = logn, ...)
        out_aic$type <- "BIC"
    }
    out_aic
}
