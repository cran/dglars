print.loglik_dglars <- function (x, digits = max(3, getOption("digits") - 3), ...){
    model <- x$object$control$method
    fml <- x$object$family
    tbl <- data.frame(g = x$g, df = x$df, loglik = x$loglik)
    cat("Sequence of log-likelihood values for the estimated", sQuote(model),"models.")
    cat("\n\nDetails:\n")
    cat("\t", sQuote(fml$family) ,"family with link function", sQuote(fml$link),"\n")
    if(!fml$family %in% c("binomial", "poisson"))
    cat("\t (dispersion parameter estimated by", sQuote(x$phi), ")\n\n")
        else
    cat("\t (dispersion parameter for", fml$family, "family taken to be 1)\n\n")
    print.data.frame(tbl, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
    invisible(tbl)
}
