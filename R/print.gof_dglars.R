print.gof_dglars <- function (x, digits = max(3, getOption("digits") - 3), ...){
    type <- x$type
    model <- x$object$control$method
    fml <- x$object$family
	phi <- x$phi
    g <- x$g
    npar <- x$npar
    val <- x$val
    rnk <- rank(val)
    rnk_min <- which.min(rnk)
    rnk <- as.character(rnk)
    rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
    rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
    tbl <- data.frame(g = g, npar = npar, val = val, Rank = rnk)
    names(tbl)[3] <- type
	cat("Sequence of", type, "values for the estimated", sQuote(model),"models.")
	cat("\n\nDetails:\n") 
	cat("\t", sQuote(fml$family) ,"family with link function", sQuote(fml$link),"\n")
    if(!fml$family %in% c("binomial", "poisson"))
        cat("\t (dispersion parameter estimated by", sQuote(phi), ")\n")
    else
        cat("\t (dispersion parameter for",fml$family, "family taken to be 1)\n")
    cat("\t", type, "values computed using k =", round(x$k, digits), "and complexity =", sQuote(x$complexity),"\n\n")
    print.data.frame(tbl, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	invisible(tbl)
}
