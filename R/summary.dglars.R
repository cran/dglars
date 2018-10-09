summary.dglars <- function(object, type = c("AIC", "BIC"), digits = max(3, getOption("digits") - 3), ...){
	if(class(object)[1] != "dglars")
        stop("summary method is not yet available for objects with class 'cvdglars'")
    if(is.character(type)) type <- match.arg(type)
    else stop("'type' argument is not a string specifying 'AIC or 'BIC'")
    tbl <- make_summary_table(object, type, ...)
    action <- object$action
    id <- which(action != "")
    n.tbl <- dim(tbl$table)[1L]
    n.space <- length(id)
    id.tbl <- vector(length = n.tbl + n.space, mode = "numeric")
    id.space <- id + seq(1, n.space)
    id.tbl[-id.space] <- seq(1:n.tbl)
    id.tbl[id.space] <- id
    tbl.format <- format(tbl$table[id.tbl, ], digits = digits)
    tbl.format[id.space - 1, 1] <- ""
    tbl.format[id.space, -1] <- ""
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    print.data.frame(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE)
    cat("\nDetails:")
    cat("\n\t", tbl$type, "values computed using k =", format(tbl$k, digits = digits), "and complexity =",sQuote(tbl$complexity))
    if(!object$family$family %in% c("binomial", "poisson"))
        cat("\n\t dispersion parameter estimated by", sQuote(tbl$phi))
    else
        cat("\n\t dispersion parameter for",object$family$family, "family taken to be 1")
    cat("\n\n===============================================================")
    cat("\n\nSummary of the Selected Model")
    if(!is.null(tbl$formula.gof)){
        cat("\n\n    Formula: ")
        print(tbl$formula.gof)
    } else cat("\n\n")
    cat("     Family:", sQuote(object$family$family))
    cat("\n       Link:", sQuote(object$family$link))
    if(!is.null(object$unpenalized))
    cat("\nUnpenalized:", paste(colnames(object$X)[object$unpenalized], collapse = ", "))
    cat("\n\nCoefficients:\n")
    coef_tbl <- as.matrix(tbl$b.gof)
    colnames(coef_tbl) <- "Estimate"
    printCoefmat(coef_tbl)
    if(!object$family$family %in% c("binomial", "poisson"))
        cat("\nDispersion parameter:", format(tbl$phi.gof, digits = digits), "(estimated by", sQuote(tbl$phi),"method)\n")
    cat("---")
    cat("\n\n                 g:", format(c(tbl$g.gof), digits = digits))
    cat("\n",
        apply(cbind(format(c("Null deviance:", "Residual deviance:", paste(tbl$type, ":", sep = "")), justify = "right"),
            format(c(tbl$nulldev, tbl$resdev.gof, min(tbl$table[5])), digits = digits),
            rep("\n", 3L)), 1L, paste, collapse = " "))
    cat("\n Algorithm", sQuote(object$control$algorithm),"( method =",sQuote(object$control$method),")\n\n")
    invisible(tbl)
}
