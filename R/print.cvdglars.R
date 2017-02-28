print.cvdglars <- function (x, digits = max(3, getOption("digits") - 3), ...){
	b <- x$beta[abs(as.vector(x$beta)) > 0, , drop = FALSE]
	colnames(b) <- "Estimate"
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
	cat("\n\nCoefficients:\n")
	printCoefmat(b)
	if(!x$family$family %in% c("binomial", "poisson"))
        cat("\nDispersion parameter:", format(x$phi, digits = digits))
    else
        cat("\ndispersion parameter for", x$family$family, "family taken to be 1")
    cat("\n\nDetails:")
    cat("\n   number of non zero estimates:", ifelse(!x$family$family %in% c("binomial", "poisson"), dim(b)[1] + 1, dim(b)[1]))
    cat("\n      cross-validation deviance:", format(min(x$dev_m), digits = digits))
    cat("\n                              g:", format(min(x$g), digits = digits))
    cat("\n                        n. fold:", x$control$nfold)
	cat("\n\nAlgorithm", sQuote(x$control$algorithm),"( method =",sQuote(x$control$method),")\n\n")
	invisible(b)
}
