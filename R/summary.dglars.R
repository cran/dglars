summary.dglars <- function (object,complexity=c("df","gdf"),digits = max(3, getOption("digits") - 3),...){
	complexity <- match.arg(complexity)
	if(object$family=="poisson" & complexity == "gdf"){
		complexity <- "df"
		warning("'complexity' was set equal to 'df' because for Poisson regression the\n model complexity can be approximated by the number of nonzero coefficients")
	}
	tbl <- make_summary_table(object,complexity)
	cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.data.frame(format(tbl$table, digits = digits), print.gap = 2, quote = FALSE,row.names=FALSE,...)
	cat("\n========================================\n")
	cat("\nBest model indentified by AIC criterion: ", paste(names(tbl$b.aic),collapse="+"),"\n")
	cat("\nCoefficients:\n\n")
	print.default(format(tbl$b.aic, digits = digits), print.gap = 2, quote = FALSE,...)
	cat("\n\nAIC: ",format(min(tbl$table$AIC), digits = digits))
	cat("\n\n========================================\n")
	cat("\nBest model indentified by BIC criterion: ", paste(names(tbl$b.bic),collapse="+"),"\n")
	cat("\nCoefficients:\n\n")
	print.default(format(tbl$b.bic, digits = digits), print.gap = 2, quote = FALSE,...)
	cat("\n\nBIC: ",format(min(tbl$table$BIC), digits = digits))	
	cat("\n\n===\n\nAlgorithm", object$control$algorithm,"( method =",object$control$method,") with exit =",object$conv,"\n\n")
	invisible(tbl)
}
