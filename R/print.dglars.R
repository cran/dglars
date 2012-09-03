print.dglars <- function (x,digits = max(3, getOption("digits") - 3), ...){
	b <- x$beta
	action <- x$action
	g <- x$g
	dev <- x$dev
	dev.ratio <- 1 - dev/dev[1]
	df <- x$df
	tbl <- data.frame(action,g,dev,dev.ratio,df)
	names(tbl) <- c("Sequence","g","Dev","%Dev","df")
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.data.frame(format(tbl, digits = digits),print.gap = 2, quote = FALSE,row.names=FALSE, ...)
	cat("\nAlgorithm", x$control$algorithm,"( method =",x$control$method,") with exit =",x$conv,"\n\n")
	invisible(tbl)
}
