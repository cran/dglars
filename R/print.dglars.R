print.dglars <- function (x, digits = max(3, getOption("digits") - 3), ...){
	action <- x$action
	g <- x$g
	dev <- x$dev
	dev.ratio <- 1 - dev/dev[1]
	nnonzero <- x$nnonzero
	tbl <- data.frame(action,g,dev,dev.ratio,nnonzero)
	names(tbl) <- c("Sequence","g","Dev","%Dev","n. non zero")
	if(x$np>1){
		id <- which(action != "")
		n.tbl <- dim(tbl)[1L]
		n.space <- length(id)
		id.tbl <- vector(length = n.tbl + n.space, mode = "numeric")
		id.space <- id + seq(1,n.space)
		id.tbl[-id.space] <- seq(1:n.tbl)
		id.tbl[id.space] <- id
		tbl.format <- format(tbl[id.tbl,], digits = digits)
		tbl.format[id.space-1,1] <- ""
		tbl.format[id.space,-1] <- ""
	} else tbl.format <- tbl
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.data.frame(tbl.format,print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	cat("\nAlgorithm", sQuote(x$control$algorithm),"( method =",sQuote(x$control$method),") with exit =",x$conv,"\n\n")
	invisible(tbl)
}
