make_coef <- function(object){
	bhat <- object$b
	np <- object$np
	p <- object$p
	bhat.mat <- Matrix(bhat, nrow = p + 1, sparse = TRUE)
	bhat.mat <- bhat.mat[, 1:np, drop=FALSE]
	vname <- c("Int.",colnames(object$X))
	rownames(bhat.mat) <- vname
	bhat.mat
}