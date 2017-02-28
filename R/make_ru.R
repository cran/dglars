make_ru <- function(object){
	ru <- object$ru
	np <- object$np
	p <- object$p
	nav <- object$nav
	A <- object$A[1:nav]
	ru <- Matrix(ru[1:(p*np)], nrow = p, ncol = np, sparse = FALSE)
	rownames(ru) <- colnames(object$X)
	ru[sort(A),]
}
