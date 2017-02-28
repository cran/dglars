make_cvdglars <- function(object, control){
	beta <- Matrix(object$b)
	var_cv <- colnames(object$X)
	rownames(beta) <- c("Int.", var_cv)
	var_cv <- var_cv[abs(beta[-1, 1]) > 0]
	phi <- object$phi
	out <- list(call = NULL, formula_cv = NULL, family = NULL, var_cv = var_cv, beta = beta, phi = phi,
				dev_m = object$dev_m, dev_v = object$dev_v, g = object$g_hat, g0 = control$g0,
				g_max = object$g[1], X = Matrix(object$X), y = object$y, w = object$w, conv = object$conv, 
				control = control)
	class(out) <- c("cvdglars", "dglars")
	out
}
