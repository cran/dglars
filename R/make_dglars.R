make_dglars <- function(object,control){
	np <- object$np
	nav <-object$nav
	beta <- make_coef(object)
	ru <- make_ru(object)
	unpenalized <- if(object$nup == 0) NULL
					else object$A[1:object$nup] 
	if(np>1) action <- make_action(beta[-1,], np)
	else action <- ""
	if(control$algorithm=="pc")
		out <- list(call = NULL, formula = NULL, family = NULL, unpenalized = unpenalized, np = np, beta = beta, 
					phi = object$phi[1:np], ru = ru, dev = object$dev[1:np], nnonzero = object$nnonzero[1:np], 
					g = object$g_seq[1:np], X = Matrix(object$X), y = object$y, w = object$w, action = action, 
					conv = object$conv, control = control)
	if(control$algorithm == "ccd"){
		out <- list(call = NULL, formula = NULL, family = NULL, unpenalized = unpenalized, np = np, beta = beta, 
					phi = object$phi[1:np], ru = ru, dev = object$dev[1:np], nnonzero = object$nnonzero[1:np], 
					g = object$g_seq[1:np], X = Matrix(object$X), y = object$y, w = object$w, action = action, 
					conv = object$conv, control = control)
	}
	class(out) <- "dglars"
	out
}
