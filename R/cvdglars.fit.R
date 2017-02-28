cvdglars.fit <- function(X, y, family = gaussian, g, unpenalized, b_wght, control = list()){
	this.call <- match.call()
	if(is.data.frame(X)) X <- as.matrix(X)
	if(is.null(colnames(X))) colnames(X) <- paste("X", 1:dim(X)[2], sep = "")
	if(missing(unpenalized)){
		unpenalized <- 0
		nup <- 0
	}
	else{
		if(!is.vector(unpenalized))	stop("''unpenalized is not a vector")
		if(is.list(unpenalized))	stop("'unpenalized' can not be a list")
		if(is.factor(unpenalized))	stop("'unpenalized' can not be a factor")
		if(is.character(unpenalized)) {
			unpenalized_id <- pmatch(unpenalized, colnames(X))
			if(any(is.na(unpenalized_id)))
			stop(gettextf("the following names are not in colnames(X): %s", paste(unpenalized[is.na(unpenalized_id)], collapse = ", ")))
			unpenalized <- sort(unpenalized_id)
		}
		else unpenalized <- sort(unpenalized)
		if(any(abs(unpenalized - round(unpenalized)) > .Machine$double.eps^0.5))
		stop("some element of 'unpenalized' is not an integers")
		if(any(unpenalized <= 0))
		stop("some element of 'unpenalized' is smaller than zero")
		if(any(unpenalized > p))
		stop("some element of 'unpenalized' is greater than the number of columns of the matrix 'X'")
		nup <- length(unpenalized)
	}
	if(missing(b_wght)) b_wght <- double(dim(X)[2] + 1)
	else {
		if(!is.vector(b_wght)) stop("argument 'b_wght' is not a vector")
		if(is.list(b_wght)) stop("argument 'b_wght' can not be a list")
		if(is.factor(b_wght)) stop("argument 'b_wght' can not be a factor")
		if(is.character(b_wght)) stop("vector 'b_wght' can not be a character")
		if(length(b_wght) != (dim(X)[2] + 1))
		stop("the length of the vector 'b_wght' is not equal to ", sQuote(dim(X)[2] + 1))
        if(any(is.nan(b_wght))) b_wght[is.nan(b_wght)] <- 0
		if(all(b_wght[-1] == 0))
		stop("all the entries of the vector 'b_wght[-1]' are equal to zero")
		if(any(b_wght[-1] == 0)){
			if(nup != 0){
				if(any(b_wght[unpenalized + 1] == 0))
				stop("inconsistency between 'b_wght' and 'unpenalized': the weight associated to the unpezalized parameter can not be zero")
			}
			zeroid <- which(b_wght[-1] == 0)
			b_wght <- b_wght[-(zeroid + 1)]
			if(nup != 0) unpenalized <- which(seq.int(1L, dim(X)[2])[-zeroid] %in% unpenalized)
			X <- X[, -zeroid]
		}
	}
	Xdim <- dim(X)
	n <- Xdim[1]
	p <- Xdim[2]
	min_np <- min(n - 1, p)
	if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) family <- do.call(family, args = list(), envir = parent.frame())
	if (is.null(family$family)){
		print(family)
		stop("'family' not recognized")
	}
	family_used <- family$family
	okFamily <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian")
	link_used <- family$link
	okLinks <- c("identity", "log", "inverse", "sqrt", "cloglog", "probit", "cauchit", "logit", "1/mu^2")
	if(!(link_used %in% okLinks))
		stop(gettextf("%s link not recognized", sQuote(link_used)), domain = NA)
	if(family_used == "gaussian" & !is.element(link_used, c("identity","log", "inverse")))
		stop(gettextf("The %s family does not accept the link %s.\n  The accepted links are: 'identity', 'log', 'inverse'", 
			sQuote(family_used), sQuote(link_used)), domain = NA)
	if(family_used == "poisson" & !is.element(link_used, c("log", "identity", "sqrt")))
		stop(gettextf("The %s family does not accept the link %s.\n  The accepted links are: 'log', 'identity', 'sqrt'", 
			sQuote(family_used), sQuote(link_used)), domain = NA)
	if(family_used == "binomial" & !is.element(link_used, c("logit", "probit", "cauchit", "log", "cloglog")))
		stop(gettextf("The %s family does not accept the link %s.\n  The accepted links are: 'logit', 'probit', 'cauchit', 'log', 'cloglog'", 
			sQuote(family_used), sQuote(link_used)), domain = NA)
	if(family_used == "Gamma" & !is.element(link_used, c("inverse", "identity", "log")))
		stop(gettextf("The %s family does not accept the link %s.\n  The accepted links are: 'inverse', 'identity', 'log'", 
			sQuote(family_used), sQuote(link_used)), domain = NA)
	if(family_used == "inverse.gaussian" & !is.element(link_used, c("1/mu^2", "inverse", "identity", "log")))
		stop(gettextf("The %s family does not accept the link %s.\n  The accepted links are: '1/mu^2', 'inverse', 'identity', 'log'", 
			sQuote(family_used), sQuote(link_used)), domain = NA)
	familyid <- pmatch(family_used, okFamily)
	linkid <- pmatch(link_used, okLinks)
	if(is.character(y) & family_used != "binomial")
		stop(gettextf("The %s family does not accept a character as responce variable", sQuote(family_used)), domain = NA)
	if(is.factor(y) & family_used != "binomial")
		stop(gettextf("The %s family does not accept a factor as responce variable", sQuote(family_used)), domain = NA)
	if(is.matrix(y) & family_used != "binomial")
		stop(gettextf("The %s family does not accept a matrix as responce variable", sQuote(family_used)), domain = NA)
	if(is.matrix(y)){
		if(is.data.frame(y))
			y <- as.matrix(y)
		if(is.character(y))
			stop("'y' can not be a character matrix")
			ly <- dim(y)[1]
			if(ly != n)
				stop("the number of the rows of the matrix 'y' is not equal to the number of rows of the matrix 'X'")
			if(dim(y)[2] > 2)
				stop("the number of columns of the matrix 'y' is greater than 2")
			if(dim(y)[2] == 1)
				y <- y[, 1, drop = FALSE]
        ny <- y
	} else {
		ly <- length(y)
		if(ly != n)
			stop("the length of the vector 'y' is not equal to the number of rows of the matrix 'X'")
		if(is.character(y))
			y <- factor(y)
		if(is.factor(y)){
			if(nlevels(y) != 2)
				stop("'y' can not be a factor with more than two levels")
			ny <- as.numeric(y) - 1
        } else ny <- y
	}
	if(family_used %in% c("binomial", "poisson") & any(abs(ny - round(ny)) > .Machine$double.eps^0.5))
		stop(gettextf("each element of 'y' must be integer when family %s is used", sQuote(family_used)), domain = NA)
	if(family_used %in% c("poisson", "Gamma", "inverse.gaussian") & any(ny < 0))
		stop(gettextf("each element of 'y' must be positive when family %s is used", sQuote(family_used)), domain = NA)
	if(family_used == "binomial"){
		if(is.vector(ny)) frq <- ny
		else frq <- ny[, 1] / (ny[, 1] + ny[, 2])
		if(any(frq < 0 | frq > 1))
			stop("some element of 'y' is outside of its range")
	}
	if(!is.list(control))
		stop("control is not a list")
	if(!missing(g)){
		if(!is.vector(g))
		stop("\n'g' is not a vector\n")
		if(is.character(g))
		stop("\n'g' can not be a character vector\n")
		if(is.factor(g))
		stop("\n'g' can not be a factor\n")
		g <- sort(g, decreasing = TRUE)
		if(any(g < 0))
		stop("\nsome value in 'g' is negative\n")
		if("np" %in% names(as.list(control)))
		stop("\nargument 'np' can not be used with 'g'\n")
		if("g0" %in% names(as.list(control)))
		stop("\nargument 'g0' can not be used with 'g'\n")
	}
	setting <- list(algorithm = "pc", method = "dgLASSO", nfold = 10, foldid = NULL, ng = 100,
					nv = min_np, np = NULL, g0 = ifelse(p < n, 1.0e-06, 0.05), dg_max = 0, nNR = 200,
					NReps = 1.0e-06, ncrct = 50, cf = 0.5, nccd = 1.0e+05, eps = 1.0e-05)
	nmsSetting <- names(setting)
	setting[(nms <- names(control))] <- control
	if(length(noNms <- nms[!nms %in% nmsSetting]))
		warning("unknow names in control: ",paste(noNms,collapse =", "))
	if(!setting$algorithm %in% c("pc","ccd"))
		stop("'algorithm' should be one of \"pc\", \"ccd\"")
	if(!missing(g) & setting$algorithm == "pc"){
		warning("\nargument 'g' can be used only with the cyclic coordinate descent algorithm\ncontrol argument 'algorithm' automatically set to \"ccd\"\n")
		setting$algorithm <- "ccd"
	}
	if(!missing(g)){
		setting$np <- length(g)
		setting$g0 <- g[setting$np]
	}
	if(!setting$method %in% c("dgLASSO","dgLARS"))
		stop("'method' should be one of \"dgLASSO\",\"dgLARS\"")
	if(setting$nfold < 1 | setting$nfold > n)
		stop("'nfold' should be an integer between 1 and n")
	if(is.null(setting$foldid)) 
		setting$foldid <- sample(n)
	else{
		if(min(setting$foldid) < 1 | max(setting$foldid) > n)
			stop("Read the documentation on 'foldid' for more details")
	}
	if(setting$ng < 1)
		stop("'ng' should be a non-negative integer")
	if(setting$nv < 0 | setting$nv > min_np) 
		stop("'nv' should be an integer between 1 and min(n-1,p)")
	if(is.null(setting$np))
		setting$np <- ifelse(setting$algorithm == "pc", min_np * 50, 100L)
	if(setting$np <= 1)
		stop("'np' should be a non-negative integer. Read the documentation more details")
	if(setting$g0 < 0)
		stop("'g0' should be a non-negative value. Read the documentation more details")
	if(setting$dg_max < 0)
		stop("'dg_max' should be a non-negative value. Read the documentation more details")
	if(setting$eps <= 0)
		stop("'eps' should be a non-negative value. Read the documentation more details")
	if(setting$ncrct <= 0)
		stop("'ncrct' should be a non-negative value")
	if(setting$NReps <= 0)
		stop("'NReps' should be a non-negative value")
	if(setting$nNR <= 0)
		stop("'nNR' should be a non-negative integer")
	if(setting$cf < 0 | setting$cf > 1)
		stop("'cf' should be a value in the interval (0,1)")
	if(setting$nccd <= 0)
		stop("'nccd' should be a non-negative integer")
	if(nup > setting$nv)
		stop("the number of unpenalized estimates can not be greater than nv")
# check for ccd algorithm
	if(setting$algorithm == "ccd"){
		if(family_used %in% okFamily[-c(2, 3)])
			stop(gettextf("'ccd' algorithm is not available for %s family; fit the model using 'pc' algorithm", sQuote(family_used)), domain = NA)
		if(family_used == "binomial" & link_used %in% okLinks[-8])
			stop(gettextf("'ccd' algorithm is not available for %s family with %s link function; fit the model using 'pc' algorithm", sQuote(family_used), sQuote(link_used)), domain = NA)
		if(family_used == "poisson" & link_used %in% okLinks[-2])
			stop(gettextf("'ccd' algorithm is not available for %s family with %s link function; fit the model using 'pc' algorithm", sQuote(family_used), sQuote(link_used)), domain = NA)
	}
	fit=switch(setting$algorithm,
			   pc = cvdglars_pc(n, p, X, ny, nup, unpenalized, b_wght, familyid, linkid, setting),
			   ccd = cvdglars_ccd(n, p, X, ny, g, nup, unpenalized, b_wght, familyid, linkid, setting)
			   )
	fit$call <- this.call
	fit$family <- family
	fit$y <- y
	if(sum(fit$conv[1:setting$nfold])!=0){
		errorfold <- which(fit$conv[1:setting$nfold] != 0)
		msg <- format(c("dglars does not converge in fold:", "error code:"), justify = "right")
		msg <- apply(cbind(msg, rbind(errorfold, fit$conv[errorfold])), 1, paste, collapse = " ")
		msg <- paste(msg, collapse = "\n")
		if(sum(fit$conv[1:setting$nfold] == 0) <= 1L) {
			msg <- paste("\n", msg, "\n\n", "impossible complete cross-validation", sep = "")
			stop(msg)
		} else {
			msg <- paste(msg, paste("cross-validation done using",sum(fit$conv[1:setting$nfold] == 0),"folds"), sep = "\n\n")
			warning(msg)
		}
	}
	if(fit$conv[setting$nfold + 1] != 0)
		stop("\n\nimpossible complete cross-validation\n'dglars' does not coverge in the last step with error code ", fit$conv[setting$nfold + 1])
	fit
}
