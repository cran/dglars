predict.dglars <- function(object, xnew, ynew, g = NULL, type = c("coefficients", "nnonzero", "predictors", "eta", "mu", "probability", "class", "deviance"), ...){
	type <- match.arg(type)
    family_used <- object$family$family
    link_used <- object$family$link
    if(type %in% c("probability", "class") & family_used != "binomial")
        stop("type ", sQuote(type), "can be used only with family 'binomial'")
    if(!missing(xnew) & type %in% c("coefficients", "nnonzero", "predictors"))
        warning("argument 'xnew' is not used with type ", sQuote(type))
    if(!missing(ynew) & type != "deviance")
        warning("argument 'ynew' is not used with type ", sQuote(type))
    if(type == "deviance" & ((missing(xnew) & !missing(ynew)) | (!missing(xnew) & missing(ynew)))) {
        if(missing(xnew))
            stop("'xnew' is required to compute the deviance when 'ynew' is available")
        else
            stop("'ynew' is required to compute the deviance when 'xnew' is available")
    }
    if(!missing(ynew)){
        if(class(ynew) != class(object$y))
        stop("'y' of class ", sQuote(class(object$y)), " while the class of the argument 'ynew' is ", sQuote(class(ynew)))
    }
    X <- object$X
    n <- dim(X)[1]
    p <- dim(X)[2]
    y <- object$y
    if(is.matrix(y)){
        if(!is.null(colnames(y))) name_class <- colnames(y)
        else name_class <- c("Y1", "Y2")
        mi <- y[, 1] + y[, 2]
        y <- y[, 1]
    } else {
        mi <- rep(1, n)
        if(is.factor(y)) {
            name_class <- levels(y)[2:1]
            y <- as.numeric(y) - 1
        } else name_class <- c(1, 0)
    }
    if(missing(xnew)) xnew <- X
    else {
        if(dim(xnew)[2] != p)
            stop("the number of columns of the matrix 'xnew' is not equal to ", p)
        if(!inherits(xnew, "Matrix"))
            xnew <- Matrix(xnew)
    }
    nnew <- dim(xnew)[1]
    if(missing(ynew)) {
        ynew <- y
        minew <- mi
    }
    else{
        if(is.character(ynew) & family_used != "binomial")
            stop(gettextf("The %s family does not accept a character as responce variable", sQuote(family_used)), domain = NA)
        if(is.factor(ynew) & family_used != "binomial")
            stop(gettextf("The %s family does not accept a factor as responce variable", sQuote(family_used)), domain = NA)
        if(is.matrix(ynew) & family_used != "binomial")
            stop(gettextf("The %s family does not accept a matrix as responce variable", sQuote(family_used)), domain = NA)
        if(is.matrix(ynew)){
            if(is.data.frame(y))
                y <- as.matrix(y)
            if(is.character(y))
                stop("'ynew' can not be a character matrix")
            if(dim(ynew)[1] != nnew)
                stop("the matrices 'xnew' and 'ynew' have a different number of rows")
            if(dim(ynew)[2] == 1){
                minew <- rep(1, length(ynew))
                ynew <- ynew[, 1]
            } else{
                if(dim(ynew)[2] > 2)
                    stop("the number of columns of the matrix 'ynew' is greater than 2")
                minew <- ynew[, 1] + ynew[, 2]
                ynew <- ynew[, 1]
            }
        } else {
            if(length(ynew) != nnew)
                stop("the length of the vector 'ynew' is different from the number of rows of the matrix 'xnew'")
            minew <- rep(1, length(ynew))
        }
        if(is.character(ynew))
            ynew <- factor(ynew)
        if(is.factor(ynew)){
            if(nlevels(ynew) != 2)
                stop("'ynew' can not be a factor with more than two levels")
            ynew <- as.numeric(ynew) - 1
        }
        if(family_used %in% c("binomial", "poisson") & any(abs(ynew - round(ynew)) > .Machine$double.eps^0.5))
            stop(gettextf("each element of 'ynew' must be integer when family %s is used", sQuote(family_used)), domain = NA)
        if(family_used %in% c("poisson", "Gamma", "inverse.gaussian") & any(ynew < 0))
            stop(gettextf("each element of 'ynew' must be positive when family %s is used", sQuote(family_used)), domain = NA)
        if(family_used == "binomial"){
            if(is.vector(ynew)) frq <- ynew
            else frq <- ynew / minew
            if(any(frq < 0 | frq > 1))
                stop("some element of 'ynew' is outside of its range")
        }
    }
    beta <- object$beta
    if(is.null(g)){
        g <- object$g
        ng <- length(g)
        phi <- object$phi
        nnonzero <- object$nnonzero
    } else {
		if(class(object)[1] != "dglars")
		stop("argument 'g' is not used for objects with class 'cvdglars'")
        okFamily <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian")
        familyid <- pmatch(family_used, okFamily)
        okLinks <- c("identity", "log", "inverse", "sqrt", "cloglog", "probit", "cauchit", "logit", "1/mu^2")
        linkid <- pmatch(link_used, okLinks)
        np <- object$np
        ng <- length(g)
        g_seq <- object$g
        if(min(g) < g_seq[np] | max(g) > g_seq[1])
        stop("Some element of 'g' is outside of the interval [", round(g_seq[np], 6), ", ", round(g_seq[1], 6), "]")
        coef <- matrix(0, nrow = p + 1, ncol = ng, dimnames = list(rownames(beta)))
        phi <- rep(1, ng)
        storage.mode(familyid) <- "integer"
        storage.mode(linkid) <- "integer"
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        X <- as(X, "matrix")
        storage.mode(X) <- "double"
        storage.mode(y) <- "double"
        storage.mode(mi) <- "double"
        storage.mode(np) <- "integer"
        beta <- as(beta, "matrix")
        storage.mode(beta) <- "double"
        storage.mode(g_seq) <- "double"
        storage.mode(ng) <- "integer"
        storage.mode(g) <- "double"
        storage.mode(coef) <- "double"
        storage.mode(phi) <- "double"
        out = .Fortran(C_predict, familyid = familyid, linkid = linkid, n = n, p = p,
                X = X, y = y, mi = mi, np = np, b = beta, g_seq = g_seq, ng = ng, g = g,
                coef = coef, phi = phi)
        beta <- as(out$coef, "dgCMatrix")
        phi <- out$phi
        nnonzero <- NULL
    }
    if(type == "coefficients"){
        coef <- list(beta = beta, phi = phi)
        return(coef)
    }
    if(type == "nnonzero"){
        if(is.null(nnonzero)) nnonzero <- apply(abs(beta) > 0, 2, sum)
        return(nnonzero)
    }
    if(type == "predictors"){
        id <- seq(p)
        nm <- rownames(beta)[-1]
        names(id) <- nm
        predictors <- apply(abs(beta[-1, , drop = FALSE]) > 0, 2,
                        function(x){
                           if(length(id[x]) > 0) id[x]
                           else {
                               out <- 0
                               names(out) <- "Int."
                               out
                           }
                        })
        if(!is.list(predictors)){
            if(is.vector(predictors)) predictors <- matrix(predictors, nrow = 1)
            predictors <- apply(predictors, 2, list)
            predictors <- lapply(predictors, function(x){x <- unlist(x); names(x) <- nm[x]; x})
        }
        names(predictors) <- g
        return(predictors)
    }
    eta <- cbind(1, xnew) %*% beta
    if(type == "eta") return(eta)
    mu <- object$family$linkinv(as(eta, "matrix"))
    if(type == "mu"){
        if(family_used != "binomial") return(as(mu, "dgeMatrix"))
        else return(as(mi * mu, "dgeMatrix"))
    }
    if(type == "probability" & family_used == "binomial") return(as(mu, "dgeMatrix"))
    if(type == "class" & family_used == "binomial"){
        class <- matrix(0, nrow = nnew, ncol = ng)
        class[mu >= 0.5] <- name_class[1]
        class[mu < 0.5] <- name_class[2]
        if(typeof(class) != "character") class <- Matrix(class)
        return(class)
    }
    dev <- apply(mu, 2, function(muh) sum(object$family$dev(ynew / minew, muh, minew)))
    dev
}
