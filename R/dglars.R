dglars <- function(formula, family = gaussian, g, unpenalized, b_wght, data, subset, contrasts = NULL, control = list()){
	this.call <- match.call()	
	if (missing(data))	data <- environment(formula)
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	if(attr(mt,"intercept") == 0)
        stop("Models without intercept are not allowed in this version of the package")
	y <- model.response(mf, "any")
	X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
            else stop("Model matrix is empty")
	X <- X[, -1, drop = FALSE]
	fit <- dglars.fit(X = X, y = y, family = family, g = g, unpenalized = unpenalized, b_wght = b_wght, control = control)
    fit$formula <- formula
	fit$call <- this.call
	fit
}
