grcv <- function(object, type = c("BIC", "AIC"), nit = 10L, trace = FALSE, control = list(), ...) {
    if(class(object) != "dglars") stop("the general refitted cross-Validation estimator is avaiable only for objects of class 'dglars'")
    if(is.element(object$family$family, c("poisson", "binomial"))) {
        phih <- 1
        names(phih) <- "phi"
        return(phih)
    }
    type <- match.arg(type)
    if(!is.logical(trace)) stop("'trace' is not an object of type 'logical'")
    X <- object$X
    y <- object$y
    id <- sample(rep.int(seq_len(2L), times = round(dim(object$X)[1L] / 2)))
    phih <- double(nit)
    for(j in seq_len(nit)) {
        if(trace) cat("Iteration number:", j, "\n")
        for(i in seq_len(2L)) {
            if(trace) cat("Cycle ", i, "\n")
            #################################
            # Step 1: splitting the dataset #
            #################################
            if(trace) cat("Step 1: splitting the dataset... ")
            # training set
            y1 <- y[id == i]
            X1 <- as.matrix(X[id == i, , drop = FALSE])
            # testing set
            y2 <- y[id != i]
            X2 <- as.matrix(X[id != i, , drop = FALSE])
            if(trace) cat("completed!\n")
            ##############################
            # Step 2: variable selection #
            ##############################
            if(trace) cat("Step 2: variable selection... ")
            data.train <- data.frame(yt = y1, Xt = X1)
            frml <- as.formula(paste("yt ~ ", paste(names(data.train)[-1L], collapse = " + ")))
            out_update <- update(object, formula. = frml, data = data.train, control = control)
            if(out_update$conv != 0) stop("error in Step 2: dglars does not converge")
            out_gof <- switch(type,
                            "AIC" = AIC(out_update, phi = "pearson", ...),
                            "BIC" = BIC(out_update, phi = "pearson", ...))
            bhat <- out_update$beta[, which.min(out_gof$val), drop = TRUE]
            A <- which(abs(bhat[-1L]) > 0)
            if(trace) cat("completed!\n")
            ############################################
            # Step 3: computing MLE using the test set #
            ############################################
            if(trace) cat("Step 3: copmuting mle... ")
            if(length(A) == 0L) out_glm <- glm(y2 ~ 1, family = object$family)
            else suppressWarnings(out_glm <- try(glm(y2 ~ X2[, A, drop = FALSE], family = object$family, start = bhat[abs(bhat) > 0]), silent = TRUE))
            if(class(out_glm)[1L] != "try-error") rp <- residuals.glm(out_glm, type = "pearson")
            else {
                data.test <- data.frame(yt = y2, Xt = X2[, A, drop = FALSE])
                frml <- as.formula(paste("yt ~ ", paste(names(data.test)[-1L], collapse = " + ")))
                out_update <- update(object, formula. = frml, data = data.test, control = list(g0 = 1.0e-4))
                if(out_update$conv != 0) stop("error in Step 3: I can not compute the maximum likelihood estimates.")
                np <- out_update$np
                mu <- as.vector(predict.dglars(out_update, type = "mu", g = out_update$g[np]))
                V <- object$family$variance(mu)
                rp <- (y2 - mu) / sqrt(V)
            }
            if(trace) cat("completed!\n\n")
            ###############################################
            # Step 4: estimating the dispersion parameter #
            ###############################################
            phih[j] <- phih[j] + 0.5 * sum(rp^2) / (dim(X2)[1L] - length(A) - 1)
        }
    }
    out <- median(phih[is.finite(phih)])
    out
}
