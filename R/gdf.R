gdf <- function(object){
    if(class(object)[1] != "dglars")
        stop("'gdf' works only with objects with class 'dglars'")
    fml <- object$family
    model_used <- paste(fml$family, fml$link, sep = "_")
    canonical <- c("gaussian_identity", "binomial_logit", "poisson_log", "Gamma_inverse", "inverse.gaussian_1/mu^2")
    if(model_used == canonical[1])
        return(object$nnonzero)
    y <- object$y
    X <- object$X
    beta <- object$beta
	beta_dim <- dim(beta)
	np <- beta_dim[2]
    start <- beta[, np]
    try(out_glm <- glm(y ~ as(X, "matrix"), family = fml, start = start), silent = TRUE)
    if(!out_glm$converged){
        warning("'glm does not converge. Complexity is set equal to 'df'")
        return(object$nnonzero)
    }
    if(is.matrix(y)){
        mi <- y[, 1] + y[, 2]
        y <- y[, 1]
    } else mi <- rep.int(1, length(y))
    X <- cbind(1, X)
    A <- predict(object, type = "predictors")
    if(fml$family == "binomial"){
        prb <- out_glm$fitted.values
        prb_A <- predict(object, type = "probability")
        muh <- mi * prb
        muh_A <- mi * prb_A
        sqrt_dmu_dth <- sqrt(mi * fml$variance(prb))
        dmu_dth_A <- apply(prb_A, 2, function(x) mi * fml$variance(x))
    } else {
        muh <- out_glm$fitted.values
        muh_A <- predict(object, type = "mu")
        if(fml$family != "gaussian"){
            sqrt_dmu_dth <- sqrt(fml$variance(muh))
            dmu_dth_A <- apply(muh_A, 2, fml$variance)
        }
    }
    if(model_used %in% canonical){
        sqrt_dmu_dth_A <- sqrt(dmu_dth_A)
        gdf_v <- sapply(1:length(A), function(id){
                            if(A[[id]][1] == 0) a <- 1
                            else a <- c(1, 1 + A[[id]])
                            X_A <- X[, a]
                            I <- crossprod(sqrt_dmu_dth * X_A)
                            J <- crossprod(sqrt_dmu_dth_A[, id] * X_A)
                            M <- solve(J, I)
                            sum(diag(M))
                        })
    } else {
        eta_A <- predict(object, type = "eta")
        r_A <- y - muh_A
        dmu_de_A <- apply(eta_A, 2, fml$mu.eta)
        if(fml$family == "binomial") dmu_de_A <- mi * dmu_de_A
        if(fml$link != "identity"){
            d2mu_de2_fun <- d2mu_de2_mk(fml$link)
            d2mu_de2_A <- apply(eta_A, 2, function(eta) d2mu_de2_fun(eta, mi))
        }
        if(fml$family != "gaussian"){
            dth_dmu_A <- dmu_dth_A^-1
            dth_de_A <- dth_dmu_A * dmu_de_A
            d2th_dmu2_fun <- d2th_dmu2_mk(fml$family)
            d2th_dmu2_A <- apply(muh_A, 2, function(mu) d2th_dmu2_fun(mu, mi))
            d2th_de2_A <- if(fml$link != "identity") d2th_dmu2_A * dmu_de_A^2 + dth_dmu_A * d2mu_de2_A
                            else    d2th_dmu2_A * dmu_de_A^2
            wght1 <- sqrt_dmu_dth * dth_de_A
            wght2 <- dth_de_A * dmu_de_A - d2th_de2_A * r_A
        } else {
            wght1 <- dmu_de_A
            wght2 <- if(fml$link != "identity") dmu_de_A^2 - d2mu_de2_A * r_A
                        else dmu_de_A^2
        }
        gdf_v <- sapply(1:length(A), function(id){
                            if(A[[id]][1] == 0) a <- 1
                            else a <- c(1, 1 + A[[id]])
                            X_A <- X[, a]
                            I <- crossprod(wght1[, id] * X_A)
                            #J <- crossprod(sqrt(wght2[, id]) * X_A)
                            J <- crossprod(X_A, wght2[, id] * X_A)
                            M <- solve(J, I)
                            sum(diag(M))
                        })
    }
    gdf_v
}
