plot.dglars <- function(x, type = c("both", "AIC", "BIC"), ...){
    type <- match.arg(type)
    dots <- list(...)
    n <- dim(x$X)[1]
    logn <- log(n)
    beta <- x$beta
    g <- x$g
    g.action <- g[x$action != ""]
    out <- if(type == "BIC") BIC(x, ...)
            else AIC(x, ...)
    loglik <- out$loglik
    comp <- out$comp
    if(type == "both") gof_name <- c("AIC", "BIC")
        else gof_name <- out$type
    setting <- switch(type,
                AIC = list(k = out$k, pch = 0),
                BIC = list(k = logn, pch = 4),
                both = list(k = c(out$k, logn), pch = c(0, 4))
                )
    gof <-  sapply(setting$k, function(k) - 2 * loglik + k * comp)
    g.gof <- g[apply(gof, 2, which.min)]
    xlab <- expression(gamma)
    xlim <- range(g)
    ylim <- range(gof)
    xlim[2] <- xlim[2] + 0.20 * (xlim[2] - xlim[1])
    g.gof <- apply(t(g.gof), 2, rep, 5)
    min.gof <- seq(from = ylim[1], to = ylim[2], length = 5)
    # Model Selection
    matplot(g, gof, type = "n", xlab = xlab, ylab = "", xlim = xlim ,main = "Model Selection Criterion", axes = F)
    abline(v = g.action, lty = 2, col = 8)
    matlines(g.gof, min.gof, pch = setting$pch, lty = 2, col = 2, type = "b")
    matpoints(g, gof, type = "o", pch = setting$pch, lty = 2, col = 1)
    legend("topright", legend = gof_name, pch = setting$pch, col = 1, bty = "n")
    axis(1)
    axis(2)
    op <- par(ask = dev.interactive())
    #Coefficients Path
    ylim <- range(beta[-1,])
    min.gof <- seq(from = ylim[1], to = ylim[2], length = 5)
    matplot(g,t(beta[-1,]), type="n", xlab = xlab, ylab = "Regression Coefficients", xlim = xlim, main = "Coefficients Path", axes = F)
    abline(v = g.action, lty = 2, col = 8)
    matlines(g.gof, min.gof, pch = setting$pch, lty = 2, col = 2, type = "b")
    matpoints(g, t(beta[-1,]), col = 1, type = "l", lty = 1)
    legend("topright", legend = gof_name, pch = setting$pch, col = 2, bty = "n")
    axis(1)
    axis(2)
    op <- par(ask=dev.interactive())
    #Rao Score Path
    ru <- t(x$ru)
    ylim <- range(ru)
    min.gof <- seq(from = ylim[1], to = ylim[2], length = 5)
    matplot(g, ru, type = "n",xlab = xlab, ylab = "Rao Score Statistics", xlim = xlim, main = "Rao Score Path", axes = F)
    abline(v = g.action,lty = 2,col = 8)
    matlines(g.gof, min.gof, pch = setting$pch, lty = 2, col = 2, type = "b")
    matpoints(g, ru, col = 1,type = "l", lty = 1)
    legend("topright", legend = gof_name, pch = setting$pch, col = 2, bty = "n")
    axis(1)
    axis(2)
    par(op)
}
