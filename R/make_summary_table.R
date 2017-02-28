make_summary_table <- function(object, type, ...){
    out_gof <- if(type == "AIC") AIC(object, ...)
                else BIC(object, ...)
    action <- object$action
    g <- out_gof$g
    dev <- object$dev
    dev.ratio <- 1 - dev/dev[1]
    compl <- out_gof$comp
    gof <- out_gof$val
    rank.gof <- rank(gof)
    best <- rank.gof == 1
    b <- object$beta
    b.gof <- b[, best]
    names(b.gof)[1] <- "1"
    b.gof <- b.gof[abs(b.gof) > 0]
    formula.gof <- if(is.null(object$formula)) NULL
                    else {
                        formula_tmp.gof <- as.formula(paste("~", paste(names(b.gof), collapse = " + ")))
                        update.formula(object$formula, formula_tmp.gof)
                    }
    names(b.gof)[1] <- "Int."
    g.gof <- g[best]
    np <- object$np
    mark <- rep("  ", np)
    mark[best] <- "<-"
    rank.gof <- paste(rank.gof, mark)
    tbl <- data.frame(Sequence = action, g = g, Dev.ratio = dev.ratio, Complexity = compl, gof = gof, Rank = rank.gof)
    names(tbl)[3] <- "%Dev"
    names(tbl)[4] <- out_gof$complexity
    names(tbl)[5] <- out_gof$type
    list(table = tbl, formula.gof = formula.gof, b.gof = b.gof, phi.gof = object$phi[best], g.gof = g.gof, nulldev = dev[1],
        resdev.gof = dev[best], type = out_gof$type, k = out_gof$k, complexity = out_gof$complexity, phi = out_gof$phi)
}
