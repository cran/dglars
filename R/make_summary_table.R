make_summary_table <- function(x,complexity){
	n <- dim(x$X)[1]
	np <- x$np
	action <- x$action
	b <- x$beta
	dev <- x$dev
	g <- x$g
	if(complexity=="df")	compl <- x$df
	else compl <- gdf(X=x$X,y=x$y,beta=x$beta)
	aic <- dev + 2 * compl
	rank.aic <- rank(aic)
	best <- rank.aic==1
	b.aic <- b[,best]
	b.aic <- b.aic[abs(b.aic)>0]
	mark <- rep("  ",np)
	mark[best] <- "<-"
	rank.aic <- paste(rank.aic,mark)
	bic <- dev + log(n) * compl
	rank.bic <- rank(bic)
	best <- rank.bic==1
	b.bic <- b[,best]
	b.bic <- b.bic[abs(b.bic)>0]
	mark <- rep("  ",np)
	mark[best] <- "->"
	rank.bic <- paste(mark,rank.bic)
	tbl <- data.frame(Sequence=action,g=g,Dev=dev,Complexity=compl,AIC=aic,Rank.AIC=rank.aic,Rank.BIC=rank.bic,BIC=bic)
	list(table=tbl,b.aic=b.aic,b.bic=b.bic)
}
