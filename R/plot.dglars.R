plot.dglars <- function(x,...){
	n <- dim(x$X)[1]
	beta <- x$beta
	dev <- x$dev
	g <- x$g
	g.action <- g[x$action!=""]
	df <- x$df
	aic <- dev + 2 * df
	bic <- dev + log(n) * df
	g.aic <- g[which.min(aic)]
	g.bic <- g[which.min(bic)]
	plot(g,bic,xlab=expression(gamma),ylab="BIC",type="n")
	abline(v=g.action,lty=2,col=8)
	abline(v=g.bic,lty=2,col=2)
	points(g,bic,xlab=expression(gamma),ylab="BIC",type="o",pch=20,lty=2,...)
	op <- par(ask=dev.interactive())
	plot(g,aic,xlab=expression(gamma),ylab="AIC",type="n")
	abline(v=g.action,lty=2,col=8)
	abline(v=g.aic,lty=2,col=2)
	points(g,aic,xlab=expression(gamma),ylab="AIC",type="o",pch=20,lty=2,...)
	matplot(g,t(beta[-1,]),col=1,type="n",xlab=expression(gamma),ylab="coef")
	abline(v=g.action,lty=2,col=8)
	abline(v=c(g.aic,g.bic),lty=2,col=2)
	matpoints(g,t(beta[-1,]),col=1,type="l",lty=1,xlab=expression(gamma),ylab="coef",...)
	axis(3,c(g.aic,g.bic),c("AIC","BIC"))
	if(x$control$algorithm=="pc"){
		ru <- x$ru
		matplot(g,t(ru),col=1,type="n",xlab=expression(gamma),ylab="ru")
		abline(v=g.action,lty=2,col=8)
		abline(v=c(g.aic,g.bic),lty=2,col=2)
		matpoints(g,t(ru),col=1,type="l",lty=1,xlab=expression(gamma),ylab="ru",...)
		if(length(dev)!=1) axis(3,c(g.aic,g.bic),c("AIC","BIC"))
	}
	par(op)
}