dglars <- function(X,y,family=c("binomial","poisson"),control=list()){
	this.call <- match.call()
	family <- match.arg(family)
	if(is.data.frame(X)) X <- as.matrix(X)
	if(is.null(colnames(X))) colnames(X) <- paste("x",1:dim(X)[2],sep="")
	Xdim <- dim(X)
	n <- Xdim[1]
	p <- Xdim[2]
	min_np <- min(n-1,p)
	if(!is.list(control))
	stop("control is not a list")
	setting <- list(algorithm="pc",method="dgLASSO",nv=min_np,np=NULL,g0=ifelse(p<n,1.0e-04,0.05),
					dg_max=0,nNR=50,NReps=1.0e-06,ncrct=50,cf=0.5,nccd=1.0e+05,eps=1.0e-05)
	nmsSetting <- names(setting)
	setting[(nms <- names(control))] <- control
	if(length(noNms <- nms[!nms %in% nmsSetting]))
	warning("unknow names in control: ",paste(noNms,collapse =", "))
	if(length(y)!=n)
	stop("length(y)!=dim(X)[1]")
	if(!setting$algorithm %in% c("pc","ccd"))
	stop("read the documentation for 'algorithm' more carefully")
	if(!setting$method %in% c("dgLASSO","dgLAR"))
	stop("read the documentation for 'method' more carefully")
	if(setting$nv<1 | setting$nv>min_np) 
	stop("read the documentation for 'nv' more carefully")
	if(is.null(setting$np))
	setting$np <- ifelse(setting$algorithm=="pc",min_np*50,100L)
	if(setting$np<=0)
	stop("read the documentation for 'np' more carefully")
	if(setting$g0<0)
	stop("read the documentation for 'g0' more carefully")
	if(setting$dg_max<0)
	stop("read the documentation for 'dg_max' more carefully")
	if(setting$eps<=0)
	stop("read the documentation for 'eps' more carefully")
	if(setting$ncrct<=0)
	stop("read the documentation for 'ncrct' more carefully")
	if(setting$NReps<=0)
	stop("read the documentation for 'NReps' more carefully")
	if(setting$nNR<=0)
	stop("read the documentation for 'nNR' more carefully")
	if(setting$cf<0 | setting$cf>1)
	stop("read the documentation for 'cf' more carefully")
	if(setting$nccd<=0)
	stop("read the documentation for 'nccd' more carefully")
	fit=switch(setting$algorithm,
			   pc=dglars_pc(n,p,X,y,family,setting),
			   ccd=dglars_ccd(n,p,X,y,family,setting)
			   )
	fit$call <- this.call
	fit$family <- family
	if(fit$conv!=0) warning("dglars with algorithm ",fit$control$algorithm," does not converge with exit ",fit$conv)
	fit
}