dglars_pc <- function(n, p, X, y, nup, up, w, fml, linkid, setting){
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(y) <- "double"
	if(fml == "binomial"){
		if(is.matrix(y)){
			mi <- y[, 1] + y[, 2]
			y <- y[, 1]
		} else mi <- rep(1, n)
		storage.mode(mi) <- "double"
	}
    storage.mode(nup) <- "integer"
	storage.mode(w) <- "double"
	storage.mode(linkid) <- "integer"
	b <- double((p+1)*setting$np)
	phi <- double(setting$np)
	ru <- double(p*setting$np)
	dev <- double(setting$np)
	A <- seq.int(from = 1, to = p)
	if(nup!=0) A <- c(A[up], A[-up])
	storage.mode(A) <- "integer"
	storage.mode(setting$nv) <- "integer"
	nav <- as.integer(0)
	nnonzero <- integer(setting$np)
	g_seq <- double(setting$np)	
	storage.mode(setting$g0) <- "double"
	g_hat <- as.double(2)
	storage.mode(setting$dg_max) <- "double"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$ncrct) <- "integer"
	storage.mode(setting$cf) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	mthd <- ifelse(setting$method=="dgLASSO",1L,0L)
	conv <- 1L
	if(fml == "gaussian"){
		if(linkid == 1)
			fit=.Fortran(C_pc_gaussian_c,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,ru=ru,dev=dev,A=A,
					nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,g0=setting$g0,g_hat=g_hat,
					eps=setting$eps,np=setting$np,conv=conv)	
		else fit=.Fortran(C_pc_gaussian_g,linkid=linkid,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,
					ru=ru,dev=dev,A=A,nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,
					g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,
					ncrct=setting$ncrct,cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
	}
	if(fml == "poisson"){
		if(linkid == 2)
			fit=.Fortran(C_pc_pois_c,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,ru=ru,dev=dev,A=A,
					nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,g0=setting$g0,g_hat=g_hat,
					dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
					cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
		else fit=.Fortran(C_pc_pois_g,linkid=linkid,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,
					ru=ru,dev=dev,A=A,nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,
					g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,
					ncrct=setting$ncrct,cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
	}	
	if(fml == "binomial"){
		if(linkid == 8)
			fit=.Fortran(C_pc_bin_c,n=n,p=p,X=X,y=y,mi=mi,nup=nup,w=w,b=b,phi=phi,ru=ru,dev=dev,A=A,
					nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,g0=setting$g0,g_hat=g_hat,
					dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
					cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
		else fit=.Fortran(C_pc_bin_g,linkid=linkid,n=n,p=p,X=X,y=y,mi=mi,nup=nup,w=w,b=b,phi=phi,
					ru=ru,dev=dev,A=A,nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,
					g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,
					ncrct=setting$ncrct,cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
	}	
	if(fml == "Gamma"){
		if(linkid == 3) {
			fit=.Fortran(C_pc_gamma_c,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,ru=ru,dev=dev,A=A,
					nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,g0=setting$g0,g_hat=g_hat,
					dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
					cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
			fit$b <- -fit$b
			fit$ru <- -fit$ru
		} else fit=.Fortran(C_pc_gamma_g,linkid=linkid,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,
					ru=ru,dev=dev,A=A,nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,
					g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,
					ncrct=setting$ncrct,cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
	}
	if(fml == "inverse.gaussian"){
		if(linkid == 9) {
			fit=.Fortran(C_pc_invgaus_c,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,ru=ru,dev=dev,A=A,
					nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,g0=setting$g0,g_hat=g_hat,
					dg_max=setting$dg_max,eps=setting$eps,np=setting$np,ncrct=setting$ncrct,
					cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
			fit$b <- -2*fit$b
			fit$ru <- -2*fit$ru
		} else fit=.Fortran(C_pc_invgaus_g,linkid=linkid,n=n,p=p,X=X,y=y,nup=nup,w=w,b=b,phi=phi,
					ru=ru,dev=dev,A=A, nv=setting$nv,nav=nav,nnonzero=nnonzero,g_seq=g_seq,mthd=mthd,
					g0=setting$g0,g_hat=g_hat,dg_max=setting$dg_max,eps=setting$eps,np=setting$np,
					ncrct=setting$ncrct,cf=setting$cf,NReps=setting$NReps,nNR=setting$nNR,conv=conv)
	}	
	fit <- make_dglars(fit, setting)
	fit
}
