dglars_ccd <- function(n, p, X, y, g, nup, up, w, fml, linkid, setting){
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
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$g0) <- "double"
	g_hat <- as.double(2)
	storage.mode(setting$nccd) <- "integer"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	mthd <- ifelse(setting$method == "dgLASSO", 1L, 0L)
	b <- double((p + 1) * setting$np)
	phi <- double(setting$np)
	ru <- double(p*setting$np)
	dev <- double(setting$np)
	if(missing(g)) g_seq <- double(setting$np)
	else {
		g_seq <- g
		storage.mode(g_seq) <- "double"
	}
	A <- seq.int(from = 1, to = p)
	if(nup!=0) A <- c(A[up], A[-up])
	A <- c(0, A)
	storage.mode(A) <- "integer"
	nnonzero <- integer(setting$np)
	nav <- as.integer(0)
	conv <- integer(1)
	if(fml %in% c("gaussian", "Gamma", "inverse.gaussian"))
		stop("'ccd' algorithm is not yet available; please, fit the model by using 'pc' algorithm ")
	if(fml == "binomial"){
		if(linkid == 8)
			fit = .Fortran(C_ccd_bin_c,n=n,p=p,X=X,y=y,mi=mi,nup=nup,w=w,np=setting$np,
							g0=setting$g0,g_hat=g_hat,nstp=setting$nccd,eps=setting$eps,
							NReps=setting$NReps,nNR=setting$nNR,mthd=mthd,b=b,phi=phi,
							ru=ru,dev=dev,g_seq=g_seq,A=A,nnonzero=nnonzero,nav=nav,
							conv=conv)
		else
			stop("'ccd' algorithm is not yet available; please, fit the model by using 'pc' algorithm ")
	}
	if(fml == "poisson"){
		if(linkid == 2)
			fit = .Fortran(C_ccd_pois_c,n=n,p=p,X=X,y=y,nup=nup,w=w,np=setting$np,
							g0=setting$g0,g_hat=g_hat,nstp=setting$nccd,eps=setting$eps,
							NReps=setting$NReps,nNR=setting$nNR,mthd=mthd,b=b,phi=phi,
							ru=ru,dev=dev,g_seq=g_seq,A=A,nnonzero=nnonzero,nav=nav,
							conv=conv)
		else
			stop("'ccd' algorithm is not yet available; please, fit the model by using 'pc' algorithm ")
	}	
	fit$A <- fit$A[-1]
	fit <- make_dglars(fit, setting)
	fit
}
