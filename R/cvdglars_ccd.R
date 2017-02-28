cvdglars_ccd <- function(n, p, X, y, g_seq, nup, up, w, familyid, linkid, setting){
	g_cv <- (setting$ng - 1):0 / (setting$ng - 1)
	storage.mode(familyid) <- "integer"
	storage.mode(linkid) <- "integer"
	storage.mode(n) <- "integer"
	storage.mode(p) <- "integer"
	storage.mode(X) <- "double"
	storage.mode(y) <- "double"
	if(familyid == 2 & is.matrix(y)){
		mi <- y[, 1] + y[, 2]
		y <- y[, 1]
	} else mi <- rep(1, n)
	storage.mode(mi) <- "double"
	storage.mode(nup) <- "integer"
	A <- seq.int(from = 1, to = p)
	if(nup != 0) A <- c(A[up], A[-up])
	A <- c(0, A)
	storage.mode(A) <- "integer"
	storage.mode(w) <- "double"
	storage.mode(setting$foldid) <- "integer"
	storage.mode(setting$nfold) <- "integer"
	storage.mode(setting$np) <- "integer"
	if(missing(g_seq)) g_seq <- double(setting$np)
	else storage.mode(g_seq) <- "double"
	storage.mode(setting$g0) <- "double"
	storage.mode(setting$ng) <- "integer"
	storage.mode(g_cv) <- "double"
	storage.mode(setting$nccd) <- "integer"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	g_hat <- double(1)
	mthd <- ifelse(setting$method == "dgLASSO", 1L, 0L)
	b <- double(p + 1)
	phi <- double(1)
	dev_m <- double(setting$ng)
	dev_v <- double(setting$ng)
	conv <- integer(setting$nfold + 1)
	fit = .Fortran(C_ccd_cvdglars, familyid = familyid, linkid = linkid, n = n, p = p, X = X, y = y, mi = mi,
					nup = nup, A = A, w = w, foldid = setting$foldid, nfold = setting$nfold, np = setting$np,
					g_seq = g_seq, g0 = setting$g0, ng = setting$ng, g_cv = g_cv, nstp = setting$nccd, eps = setting$eps,
					NReps = setting$NReps, nNR = setting$nNR, g_hat = g_hat, mthd = mthd, b = b, phi = phi,
					dev_m = dev_m, dev_v = dev_v, conv = conv)
	fit <- make_cvdglars(fit, setting)
	fit
}
