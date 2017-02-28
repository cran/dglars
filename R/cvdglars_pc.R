cvdglars_pc <- function(n, p, X, y, nup, up, w, familyid, linkid, setting){
	g <- (setting$ng - 1):0 / (setting$ng - 1)
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
	storage.mode(w) <- "double"
	A <- seq.int(from = 1, to = p)
	if(nup != 0) A <- c(A[up], A[-up])
	storage.mode(A) <- "integer"
	storage.mode(w) <- "double"
	storage.mode(setting$foldid) <- "integer"
	storage.mode(setting$nfold) <- "integer"
	storage.mode(setting$ng) <- "integer"
	storage.mode(g) <- "double"
	b <- double(p + 1)
	phi <- double(1)
	dev_m <- double(setting$ng)
	dev_v <- double(setting$ng)
	g_hat <- double(1)
	storage.mode(setting$nv) <- "integer"
	mthd <- ifelse(setting$method == "dgLASSO", 1L, 0L)
	storage.mode(setting$g0) <- "double"
	storage.mode(setting$dg_max) <- "double"
	storage.mode(setting$eps) <- "double"
	storage.mode(setting$np) <- "integer"
	storage.mode(setting$ncrct) <- "integer"
	storage.mode(setting$cf) <- "double"
	storage.mode(setting$NReps) <- "double"
	storage.mode(setting$nNR) <- "integer"
	conv <- integer(setting$nfold + 1)
	fit = .Fortran(C_pc_cvdglars, familyid = familyid, linkid = linkid, n = n, p = p, X = X, y = y, mi = mi,
										nup = nup, A = A, w = w, foldid = setting$foldid, nfold = setting$nfold,ng = setting$ng,
										g = g, b = b, phi = phi, dev_m = dev_m, dev_v = dev_v, g_hat = g_hat, nv = setting$nv,
										mthd = mthd, g0 = setting$g0, dg_max = setting$dg_max, eps = setting$eps,
										np = setting$np, ncrct = setting$ncrct, cf = setting$cf, NReps = setting$NReps,
										nNR = setting$nNR, conv = conv)
	fit <- make_cvdglars(fit, setting)
	fit
}
