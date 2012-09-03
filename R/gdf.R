gdf <- function(X,y,beta=NULL,control=list()){
	if(!is.list(control))
	stop("control is not a list")
	if(is.null(beta)){
		out.dglars <- dglars(X,y,family="binomial",control=control)
		if(out.dglars$conv!=0) stop("")
		beta <- out.dglars$beta
	} else {
		if(is.vector(beta)) beta <- as.matrix(beta)
	}
	beta_dim <- dim(beta)
	p <- beta_dim[1]-1
	np <- beta_dim[2]
	if(dim(X)[2]!=p)
	stop("dim(beta)[1]!=p")
	gdf_v <- vector(length=np)
	out.glm <- glm(y~X,family="binomial")
	if(!out.glm$converged){
		warning("complexity was set equal to 'df'")
		gdf_v <- colSums(abs(beta)>0)
	} else {
		mu <- out.glm$fit
		V <- binomial()$variance(mu)
		for(i in 1:np){
			bh <- beta[,i,drop=TRUE]
			A <- which(abs(bh[-1])>0)
			if(length(A)==0){
				Xa <- rep(1,length(y))
				eta <- rep(bh[1],length(y))
			} else {
				Xa <- cbind(1,X[,A])
				ba <- bh[c(1,A+1)]
				eta <- drop(tcrossprod(ba,Xa))
			}
			mu_g <- binomial()$linkinv(eta)
			V_g <- binomial()$variance(mu_g)
			Ib <- crossprod(sqrt(V)*Xa)
			invIb_g <- solve(crossprod(sqrt(V_g)*Xa))
			gdf_v[i] <- drop(crossprod(as.vector(invIb_g),as.vector(Ib)))			
		}
	}
	gdf_v
}