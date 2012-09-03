!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!						MAIN SUBROUTNES USED TO COMPUTE THE dgLARS SOLUTION CURVE
!
!										
!
!	DESCRIPTION:	dglars_pc_b and dglars_pc_p are the main subroutines used to compute the dgLARS solution 
!					curve proposed in Augugliaro L., Mineo A.M., and Wit E.C. (JRSSB, accepted) by the predictor
!					corrector method. dglars_pc_b is specific for the binomial family while dglars_pc_p is
!					specific for the Poisson family.
!
!
!	INPUT:
!				n = sample size
!				p = number of predictors
!				X = (n x p) design matrix
!				y = n dimentional responce vector
!				nv = maximum number of predictors included in the final model
!				np = maximum number of the points of the solution curve
!				mthd = method of dglars used to compute the solution curve
!						0 => dgLAR
!						1 => dgLASSO
!				g0 = minimum value of the gamma parameter
!				g_hat = internal flag
!				dg_max = maximum step size allowed
!				eps =
!				ncrct = maximum number of trials in the corrector step
!				cf = corrector factor
!				NReps =
!				nNR =  maximum number of trials for the Newton-Raphson algorithm
!
!
!	OUTPUT:
!				np = number of the points of the solution curve
!				b = (nv+1 x np) coefficient matrix used to store the solution curve
!				ru = (p x np) matrix used to store the path fo the Rao's score test statistics
!				dev = np dimentional vector of the deviance
!				A = p dimentional vector used to identify the sequence of predictors included in the model
!				nav = number of predictor included in the final model
!				df = np dimentional vector of the non-zero coefficients
!				g_seq = np dimentional vector of the gamma values used to compute the solution curve
!				conv = integer used to code the warnings and the errors, i.e.
!						0 => convergence of the algorithm is been achieved
!						1 => error in predictor step
!						2 => error in corrector step
!						3 => maximum number of iterations is been reached
!						4 => error in dynamic allocation memory
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine dglars_pc_b(n,p,X,y,b,ru,dev,A,nv,nav,df,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	double precision,	parameter	:: dzero=1.0e-8,zero=1.0e-5
	integer	::	n,p,nv,nav,np,ncrct,nNR,conv,A(p),df(np),mthd,ai,i,j,k,nstp,nsbstp
	double precision	::	X(n,p),y(n),dev(np),g_seq(np),eps,b(0:p,np),ru(p,np),dg_max,g0,g_hat,NReps,cf,dg,g
	double precision,	dimension(:), allocatable	::	eta,mu,ruv,dmu_de,sqrt_ib,ba,dba,ba_crct,dabsruac
	double precision,	dimension(:,:), allocatable ::	X2
	logical	:: test(3),action,final
	allocate(eta(1:n),mu(1:n),ruv(1:p),dmu_de(1:n),sqrt_ib(1:p),X2(1:n,1:p),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	final= .false.
	nstp=1
	nav=0
	dg=0.
	X2=X*X
	mu=sum(y)/n
	b(0,nstp)=log(mu(1)/(1-mu(1)))
	if(b(0,nstp).eq.0.) then
		df(nstp)=nav
	else
		df(nstp)=nav+1
	end if
	dev(nstp)=-2*(sum(y)*log(mu(1))+((n-sum(y))*log(1-mu(1))))
	eta=b(0,nstp)
	dmu_de=mu*(1-mu)
	forall(j=1:p)	sqrt_ib(j)=sqrt(dmu_de(1)*sum(X2(:,j)))
	call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
	ru(:,nstp)=ruv
	ai=maxloc(abs(ruv),dim=1)
	call shift_A(p,A,nav,ai,1)
	g=abs(ruv(A(1)))
	g_seq(nstp)=g
	if(g.le.g0.or.g_hat.eq.1.) then
		np=nstp
		return
	end if	
	if(g_hat.ne.2.) g0=g0+g_hat*(g-g0)
	nav=nav+1
	action=.true.
	do k=2,np
		test=.false.
		if(action) then
			allocate(ba(0:nav),dba(0:nav),ba_crct(0:nav),dabsruac(1:(p-nav)),stat=conv)
			if(conv.ne.0) then
				conv=4
				np=nstp
				return
			end if
			action=.false.
		end if
		ba(0)=b(0,nstp)
		ba(1:nav)=b(A(1:nav),nstp)
		call prd_b(mthd,g,g0,n,p,X,X2,A,nav,ba,mu,dmu_de,sqrt_ib,ruv,dg_max,dba,dg,conv,ai,final)
		if(conv.ne.0) then
			np=nstp
			return
		end if
		do nsbstp=1,ncrct
			call crct_b(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),y,ba,dba,g,dg,ruv(A(1:nav)),NReps,nNR,mu,dmu_de,ba_crct,conv)
			if(conv.eq.4) then
				np=nstp
				return
			end if
			if(conv.ne.0)	then
				call mu_mk_b(n,eta,mu)
				call dmu_de_mk_b(n,mu,dmu_de)
				conv=0
				dg=dg*cf
				if(dg.le.dzero) then
					conv=2
					np=nstp
					return
				end if
			else
				call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
				call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
				if(.not.final) then
					dabsruac=abs(ruv(A(nav+1:p)))-g+dg
					test(1)=any(dabsruac.gt.eps)
				else
					test(1)=.false.
				end if
				if(mthd.eq.1) then
					test(2)=any(ba_crct(1:nav)*ruv(A(1:nav)).lt.0.)
				else
					test(2)=.false.
				end if
				if(test(1).or.test(2)) then
					call mu_mk_b(n,eta,mu)
					call dmu_de_mk_b(n,mu,dmu_de)
					call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
					call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
					dg=dg*cf
					if(dg.le.dzero) exit
				else
					exit
				end if
			end if
		end do
		if(nsbstp.ge.ncrct) then
			conv=2
			np=nstp
			return
		end if
		test(3)=abs(abs(ruv(A(1)))-g0).le.eps
		if(.not.test(1).and..not.test(2)) then		
			nstp=nstp+1
			b(0,nstp)=ba_crct(0)
			b(A(1:nav),nstp)=ba_crct(1:nav)
			call eta_mk(n,nav,X(:,A(1:nav)),ba_crct,eta)
			call mu_mk_b(n,eta,mu)
			call dmu_de_mk_b(n,mu,dmu_de)
			call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
			call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
			call deviance_b(n,y,mu,dev(nstp))
			df(nstp)=nav+1
			ru(:,nstp)=ruv
			g=abs(ruv(A(1)))
			g_seq(nstp)=g		
			if(minval(abs(dabsruac)).le.eps .and..not.final.and. ai.gt.0) then
				ai=minloc(abs(dabsruac),dim=1)
				call shift_A(p,A,nav,ai,1)
				nav=nav+1
				action=.true.
				if(nav.ge.nv) final=.true.
			end if
			if(any(abs(ba_crct(1:nav)).le.zero).and.mthd.eq.1 .and. ai.lt.0) then
				ai=abs(ai)
				ruv(A(ai))=g
				call shift_A(p,A,nav,ai,-1)
				nav=nav-1
				if(abs(ba_crct(ai)).eq.0.)	df(nstp)=nav+1
				action=.true.
				final=.false.
			end if
		else
			action=.true.
			j=p-nav
			i=nav
			if(test(1)) then
				do ai=1,j
					if(dabsruac(ai).ge.eps) then
						call shift_A(p,A,nav,ai+i-nav,1)
						nav=nav+1
						if(nav.ge.nv) then
							final=.true.
							exit
						end if
					end if
				end do
			end if
			if(test(2)) then
				do ai=1,i
					if(ba_crct(ai)*ruv(A(ai)).lt.0) then
						call shift_A(p,A,nav,ai,-1)
						nav=nav-1
						final=.false.
					end if
				end do
			end if
		end if
		if(action) then
			deallocate(ba,dba,ba_crct,dabsruac,stat=conv)
			if(conv.ne.0) then
				conv=4
				np=nstp
				return
			end if
		end if
		if(test(3)) exit
	end do
	if(k.ge.np) then
		conv=3
		np=nstp
		return
	end if
	np=nstp
	deallocate(eta,mu,ruv,dmu_de,sqrt_ib,X2,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine dglars_pc_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine dglars_pc_p(n,p,X,y,b,ru,dev,A,nv,nav,df,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	double precision,	parameter	:: dzero=1.0e-8,zero=1.0e-5
	integer	::	n,p,nv,nav,np,ncrct,nNR,conv,A(p),df(np),mthd,ai,i,j,k,nstp,nsbstp
	double precision	::	X(n,p),y(n),dev(np),g_seq(np),eps,b(0:p,np),ru(p,np),dg_max,g0,g_hat,NReps,cf,dg,g
	double precision,	dimension(:), allocatable	::	eta,mu,ruv,dmu_de,sqrt_ib,ba,dba,ba_crct,dabsruac
	double precision,	dimension(:,:), allocatable ::	X2
	logical	:: test(3),action,final
	allocate(eta(1:n),mu(1:n),ruv(1:p),dmu_de(1:n),sqrt_ib(1:p),X2(1:n,1:p),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	final= .false.
	nstp=1
	nav=0
	dg=0.
	X2=X*X
	mu=sum(y)/n
	b(0,nstp)=log(mu(1))
	df(nstp)=nav+1
	dev(nstp)=2*sum(y*log(y/mu(1)),mask=y.gt.0.)
	eta=b(0,nstp)
	dmu_de=mu
	forall(j=1:p)	sqrt_ib(j)=sqrt(dmu_de(1)*sum(X2(:,j)))
	call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
	ru(:,nstp)=ruv
	ai=maxloc(abs(ruv),dim=1)
	call shift_A(p,A,nav,ai,1)
	g=abs(ruv(A(1)))
	g_seq(nstp)=g
	if(g.le.g0.or.g_hat.eq.1.) then
		np=nstp
		return
	end if	
	if(g_hat.ne.2.) g0=g0+g_hat*(g-g0)
	nav=nav+1
	action=.true.
	do k=2,np
		test=.false.
		if(action) then
			allocate(ba(0:nav),dba(0:nav),ba_crct(0:nav),dabsruac(1:(p-nav)),stat=conv)
			if(conv.ne.0) then
				conv=4
				np=nstp
				return
			end if
			action=.false.
		end if
		ba(0)=b(0,nstp)
		ba(1:nav)=b(A(1:nav),nstp)
		call prd_p(mthd,g,g0,n,p,X,X2,A,nav,ba,mu,dmu_de,sqrt_ib,ruv,dg_max,dba,dg,conv,ai,final)
		if(conv.ne.0) then
			np=nstp
			return
		end if
		do nsbstp=1,ncrct
			call crct_p(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),y,ba,dba,g,dg,ruv(A(1:nav)),NReps,nNR,mu,dmu_de,ba_crct,conv)
			if(conv.eq.4) then
				np=nstp
				return
			end if
			if(conv.ne.0)	then
				call mu_mk_p(n,eta,mu)
				call dmu_de_mk_p(n,mu,dmu_de)
				conv=0
				dg=dg*cf
				if(dg.le.dzero) then
					conv=2
					np=nstp
					return
				end if
			else
				call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
				call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
				if(.not.final) then
					dabsruac=abs(ruv(A(nav+1:p)))-g+dg
					test(1)=any(dabsruac.gt.eps)
				else
					test(1)=.false.
				end if
				if(mthd.eq.1) then
					test(2)=any(ba_crct(1:nav)*ruv(A(1:nav)).lt.0.)
				else
					test(2)=.false.
				end if
				if(test(1).or.test(2)) then
					call mu_mk_p(n,eta,mu)
					call dmu_de_mk_p(n,mu,dmu_de)
					call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
					call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
					dg=dg*cf
					if(dg.le.dzero) exit
				else
					exit
				end if
			end if
		end do
		if(nsbstp.ge.ncrct) then
			conv=2
			np=nstp
			return
		end if
		test(3)=abs(abs(ruv(A(1)))-g0).le.eps
		if(.not.test(1).and..not.test(2)) then		
			nstp=nstp+1
			b(0,nstp)=ba_crct(0)
			b(A(1:nav),nstp)=ba_crct(1:nav)
			call eta_mk(n,nav,X(:,A(1:nav)),ba_crct,eta)
			call mu_mk_p(n,eta,mu)
			call dmu_de_mk_p(n,mu,dmu_de)
			call sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
			call rao_score(n,p,X,y,mu,sqrt_ib,ruv)
			call deviance_p(n,y,mu,dev(nstp))
			df(nstp)=nav+1
			ru(:,nstp)=ruv
			g=abs(ruv(A(1)))
			g_seq(nstp)=g		
			if(minval(abs(dabsruac)).le.eps .and..not.final.and. ai.gt.0) then
				ai=minloc(abs(dabsruac),dim=1)
				call shift_A(p,A,nav,ai,1)
				nav=nav+1
				action=.true.
				if(nav.ge.nv) final=.true.
			end if
			if(any(abs(ba_crct(1:nav)).le.zero).and.mthd.eq.1 .and. ai.lt.0) then
				ai=abs(ai)
				ruv(A(ai))=g
				call shift_A(p,A,nav,ai,-1)
				nav=nav-1
				if(abs(ba_crct(ai)).eq.0.)	df(nstp)=nav+1
				action=.true.
				final=.false.
			end if
		else
			action=.true.
			j=p-nav
			i=nav
			if(test(1)) then
				do ai=1,j
					if(dabsruac(ai).ge.eps) then
						call shift_A(p,A,nav,ai+i-nav,1)
						nav=nav+1
						if(nav.ge.nv) then
							final=.true.
							exit
						end if
					end if
				end do
			end if
			if(test(2)) then
				do ai=1,i
					if(ba_crct(ai)*ruv(A(ai)).lt.0) then
						call shift_A(p,A,nav,ai,-1)
						nav=nav-1
						final=.false.
					end if
				end do
			end if
		end if
		if(action) then
			deallocate(ba,dba,ba_crct,dabsruac,stat=conv)
			if(conv.ne.0) then
				conv=4
				np=nstp
				return
			end if
		end if
		if(test(3)) exit
	end do
	if(k.ge.np) then
		conv=3
		np=nstp
		return
	end if
	np=nstp
	deallocate(eta,mu,ruv,dmu_de,sqrt_ib,X2,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine dglars_pc_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the prediction step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine prd_b(mthd,g,g0,n,p,X,X2,A,nav,ba,mu,dmu_de,sqrt_ib,ruv,dg_max,dba,dg,conv,ai,final)
	integer	::	mthd,n,p,nav,conv,A(p),ai,j
	double precision	::	g,g0,X(n,p),X2(n,p),ba(0:nav),mu(n),dmu_de(n),sqrt_ib(p),ruv(p),dg_max,dg,dg_out,dba(0:nav)
	logical	:: final
	double precision,	dimension(:),	allocatable	::	d2mu_de2
	double precision,	dimension(:,:), allocatable	::	Drua
	allocate(d2mu_de2(1:n),Drua(0:nav,0:nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	dba(0)=0.
	dba(1:nav)=sign(1.,real(ruv(A(1:nav))))
	call d2mu_de2_mk_b(n,mu,dmu_de,d2mu_de2)	
	call jacob(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),mu,dmu_de,d2mu_de2,sqrt_ib(A(1:nav)),ruv(A(1:nav)),Drua)
	call solve(nav+1,-Drua,dba,conv)
	if(conv.ne.0) then
		conv=2
		return
	end if
	if(final) then
		if(dg_max.gt.0.) then
			dg=min(dg_max,g-g0)
		else
			dg=g-g0
		end if
	else
		call step_size(n,g,g0,p,nav,X(:,A(1:nav)),X(:,A((nav+1):p)),X2(:,A((nav+1):p)),dba, &
			dmu_de,d2mu_de2,sqrt_ib(A((nav+1):p)),ruv(A((nav+1):p)),dg_max,ai,dg,conv)
		if(conv.ne.0) return
	end if
	if(mthd.eq.1) then
		do j=1,nav
			dg_out=ba(j)/dba(j)
			if(dg_out.gt.0..and.dg_out.le.dg) then
				dg=dg_out
				ai=-j
			end if
		end do
	end if
	deallocate(d2mu_de2,Drua,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine prd_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine prd_p(mthd,g,g0,n,p,X,X2,A,nav,ba,mu,dmu_de,sqrt_ib,ruv,dg_max,dba,dg,conv,ai,final)
	integer	::	mthd,n,p,nav,conv,A(p),ai,j
	double precision	::	g,g0,X(n,p),X2(n,p),ba(0:nav),mu(n),dmu_de(n),sqrt_ib(p),ruv(p),dg_max,dg,dg_out,dba(0:nav)
	logical	:: final
	double precision,	dimension(:),	allocatable	::	d2mu_de2
	double precision,	dimension(:,:), allocatable	::	Drua
	allocate(d2mu_de2(1:n),Drua(0:nav,0:nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	dba(0)=0.
	dba(1:nav)=sign(1.,real(ruv(A(1:nav))))
	call d2mu_de2_mk_p(n,mu,d2mu_de2)
	call jacob(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),mu,dmu_de,d2mu_de2,sqrt_ib(A(1:nav)),ruv(A(1:nav)),Drua)
	call solve(nav+1,-Drua,dba,conv)
	if(conv.ne.0)	then
		conv=2
		return
	end if
	if(final) then
		if(dg_max.gt.0.) then
			dg=min(dg_max,g-g0)
		else
			dg=g-g0
		end if
	else
		call step_size(n,g,g0,p,nav,X(:,A(1:nav)),X(:,A((nav+1):p)),X2(:,A((nav+1):p)),dba, &
			dmu_de,d2mu_de2,sqrt_ib(A((nav+1):p)),ruv(A((nav+1):p)),dg_max,ai,dg,conv)
		if(conv.ne.0) return
	end if
	if(mthd.eq.1) then
		do j=1,nav
			dg_out=ba(j)/dba(j)
			if(dg_out.gt.0..and.dg_out.le.dg) then
				dg=dg_out
				ai=-j
			end if
		end do
	end if
	deallocate(d2mu_de2,Drua,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine prd_p
! generic subroutine
subroutine step_size(n,g,g0,p,nav,Xa,Xac,X2ac,dba,dmu_de,d2mu_de2,sqrt_i_bac,ruac,dg_max,ai,dg,conv)
	integer	::	n,p,nav,conv,k,h,ai
	double precision	::	g,g0,dg,Xa(n,nav),Xac(n,p-nav),dmu_de(n),d2mu_de2(n),sqrt_i_bac(p-nav)
	double precision	::	X2ac(n,p-nav),ruac(p-nav),dg_max,druac,druac_tmp,dg_opt,dba(0:nav)
	double precision,	dimension(:),	allocatable	::	i_bac
	allocate(i_bac(1:(p-nav)),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	i_bac=sqrt_i_bac**2
	dg=g
	do	k=1,p-nav
		if(abs(ruac(k)).ne.g) then
			druac_tmp=0.5*ruac(k)/i_bac(k)
			druac=-dba(0)*(dot_product(Xac(:,k),dmu_de)/sqrt_i_bac(k) &
				+druac_tmp*dot_product(X2ac(:,k),d2mu_de2))
			do h=1,nav
				druac=druac-dba(h)*(dot_product(Xac(:,k),dmu_de*Xa(:,h))/sqrt_i_bac(k) &
					+druac_tmp*dot_product(X2ac(:,k),d2mu_de2*Xa(:,h)))
			end do
			dg_opt=(g-ruac(k))/(1.-druac)
			if(dg_opt<0.or.dg_opt>g) then
				dg_opt=(g+ruac(k))/(1.+druac)
			end if
			if(dg_opt<dg.and.dg_opt.gt.0.) then
				dg=dg_opt
				ai=k
			end if
		end if
	end do
	if(dg_max>0. .and. dg>dg_max) then
		dg=dg_max
		ai=0
	end if
	if(dg>g-g0) then
		dg=g-g0
		ai=0
	end if
	deallocate(i_bac,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if	
end subroutine step_size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the corrector step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine crct_b(n,nav,Xa,X2a,y,ba,dba,g,dg,rua,NReps,nNR,mu,dmu_de,ba_crct,conv)
	integer	::	n,nav,nNR,conv
	double precision	::	Xa(n,nav),X2a(n,nav),y(n),g,dg,rua(nav),NReps,mu(n),dmu_de(n),va(nav)
	double precision,	dimension(0:nav)	:: ba,dba,ba_crct,ba_prd
	va=(g-dg)*sign(1.,real(rua))
	ba_prd=ba-dg*dba
	call newt_b(n,nav,va,Xa,X2a,y,NReps,nNR,mu,dmu_de,ba_prd,conv)
	if(conv.ne.0) return
	ba_crct=ba_prd
end subroutine crct_b
subroutine newt_b(n,nav,va,Xa,X2a,y,NReps,n_step,mu,dmu_de,ba_crct,conv)
	integer	:: n,nav,n_step,conv,i
	double precision	:: ba_crct(0:nav),va(nav),Xa(n,nav),X2a(n,nav),y(n),NReps,mu(n),dmu_de(n),sum_abs_f,sum_abs_db
	double precision,	dimension(:),	allocatable	::	dba,eta,d2mu_de2,sqrt_i_ba,rua
	double precision,	dimension(:,:),	allocatable	::	Drua
	allocate(dba(0:nav),eta(n),d2mu_de2(n),sqrt_i_ba(nav),rua(nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	allocate(Drua(0:nav,0:nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	!
	do	i=1, n_step
		call eta_mk(n,nav,Xa,ba_crct,eta)
		call mu_mk_b(n,eta,mu)
		call dmu_de_mk_b(n,mu,dmu_de)
		call sqrt_i_b_mk(n,nav,X2a,dmu_de,sqrt_i_ba)
		call rao_score(n,nav,Xa,y,mu,sqrt_i_ba,rua)
		dba(0)=sum(y-mu)
		dba(1:nav)=rua-va
		sum_abs_f=sum(abs(dba))
		if(sum_abs_f.le.NReps) exit
		call d2mu_de2_mk_b(n,mu,dmu_de,d2mu_de2)
		call jacob(n,nav,Xa,X2a,mu,dmu_de,d2mu_de2,sqrt_i_ba,rua,Drua)
		call solve(nav+1,Drua,dba,conv)
		if(conv.ne.0) then
			conv=2
			return
		end if
		sum_abs_db=sum(abs(dba))
		if(sum_abs_db.ne.sum_abs_db) then
			conv=2
			return
		end if
		ba_crct=ba_crct+dba
	end do
	!
	if(i.eq.n_step) then
		conv=2
		return
	end if
	!
	deallocate(dba,eta,d2mu_de2,sqrt_i_ba,rua,Drua,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine newt_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine crct_p(n,nav,Xa,X2a,y,ba,dba,g,dg,rua,NReps,nNR,mu,dmu_de,ba_crct,conv)
	integer	::	n,nav,nNR,conv
	double precision	::	Xa(n,nav),X2a(n,nav),y(n),g,dg,rua(nav),NReps,mu(n),dmu_de(n),va(nav)
	double precision,	dimension(0:nav)	:: ba,dba,ba_crct,ba_prd
	va=(g-dg)*sign(1.,real(rua))
	ba_prd=ba-dg*dba
	call newt_p(n,nav,va,Xa,X2a,y,NReps,nNR,mu,dmu_de,ba_prd,conv)
	if(conv.ne.0) return
	ba_crct=ba_prd
end subroutine crct_p
subroutine newt_p(n,nav,va,Xa,X2a,y,NReps,n_step,mu,dmu_de,ba_crct,conv)
	integer	:: n,nav,n_step,conv,i
	double precision	:: ba_crct(0:nav),va(nav),Xa(n,nav),X2a(n,nav),y(n),NReps,mu(n),dmu_de(n),sum_abs_f,sum_abs_db
	double precision,	dimension(:),	allocatable	::	dba,eta,d2mu_de2,sqrt_i_ba,rua
	double precision,	dimension(:,:),	allocatable	::	Drua
	allocate(dba(0:nav),eta(n),d2mu_de2(n),sqrt_i_ba(nav),rua(nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	allocate(Drua(0:nav,0:nav),stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
	!
	do	i=1, n_step
		call eta_mk(n,nav,Xa,ba_crct,eta)
		call mu_mk_p(n,eta,mu)
		call dmu_de_mk_p(n,mu,dmu_de)
		call sqrt_i_b_mk(n,nav,X2a,dmu_de,sqrt_i_ba)
		call rao_score(n,nav,Xa,y,mu,sqrt_i_ba,rua)
		dba(0)=sum(y-mu)
		dba(1:nav)=rua-va
		sum_abs_f=sum(abs(dba))
		if(sum_abs_f.le.NReps) exit
		call d2mu_de2_mk_p(n,mu,d2mu_de2)
		call jacob(n,nav,Xa,X2a,mu,dmu_de,d2mu_de2,sqrt_i_ba,rua,Drua)
		call solve(nav+1,Drua,dba,conv)
		if(conv.ne.0) then
			conv=2
			return
		end if
		sum_abs_db=sum(abs(dba))
		if(sum_abs_db.ne.sum_abs_db) then
			conv=2
			return
		end if
		ba_crct=ba_crct+dba
	end do
	!
	if(i.eq.n_step) then
		conv=2
		return
	end if
	!
	deallocate(dba,eta,d2mu_de2,sqrt_i_ba,rua,Drua,stat=conv)
	if(conv.ne.0) then
		conv=4
		return
	end if
end subroutine newt_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used for the cross-validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine cvdglars_pc_b(n,p,X,y,foldid,nfold,ng,g,b,dev_m,dev_v,g_hat,nv,mthd,g0,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	integer	:: n,p,foldid(n),nfold,ng,nv,mthd,np,ncrct,nNR,conv
	double precision	:: X(n,p),y(n),g(ng),b(0:p),dev_m(ng),dev_v(ng),g_hat,g0,dg_max,eps,cf,NReps
	! internal variables
	integer	:: i,j,lfold,np_cv,A(p),nav,g_id,df(np)
	double precision	:: b_cv(0:p,np),ru_cv(p,np),dev(np),g_seq(np),dev_cv(ng,nfold)
	lfold=n/nfold
	do i=1,nfold
		b_cv=0.
		ru_cv=0.
		dev=0.
		A=(/ (j,j=1,p) /)
		nav=0
		df=0
		g_seq=0.
		np_cv=np
		call dglars_pc_b(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),b_cv,ru_cv,dev,A,nv,nav,df,g_seq,mthd,g0,&
		dble(2.),dg_max,eps,np_cv,ncrct,cf,NReps,nNR,conv)
		if(conv.ne.0) return
		call predict_b(lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),np_cv,b_cv,g_seq,ng,g,dev_cv(:,i))
		foldid=cshift(foldid,lfold)
	end do
	dev_m=sum(dev_cv,dim=2)/nfold
	dev_v=sum(dev_cv**2,dim=2)/nfold-dev_m**2
	dev_v=dev_v*nfold/(nfold-1)
	b_cv=0.
	ru_cv=0.
	dev=0.
	A=(/ (j,j=1,p) /)
	nav=0
	df=0
	g_seq=0.
	g_id=minloc(dev_m,dim=1)
	g_hat=g(g_id)
	call dglars_pc_b(n,p,X,y,b_cv,ru_cv,dev,A,nv,nav,df,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	if(conv.ne.0) return
	b=b_cv(:,np)
	g_hat=g0
end subroutine cvdglars_pc_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine cvdglars_pc_p(n,p,X,y,foldid,nfold,ng,g,b,dev_m,dev_v,g_hat,nv,mthd,g0,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	integer	:: n,p,foldid(n),nfold,ng,nv,mthd,np,ncrct,nNR,conv
	double precision	:: X(n,p),y(n),g(ng),b(0:p),dev_m(ng),dev_v(ng),g_hat,g0,dg_max,eps,cf,NReps
	! internal variables
	integer	:: i,j,lfold,np_cv,A(p),nav,g_id,df(np)
	double precision	:: b_cv(0:p,np),ru_cv(p,np),dev(np),g_seq(np),dev_cv(ng,nfold)
	lfold=n/nfold
	do i=1,nfold
		b_cv=0.
		ru_cv=0.
		dev=0.
		A=(/ (j,j=1,p) /)
		nav=0
		df=0
		g_seq=0.
		np_cv=np
		call dglars_pc_p(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),b_cv,ru_cv,dev,A,nv,nav,df,g_seq,mthd,g0,&
		dble(2.),dg_max,eps,np_cv,ncrct,cf,NReps,nNR,conv)
		if(conv.ne.0) return
		call predict_p(lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),np_cv,b_cv,g_seq,ng,g,dev_cv(:,i))
		foldid=cshift(foldid,lfold)
	end do
	dev_m=sum(dev_cv,dim=2)/nfold
	dev_v=sum(dev_cv**2,dim=2)/nfold-dev_m**2
	dev_v=dev_v*nfold/(nfold-1)
	b_cv=0.
	ru_cv=0.
	dev=0.
	A=(/ (j,j=1,p) /)
	nav=0
	df=0
	g_seq=0.
	g_id=minloc(dev_m,dim=1)
	g_hat=g(g_id)
	call dglars_pc_p(n,p,X,y,b_cv,ru_cv,dev,A,nv,nav,df,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
	if(conv.ne.0) return
	b=b_cv(:,np)
	g_hat=g0
end subroutine cvdglars_pc_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generic subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eta_mk(n,nav,Xa,ba,eta)
	integer	::	n,nav,j
	double precision	::	Xa(n,nav),eta(n)
	double precision, dimension(0:nav)	:: ba
	eta=ba(0)
	do j=1,nav
		eta=eta+Xa(:,j)*ba(j)
	end do
end subroutine eta_mk
subroutine sqrt_i_b_mk(n,p,X2,dmu_de,sqrt_ib)
	integer	::	n,p,j
	double precision	::	X2(n,p),dmu_de(n),sqrt_ib(p)
	do j=1,p
		sqrt_ib(j)=sqrt(sum(dmu_de*X2(:,j)))
	end do
end subroutine sqrt_i_b_mk
subroutine rao_score(n,p,X,y,mu,sqrt_ib,ruv)
	integer	::	n,p,j
	double precision	::	X(n,p),y(n),mu(n),sqrt_ib(p),ruv(p),r(n)
	r=y-mu
	do j=1,p
		ruv(j)=dot_product(r,X(:,j))/sqrt_ib(j)
	end do
end subroutine rao_score
subroutine shift_A(p,A,nav,ai,action)
	integer	:: nav,ai,p,tmp,A(p),action
	if(action.eq.1) then
		tmp = A(nav+1)
		A(nav+1) = A(nav+ai)
		A(nav+ai) = tmp
	end if
	if(action.eq.-1) then
		tmp=A(ai)
		A(ai)=A(nav)
		A(nav)=tmp
	end if
end subroutine shift_A
subroutine	solve(nba,Drua,dba,conv)
	integer	::	nba,ipiv(nba),conv
	integer,	parameter	::	p=1
	double precision	:: Drua(nba,nba),dba(nba,p)
	call dgesv(nba,p,Drua,nba,ipiv,dba,nba,conv)
	if(conv.ne.0) then
		conv=1
		return
	end if
end subroutine
subroutine jacob(n,nav,Xa,X2a,mu,dmu_de,d2mu_de2,sqrt_i_ba,rua,Drua)
	integer	::	n,nav,h,k
	double precision	::	Xa(n,nav),X2a(n,nav),mu(n),dmu_de(n),d2mu_de2(n),sqrt_i_ba(nav),rua(nav),Drua(0:nav,0:nav),xkxh,i_ba(1:nav)
	i_ba=sqrt_i_ba**2
	Drua(0,0)=sum(dmu_de)
	do h=1,nav
		Drua(0,h)= sum(dmu_de*Xa(:,h))
		Drua(h,0)= Drua(0,h)/sqrt_i_ba(h) +0.5*rua(h)/i_ba(h)*dot_product(X2a(:,h),d2mu_de2)
	end do
	if (nav>1) then
		do k=1,nav-1
			Drua(k,k)=dot_product(Xa(:,k),dmu_de*Xa(:,k))/sqrt_i_ba(k) &
			+0.5*rua(k)/i_ba(k)*dot_product(X2a(:,k),d2mu_de2*Xa(:,k))
			do h=k+1,nav
				xkxh=dot_product(Xa(:,k),dmu_de*Xa(:,h))
				Drua(k,h)=xkxh/sqrt_i_ba(k) &
				+0.5*rua(k)/i_ba(k)*dot_product(X2a(:,k),d2mu_de2*Xa(:,h))
				Drua(h,k)=xkxh/sqrt_i_ba(h) &
				+0.5*rua(h)/i_ba(h)*dot_product(X2a(:,h),d2mu_de2*Xa(:,k))
			end do
		end do
	end if
	Drua(nav,nav)=dot_product(Xa(:,nav),dmu_de*Xa(:,nav))/sqrt_i_ba(nav) &
			+0.5*rua(nav)/i_ba(nav)*dot_product(X2a(:,nav),d2mu_de2*Xa(:,nav))
end subroutine jacob
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines related to the exponential family
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY!
!!!!!!!!!!!!!!!!!!
subroutine mu_mk_b(n,eta,mu)
	integer	::	n
	double precision	::	eta(n),mu(n)
	mu = 1/(1+exp(-eta))
end subroutine mu_mk_b
subroutine dmu_de_mk_b(n,mu,dmu_de)
	integer	::	n
	double precision	::	mu(n),dmu_de(n)
	dmu_de=mu*(1-mu)
end subroutine dmu_de_mk_b
subroutine	d2mu_de2_mk_b(n,mu,dmu_de,d2mu_de2)
	integer	::	n
	double precision	::	mu(n),dmu_de(n),d2mu_de2(n)
	d2mu_de2=dmu_de*(1-2*mu)
end subroutine	d2mu_de2_mk_b
subroutine deviance_b(n,y,mu,dev)
	integer	::	n
	double precision	:: y(n),mu(n),dev
	dev=-2*sum(log(mu),mask=y>0.5)-2*sum(log(1-mu),mask=y<0.5)
end subroutine deviance_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine mu_mk_p(n,eta,mu)
	integer	::	n
	double precision	::	eta(n),mu(n)
	mu = exp(eta)
end subroutine mu_mk_p
subroutine dmu_de_mk_p(n,mu,dmu_de)
	integer	::	n
	double precision	::	mu(n),dmu_de(n)
	dmu_de=mu
end subroutine dmu_de_mk_p
subroutine	d2mu_de2_mk_p(n,mu,d2mu_de2)
	integer	::	n
	double precision	::	mu(n),d2mu_de2(n)
	d2mu_de2=mu
end subroutine	d2mu_de2_mk_p
subroutine deviance_p(n,y,mu,dev)
	integer	::	n
	double precision	:: y(n),mu(n),dev
	dev=2*(sum(y*log(y/mu),mask=y>0.)-sum(y-mu))
end subroutine deviance_p
