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