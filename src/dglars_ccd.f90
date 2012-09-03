subroutine dglars_ccd_b(n,p,X,y,np,g0,g_hat,nstp,eps,mthd,b,dev,g_seq,A,df,nav,conv)
	double precision,	parameter :: big=9.9e30
	integer	:: n,p,np,nstp,mthd,A(p),df(np),nav,conv
	double precision	:: X(n,p),y(n),g0,g_hat,eps,b(0:p,np),dev(np),g_seq(np)
	! internal variables
	integer	:: step,m,ax(p),ax_old(p),A_old(p),nav_old,step_ccd,ai
	double precision	:: ob(0:p),nb(0:p),eta(n),mu(n),r(n),V(n),sV,g,Xai(n),dlai,Iai,I(p),cf,db_max,zai,adlai,dbai,adlai_max
	logical	:: new_step
	step_ccd=0
	ax=0
	nav=0
	A=0
	ob=0.
	nb=0.
	mu=sum(y)/n
	r=y-mu
	V=mu*(1.-mu)
	sV=n*V(1)
	ob(0)=log(mu(1)/(1.-mu(1)))	
	if(np.gt.1) then
		g=0.
		step=1
		do m=1,p
			Xai=X(:,m)
			dlai=dot_product(Xai,r)
			Iai=dot_product(V,Xai**2)
			I(m)=Iai
			g=max(g,dlai/sqrt(Iai))
		end do
		if(g_hat.ne.2.) g0=g0+g_hat*(g-g0)
		cf=exp((log(g0)-log(g))/(np-1))
		b(0,step)=ob(0)
		dev(step)=-2.*(sum(y)*log(mu(1))+((n-sum(y))*log(1.-mu(1))))
		g_seq(step)=g
		df(step)=1
		else
			cf=1.
			step=0
	end if
	do
		new_step=.true.
		g=g*cf
		do m=1,p
			if(ax(m).eq.0) then
				Xai=X(:,m)
				dlai=dot_product(Xai,r)
				Iai=dot_product(V,Xai**2)
				I(m)=Iai
				if(step.gt.1.and.abs(dlai).ge.(sqrt(Iai)*g_seq(step))) then
					new_step=.false.
					ax_old(m)=1
					nav_old=nav_old+1
					A_old(nav_old)=m
				end if
				if(new_step.and.abs(dlai).ge.(sqrt(Iai)*g)) then
					ax(m)=1
					nav=nav+1
					A(nav)=m
				end if
			end if
		end do
		if(new_step) then
			step=step+1
			else
				g=g/cf
				ax=ax_old
				nav=nav_old
				A=A_old
		end if
		do
			step_ccd=step_ccd+1
			if(step_ccd.gt.nstp) then
				np=step-1
				conv=3
				return
			end if
			db_max=0.
			do m=1,nav
				ai=A(m)
				Xai=X(:,ai)
				Iai=I(ai)
				dlai=dot_product(Xai,r)
				zai=dlai+Iai*ob(ai)
				adlai=abs(zai)-sqrt(Iai)*g
				nb(ai)=sign(1.,real(zai))*adlai/Iai			
				if(dlai*nb(ai).le.0..and.mthd.eq.1) then
					nb(ai)=0.
					ax(ai)=0
					A(m)=A(nav)
					nav=nav-1
					db_max=big
					exit
				end if
				dbai=nb(ai)-ob(ai)
				db_max=max(db_max,abs(dbai))
				r=r-dbai*V*Xai
				ob(ai)=nb(ai)
			end do
			dbai=sum(r)/sV
			nb(0)=ob(0)+dbai
			r=r-dbai*V
			ob(0)=nb(0)
			db_max=max(db_max,abs(dbai))
			if(db_max.le.eps) then
				eta=nb(0)
				do m=1,nav
					ai=A(m)
					eta=eta+X(:,ai)*nb(ai)
				end do
				mu=1./(1.+exp(-eta))
				r=y-mu
				V=mu*(1.-mu)
				sV=sum(V)
				adlai_max=0.
				adlai_max=max(adlai_max,abs(sum(r)))
				do m=1,nav
					ai=A(m)
					Xai=X(:,ai)
					dlai=dot_product(Xai,r)
					Iai=dot_product(V,Xai**2)
					I(ai)=Iai
					adlai_max=max(adlai_max,abs(abs(dlai)-sqrt(Iai)*g))
				end do
				if(adlai_max.le.eps) exit
			end if
		end do
		b(0,step)=nb(0)
		b(A(1:nav),step)=nb(A(1:nav))
		df(step)=nav+1
		dev(step)=-2*sum(log(mu),mask=y>0.5)-2*sum(log(1-mu),mask=y<0.5)
		g_seq(step)=g
		ax_old=ax
		nav_old=nav
		A_old=A
		if(step.eq.np) exit
	end do
end subroutine dglars_ccd_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine dglars_ccd_p(n,p,X,y,np,g0,g_hat,nstp,eps,mthd,b,dev,g_seq,A,df,nav,conv)
	double precision,	parameter :: big=9.9e30
	integer	:: n,p,np,nstp,mthd,A(p),df(np),nav,conv
	double precision	:: X(n,p),y(n),g0,g_hat,eps,b(0:p,np),dev(np),g_seq(np)
	! internal variables
	integer	:: step,m,ax(p),ax_old(p),A_old(p),nav_old,step_ccd,ai
	double precision	:: ob(0:p),nb(0:p),eta(n),mu(n),r(n),V(n),sV,g,Xai(n),dlai,Iai,I(p),cf,db_max,zai,adlai,dbai,adlai_max
	logical	:: new_step
	step_ccd=0
	ax=0
	nav=0
	A=0
	ob=0.
	nb=0.
	mu=sum(y)/n
	r=y-mu
	V=mu
	sV=n*V(1)
	ob(0)=log(mu(1))
	if(np.gt.1) then
		g=0.
		step=1
		do m=1,p
			Xai=X(:,m)
			dlai=dot_product(Xai,r)
			Iai=dot_product(V,Xai**2)
			I(m)=Iai
			g=max(g,dlai/sqrt(Iai))
		end do
		if(g_hat.ne.2.) g0=g0+g_hat*(g-g0)
		cf=exp((log(g0)-log(g))/(np-1))
		b(0,step)=ob(0)
		dev(step)=2*sum(y*log(y/mu),mask=y.gt.0.5)
		g_seq(step)=g
		df(step)=1
		else
			cf=1.
			step=0
	end if
	do
		new_step=.true.
		g=g*cf
		do m=1,p
			if(ax(m).eq.0) then
				Xai=X(:,m)
				dlai=dot_product(Xai,r)
				Iai=dot_product(V,Xai**2)
				I(m)=Iai
				if(step.gt.1.and.abs(dlai).ge.(sqrt(Iai)*g_seq(step))) then
					new_step=.false.
					ax_old(m)=1
					nav_old=nav_old+1
					A_old(nav_old)=m
				end if
				if(new_step.and.abs(dlai).ge.(sqrt(Iai)*g)) then
					ax(m)=1
					nav=nav+1
					A(nav)=m
				end if
			end if
		end do
		if(new_step) then
			step=step+1
			else
				g=g/cf
				ax=ax_old
				nav=nav_old
				A=A_old
		end if
		do
			step_ccd=step_ccd+1
			if(step_ccd.gt.nstp) then
				np=step-1
				conv=3
				return
			end if
			db_max=0.
			do m=1,nav
				ai=A(m)
				Xai=X(:,ai)
				Iai=I(ai)
				dlai=dot_product(Xai,r)
				zai=dlai+Iai*ob(ai)
				adlai=abs(zai)-sqrt(Iai)*g
				nb(ai)=sign(1.,real(zai))*adlai/Iai
				if(dlai*nb(ai).le.0..and.mthd.eq.1) then
					nb(ai)=0.
					ax(ai)=0
					A(m)=A(nav)
					nav=nav-1
					db_max=big
					exit
				end if
				dbai=nb(ai)-ob(ai)
				db_max=max(db_max,abs(dbai))
				r=r-dbai*V*Xai
				ob(ai)=nb(ai)
			end do
			dbai=sum(r)/sV
			nb(0)=ob(0)+dbai
			r=r-dbai*V
			ob(0)=nb(0)
			db_max=max(db_max,abs(dbai))
			if(db_max.le.eps) then
				eta=nb(0)
				do m=1,nav
					ai=A(m)
					eta=eta+X(:,ai)*nb(ai)
				end do
				mu=exp(eta)
				r=y-mu
				V=mu
				sV=sum(V)
				adlai_max=0.
				adlai_max=max(adlai_max,abs(sum(r)))
				do m=1,nav
					ai=A(m)
					Xai=X(:,ai)
					dlai=dot_product(Xai,r)
					Iai=dot_product(V,Xai**2)
					I(ai)=Iai
					adlai_max=max(adlai_max,abs(abs(dlai)-sqrt(Iai)*g))
				end do
				if(adlai_max.le.eps) exit
			end if
		end do
		b(0,step)=nb(0)
		b(A(1:nav),step)=nb(A(1:nav))
		df(step)=nav+1
		dev(step)=2*(sum(y*log(y/mu),mask=y.gt.0.5)-sum(r))
		g_seq(step)=g
		ax_old=ax
		nav_old=nav
		A_old=A
		if(step.eq.np) exit
	end do
end subroutine dglars_ccd_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used for the cross-validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine cvdglars_ccd_b(n,p,X,y,foldid,nfold,ng,g,b,dev_m,dev_v,g_hat,mthd,g0,eps,np,nstp,conv)
	integer	:: n,p,foldid(n),nfold,ng,mthd,np,nstp,conv
	double precision	:: X(n,p),y(n),g(ng),b(0:p),dev_m(ng),dev_v(ng),g_hat,g0,eps
	! internal variables
	integer	:: i,j,lfold,A(p),nav,g_id,df(np)
	double precision	:: b_cv(0:p,np),dev(np),g_seq(np),dev_cv(ng,nfold)
	lfold=n/nfold
	do i=1,nfold
		b_cv=0.
		dev=0.
		A=(/ (j,j=1,p) /)
		nav=0
		df=0
		g_seq=0.
		call dglars_ccd_b(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),np,g0,dble(2.),nstp,eps,mthd,b_cv,dev,g_seq,A,df,nav,conv)
		if(conv.ne.0) return
		call predict_b(lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),np,b_cv,g_seq,ng,g,dev_cv(:,i))
		foldid=cshift(foldid,lfold)
	end do
	dev_m=sum(dev_cv,dim=2)/nfold
	dev_v=sum(dev_cv**2,dim=2)/nfold-dev_m**2
	dev_v=dev_v*nfold/(nfold-1)
	b_cv=0.
	dev=0.
	A=(/ (j,j=1,p) /)
	nav=0
	df=0
	g_seq=0.
	g_id=minloc(dev_m,dim=1)
	g_hat=g(g_id)
	call dglars_ccd_b(n,p,X,y,np,g0,g_hat,nstp,eps,mthd,b_cv,dev,g_seq,A,df,nav,conv)
	if(conv.ne.0) return
	b=b_cv(:,np)
	g_hat=g0
end subroutine cvdglars_ccd_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine cvdglars_ccd_p(n,p,X,y,foldid,nfold,ng,g,b,dev_m,dev_v,g_hat,mthd,g0,eps,np,nstp,conv)
	integer	:: n,p,foldid(n),nfold,ng,mthd,np,nstp,conv
	double precision	:: X(n,p),y(n),g(ng),b(0:p),dev_m(ng),dev_v(ng),g_hat,g0,eps
	! internal variables
	integer	:: i,j,lfold,A(p),nav,g_id,df(np)
	double precision	:: b_cv(0:p,np),dev(np),g_seq(np),dev_cv(ng,nfold)
	lfold=n/nfold
	do i=1,nfold
		b_cv=0.
		dev=0.
		A=(/ (j,j=1,p) /)
		nav=0
		df=0
		g_seq=0.
		call dglars_ccd_p(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),np,g0,dble(2.),nstp,eps,mthd,b_cv,dev,g_seq,A,df,nav,conv)
		if(conv.ne.0) return
		call predict_p(lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),np,b_cv,g_seq,ng,g,dev_cv(:,i))
		foldid=cshift(foldid,lfold)
	end do
	dev_m=sum(dev_cv,dim=2)/nfold
	dev_v=sum(dev_cv**2,dim=2)/nfold-dev_m**2
	dev_v=dev_v*nfold/(nfold-1)
	b_cv=0.
	dev=0.
	A=(/ (j,j=1,p) /)
	nav=0
	df=0
	g_seq=0.
	g_id=minloc(dev_m,dim=1)
	g_hat=g(g_id)
	call dglars_ccd_p(n,p,X,y,np,g0,g_hat,nstp,eps,mthd,b_cv,dev,g_seq,A,df,nav,conv)
	if(conv.ne.0) return
	b=b_cv(:,np)
	g_hat=g0
end subroutine cvdglars_ccd_p
