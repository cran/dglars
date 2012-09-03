!!!!!!!!!!!!!!!!!!!
! BINOMIAL FAMILY !
!!!!!!!!!!!!!!!!!!!
subroutine predict_b(n,p,X,y,np,b,g_seq,ng,g,dev)
	integer	:: n,p,np,ng,conv
	double precision	:: X(n,p),y(n),b(0:p,np),g_seq(np),g(ng),dev(ng)
	! internal variables
	integer	:: i,m,right,left
	double precision	:: g_min,g_max,g_frac(np),b_prd(0:p),dlt,eta(n),mu(n)
	g_min=minval(g_seq)
	g_max=maxval(g_seq)
	g_frac=(g_seq-g_min)/(g_max-g_min)
	do i=1,ng
		if(any(g_frac.eq.g(i))) then		
			right=count(g_frac.ge.g(i))
			b_prd=b(:,right)
		else
			right=count(g_frac.gt.g(i))
			left=right+1
			dlt=(g(i)-g_frac(right))/(g_frac(left)-g_frac(right))
			b_prd=b(:,right)+dlt*(b(:,left)-b(:,right))
		end if
		eta=b_prd(0)
		do m=1,p
			if(abs(b_prd(m)).gt.0.) eta=eta+X(:,m)*b_prd(m)
		end do
		call mu_mk_b(n,eta,mu)
		call deviance_b(n,y,mu,dev(i))
	end do
end subroutine predict_b
!!!!!!!!!!!!!!!!!!
! POISSON FAMILY !
!!!!!!!!!!!!!!!!!!
subroutine predict_p(n,p,X,y,np,b,g_seq,ng,g,dev)
	integer	:: n,p,np,ng,conv
	double precision	:: X(n,p),y(n),b(0:p,np),g_seq(np),g(ng),dev(ng)
	! internal variables
	integer	:: i,m,right,left
	double precision	:: g_min,g_max,g_frac(np),b_prd(0:p),dlt,eta(n),mu(n)
	g_min=minval(g_seq)
	g_max=maxval(g_seq)
	g_frac=(g_seq-g_min)/(g_max-g_min)
	do i=1,ng
		if(any(g_frac.eq.g(i))) then		
			right=count(g_frac.ge.g(i))
			b_prd=b(:,right)
		else
			right=count(g_frac.gt.g(i))
			left=right+1
			dlt=(g(i)-g_frac(right))/(g_frac(left)-g_frac(right))
			b_prd=b(:,right)+dlt*(b(:,left)-b(:,right))
		end if
		eta=b_prd(0)
		do m=1,p
			if(abs(b_prd(m)).gt.0.) eta=eta+X(:,m)*b_prd(m)
		end do
		call mu_mk_p(n,eta,mu)
		call deviance_p(n,y,mu,dev(i))
	end do
end subroutine predict_p
