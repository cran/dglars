!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eta_mk(n, nav, Xa, ba, eta)
integer :: n,nav,j
double precision :: Xa(n,nav),eta(n)
double precision, dimension(0:nav) :: ba
eta = ba(0)
do j = 1, nav
    eta = eta + Xa(:,j) * ba(j)
end do
end subroutine eta_mk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutines mu_mk, dmu_de_mk and d2mu_de2_mk are used to compute the
! inverse of the link fucntion, its first and second derivative, respectively.
! Arguments:
! linkid: integer used to specify the link function:
!       linkid = 1 => identity
!       linkid = 2 => log
!       linkid = 3 => inverse
!       linkid = 4 => sqrt
!       linkid = 5 => cloglog
!       linkid = 6 => probit
!       linkid = 7 => cauchit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linkfun(linkid, mu, eta)
integer :: linkid
double precision :: mu,eta
double precision, external :: qnorm,qcauchy
select case (linkid)
    case (1)
        eta=mu
    case (2)
        eta=log(mu)
    case (3)
        eta=1.d0/mu
    case (4)
        eta=sqrt(mu)
    case (5)
        eta=log(-log(1.d0-mu))
    case (6)
        eta=qnorm(mu)
    case (7)
        eta=qcauchy(mu)
end select
end subroutine linkfun
!
subroutine mu_mk(linkid,n,eta,mi,mu)
integer :: linkid,n,i
double precision :: eta(n),mi(n),mu(n),thresh,etai
double precision, parameter :: eps=2.d0**(-52)
double precision, external :: qnorm,qcauchy,pcauchy
select case (linkid)
    case (1)
        mu=eta
    case (2)
        forall(i=1:n) mu(i)=mi(i)*max(exp(eta(i)),eps)
    case (3)
        mu=1.d0/eta
    case (4)
        mu=eta**2
    case (5)
        forall(i=1:n) mu(i)=mi(i)*max(min(1.d0-exp(-exp(eta(i))),1.d0-eps),eps)
    case (6)
        do i=1,n
            thresh=-qnorm(eps)
            etai=min(max(eta(i),-thresh),thresh)
            mu(i)=mi(i)/2.d0*erfc(-etai/sqrt(2.d0))
        end do
    case (7)
        do i=1,n
            thresh=-qcauchy(eps)
            etai=min(max(eta(i),-thresh),thresh)
            mu(i)=mi(i)*pcauchy(etai)
        end do
end select
end subroutine mu_mk
!
subroutine dmu_de_mk(linkid,n,mi,eta,dmu_de)
integer :: linkid,n,i
double precision :: mi(n),eta(n),dmu_de(n),etai
double precision, parameter :: eps=2.d0**(-52)
double precision, external :: dnorm,dcauchy
select case (linkid)
    case (1)
        dmu_de=1.d0
    case (2)
        forall(i=1:n) dmu_de(i)=mi(i)*max(exp(eta(i)),eps)
    case (3)
        dmu_de=-1.d0/eta**2
    case (4)
        dmu_de=2.d0*eta
    case (5)
        do i=1,n
            etai=min(eta(i),700.d0)
            dmu_de(i)=mi(i)*max(exp(etai-exp(etai)),eps)
        end do
    case (6)
        do i=1,n
            dmu_de(i)=mi(i)*max(dnorm(eta(i)),eps)
        end do
    case (7)
        do i=1,n
            dmu_de(i)=mi(i)*max(dcauchy(eta(i)),eps)
        end do
end select
end subroutine dmu_de_mk
!
subroutine d2mu_de2_mk(linkid,n,mi,eta,d2mu_de2)
integer :: linkid,n,i
double precision :: mi(n),eta(n),d2mu_de2(n)
double precision, parameter :: eps=2.d0**(-52)
double precision, external :: dnorm,dcauchy
select case (linkid)
    case (1)
        d2mu_de2=0.d0
    case (2)
        forall(i=1:n) d2mu_de2(i)=mi(i)*max(exp(eta(i)),eps)
    case (3)
        d2mu_de2=2.d0/eta**3
    case (4)
        d2mu_de2=2.d0
    case (5)
        d2mu_de2=mi*(1.d0-exp(eta))*exp(eta-exp(eta))
    case (6)
        do i=1,n
            d2mu_de2(i)=-mi(i)*eta(i)*max(dnorm(eta(i)),eps)
        end do
    case (7)
        do i=1,n
            d2mu_de2(i)=-2.d0*mi(i)*eta(i)*max(dcauchy(eta(i))/(1.d0+eta(i)**2),eps)
        end do
end select
end subroutine d2mu_de2_mk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! step_size subroutines
! canonical link function -> step_size_c(n,g,g0,p,nav,Xa,Xac,X2ac,dba,dmu_dth,d2mu_dth2,sqrt_i_bac,wac,ruac,dg_max,ai,dg,conv)
! generic link function   -> step_size_g(n,g,g0,p,nav,Xa,Xac,X2ac,dba,wght1,wght2,sqrt_i_bac,wac,ruac,dg_max,ai,dg)
subroutine step_size_c(n,g,g0,p,nav,Xa,Xac,X2ac,dba,dmu_dth,d2mu_dth2,sqrt_i_bac,wac,ruac,dg_max,ai,dg)
integer    ::    n,p,nav,k,h,ai
double precision :: g,g0,dg,Xa(n,nav),Xac(n,p-nav),dmu_dth(n),d2mu_dth2(n),sqrt_i_bac(p-nav),i_bac(p-nav)
double precision :: X2ac(n,p-nav),wac(p-nav),ruac(p-nav),dg_max,druac,druac_tmp,dg_opt,dba(0:nav)
i_bac = sqrt_i_bac**2
dg = g
do    k=1, p-nav
!    if(abs(ruac(k)).ne.g) then
        druac_tmp = 0.5d0 * ruac(k) / i_bac(k)
        druac = -dba(0)*(wac(k)*dot_product(Xac(:,k),dmu_dth)/sqrt_i_bac(k) &
            +druac_tmp*dot_product(X2ac(:,k),d2mu_dth2))
        do h = 1, nav
            druac = druac-dba(h)*(wac(k)*dot_product(Xac(:,k),dmu_dth*Xa(:,h))/sqrt_i_bac(k) &
                +druac_tmp*dot_product(X2ac(:,k),d2mu_dth2*Xa(:,h)))
        end do
        dg_opt = (g - ruac(k)) / (1.d0 - druac)
        if(dg_opt.le.0.d0 .or. dg_opt.ge.g) dg_opt = (g + ruac(k) )/ (1.d0 + druac)
!        dg_opt = (g - ruac(k))/ (1.d0 - druac)
!        if(dg_opt<0.d0 .or. dg_opt>g) then
!            dg_opt = (g + ruac(k) )/ (1.d0 + druac)
!        end if
        if(dg_opt.lt.dg .and. dg_opt .gt. 0.d0) then
            dg = dg_opt
            ai = k
        end if
!    end if
end do
if(dg_max>0.d0 .and. dg>dg_max) then
    dg = dg_max
    ai = 0
end if
if(dg > g - g0) then
    dg = g - g0
    ai = 0
end if
end subroutine step_size_c
!
subroutine step_size_g(n,g,g0,p,nav,Xa,Xac,X2ac,dba,wght1,wght2,sqrt_i_bac,wac,ruac,dg_max,ai,dg)
integer :: n,p,nav,k,h,ai
double precision :: g,g0,dg,Xa(n,nav),Xac(n,p-nav),wght1(n),wght2(n),sqrt_i_bac(p-nav),i_bac(p-nav)
double precision :: X2ac(n,p-nav),wac(p-nav),ruac(p-nav),dg_max,druac,druac_tmp,dg_opt,dba(0:nav)
i_bac = sqrt_i_bac**2
dg = g
do    k = 1, p-nav
!    if(abs(ruac(k)).ne.g) then
        druac_tmp = 0.5d0*ruac(k)/i_bac(k)
        druac = -dba(0)*(wac(k)*dot_product(Xac(:,k),wght1)/sqrt_i_bac(k) &
            +druac_tmp*dot_product(X2ac(:,k),wght2))
        do h = 1, nav
            druac = druac-dba(h)*(wac(k)*dot_product(Xac(:,k),wght1*Xa(:,h))/sqrt_i_bac(k) &
            +druac_tmp*dot_product(X2ac(:,k),wght2*Xa(:,h)))
        end do
        dg_opt = (g - ruac(k)) / (1.d0 - druac)
        if(dg_opt.le.0.d0 .or. dg_opt.ge.g) dg_opt=(g+ruac(k))/(1.d0+druac)
        if(dg_opt.lt.dg .and. dg_opt .gt. 0.d0) then
            dg = dg_opt
            ai = k
        end if
!    end if
end do
if(dg_max>0.d0 .and. dg>dg_max) then
    dg = dg_max
    ai = 0
end if
if(dg > g - g0) then
    dg = g - g0
    ai = 0
end if
end subroutine step_size_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fisher information subroutines
subroutine sqrt_i_b_mk(n,p,X2,wgh,sqrt_ib)
integer    ::    n,p,j
double precision :: X2(n,p),wgh(n),sqrt_ib(p)
forall(j=1:p) sqrt_ib(j)=sqrt(sum(wgh*X2(:,j)))
end subroutine sqrt_i_b_mk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rao subroutines:
! canonical link function -> rao_c(n,p,X,y,w,mu,sqrt_ib,ruv)
! generic link function   -> rao_g(n,p,X,y,w,mu,dth_de,sqrt_ib,ruv)
subroutine rao_c(n,p,X,y,w,mu,sqrt_ib,ruv)
integer :: n,p,j
double precision :: X(n,p),y(n),w(p),mu(n),sqrt_ib(p),ruv(p),r(n)
r=y-mu
forall(j=1:p) ruv(j)=w(j)*dot_product(r,X(:,j))/sqrt_ib(j)
end subroutine rao_c
!
subroutine rao_g(n,p,X,y,w,mu,dth_de,sqrt_ib,ruv)
integer :: n,p,j
double precision :: X(n,p),y(n),w(p),mu(n),dth_de(n),sqrt_ib(p),ruv(p),dl_de(n)
dl_de=dth_de*(y-mu)
forall(j=1:p) ruv(j)=w(j)*dot_product(dl_de,X(:,j))/sqrt_ib(j)
end subroutine rao_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shift_A(p,A,nav,ai,action)
integer    :: nav,ai,p,tmp,A(p),action
if(action.eq.1) then
    tmp = A(nav+1)
    A(nav+1) = A(nav+ai)
    A(nav+ai) = tmp
end if
if(action .eq. -1) then
    tmp = A(ai)
    A(ai) = A(nav)
    A(nav) = tmp
end if
end subroutine shift_A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine    solve(nba,Drua,dba,conv)
integer    ::    nba,ipiv(nba),conv
integer,    parameter    ::    p=1
double precision    :: Drua(nba,nba),dba(nba,p)
call dgesv(nba,p,Drua,nba,ipiv,dba,nba,conv)
if(conv.ne.0) then
    conv=1
    return
end if
end subroutine solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! jacob subroutines:
! canonical link function -> jacob_c(n,nav,Xa,X2a,nup,dmu_dth,d2mu_dth2,sqrt_i_ba,wa,rua,Drua)
! generic link function   -> jacob_g(n,nav,Xa,X2a,nup,wght1,wght2,sqrt_i_ba,wa,rua,Drua)
subroutine jacob_c(n,nav,Xa,X2a,nup,dmu_dth,d2mu_dth2,sqrt_i_ba,wa,rua,Drua)
integer :: n,nav,nup,h,k
double precision :: Xa(n,nav),X2a(n,nav),dmu_dth(n),d2mu_dth2(n),sqrt_i_ba(nav),wa(nav)
double precision :: rua(nav),Drua(0:nav,0:nav),Ikh,i_ba(1:nav)
i_ba=sqrt_i_ba**2
Drua(0,0)=sum(dmu_dth)
do h=1,nav
    Drua(0,h)=sum(dmu_dth*Xa(:,h))
    Drua(h,0)=Drua(0,h)
    if(h.gt.nup) Drua(h,0)=wa(h)*Drua(h,0)/sqrt_i_ba(h)+0.5d0*rua(h)/i_ba(h)*dot_product(X2a(:,h),d2mu_dth2)
end do
if (nav>1) then
    do k=1,nav-1
        Drua(k,k)=sum(dmu_dth*X2a(:,k))
        if(k.gt.nup) Drua(k,k)=wa(k)*Drua(k,k)/sqrt_i_ba(k)+0.5d0*rua(k)/i_ba(k)*dot_product(X2a(:,k),d2mu_dth2*Xa(:,k))
        do h=k+1,nav
            Ikh=dot_product(Xa(:,k),dmu_dth*Xa(:,h))
            if(k.gt.nup) then
                Drua(k,h)=wa(k)*Ikh/sqrt_i_ba(k)+0.5d0*rua(k)/i_ba(k)*dot_product(X2a(:,k),d2mu_dth2*Xa(:,h))
            else
                Drua(k,h)=Ikh
            end if
            if(h.gt.nup) then
                Drua(h,k)=wa(h)*Ikh/sqrt_i_ba(h)+0.5d0*rua(h)/i_ba(h)*dot_product(X2a(:,h),d2mu_dth2*Xa(:,k))
            else
                Drua(h,k)=Ikh
            end if
        end do
    end do
end if
Drua(nav,nav)=sum(dmu_dth*X2a(:,nav))
if(nav.gt.nup) Drua(nav,nav)=wa(nav)*Drua(nav,nav)/sqrt_i_ba(nav)+0.5d0*rua(nav)/i_ba(nav)*dot_product(X2a(:,nav),&
    d2mu_dth2*Xa(:,nav))
end subroutine jacob_c
!
subroutine jacob_g(n,nav,Xa,X2a,nup,wght1,wght2,sqrt_i_ba,wa,rua,Drua)
integer :: n,nav,nup,h,k
double precision :: Xa(n,nav),X2a(n,nav),sqrt_i_ba(nav),wa(nav),rua(nav),Drua(0:nav,0:nav)
double precision :: wght1(n),wght2(n),Ikh,i_ba(1:nav)
i_ba=sqrt_i_ba**2
Drua(0,0)=sum(wght1)
do h=1,nav
    Drua(0,h)=sum(wght1*Xa(:,h))
    Drua(h,0)=Drua(0,h)
    if(h.gt.nup) Drua(h,0)=wa(h)*Drua(h,0)/sqrt_i_ba(h)+0.5d0*rua(h)/i_ba(h)*dot_product(X2a(:,h),wght2)
end do
if (nav>1) then
    do k=1,nav-1
        Drua(k,k)=sum(wght1*X2a(:,k))
        if(k.gt.nup) Drua(k,k)=wa(k)*Drua(k,k)/sqrt_i_ba(k)+0.5d0*rua(k)/i_ba(k)*dot_product(X2a(:,k),wght2*Xa(:,k))
        do h=k+1,nav
            Ikh=dot_product(Xa(:,k),wght1*Xa(:,h))
            if(k.gt.nup) then
                Drua(k,h)=wa(k)*Ikh/sqrt_i_ba(k)+0.5d0*rua(k)/i_ba(k)*dot_product(X2a(:,k),wght2*Xa(:,h))
            else
                Drua(k,h)=Ikh
            end if
            if(h.gt.nup) then
                Drua(h,k)=wa(h)*Ikh/sqrt_i_ba(h)+0.5d0*rua(h)/i_ba(h)*dot_product(X2a(:,h),wght2*Xa(:,k))
            else
                Drua(h,k)=Ikh
            end if
        end do
    end do
end if
Drua(nav,nav)=sum(wght1*X2a(:,nav))
if(nav.gt.nup) Drua(nav,nav)=wa(nav)*Drua(nav,nav)/sqrt_i_ba(nav)+0.5d0*rua(nav)/i_ba(nav)*dot_product(X2a(:,nav),&
    wght2*Xa(:,nav))
end subroutine jacob_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine phi_hat(n,mu,y,dmu_dth,nnonzero,phi)
integer :: n,nnonzero
double precision :: mu(n),dmu_dth(n),y(n),r(n),phi
r=y-mu
phi=sum(r**2/dmu_dth)/(n-nnonzero)
end subroutine phi_hat
