subroutine pc_bin_c(n,p,X,y,mi,nup,w,b,phi,ru,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
double precision, parameter :: dzero=1.0e-8,zero=1.0e-5
integer :: n,p,nup,nv,nav,np,ncrct,nNR,conv,A(p),nnonzero(np),mthd,ai,ii,i,j,k,nstp,nsbstp
double precision :: X(n,p),y(n),mi(n),w(0:p),dev(np),g_seq(np),eps,b(0:p,np),phi(np),ru(p,np),dg_max,g0,g_hat,NReps,cf
double precision :: dg,dg_tmp1,dg_tmp2,g
double precision :: eta(n),mu(n),ruv_old(p),ruv_new(p),dmu_dth(n),sqrt_ib(p),X2(n,p)
double precision, dimension(:), allocatable :: ba,dba,ba_crct,dabsruac
logical :: test(3),action,final
final = .false.
dg = 0.d0
X2 = X*X
call w_mk_bin_c(n,p,mi,X,X2,w)
nstp = 1
mu = sum(y)/sum(mi)
b(0,nstp) = log(mu(1) / (1.d0 - mu(1)))
if(b(0,nstp).eq.0.d0) then
    nnonzero(nstp)=0
else
    nnonzero(nstp)=1
end if
mu = mi * mu
eta=b(0,nstp)
call deviance_bin(n,y,mi,mu,dev(nstp))
if(nup.ne.0) then
    allocate(ba(0:nup),stat=conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
    ba=0.d0
    ba(0)=b(0,nstp)
    call bastart_bin_c(n,nup,X(:,A(1:nup)),X2(:,A(1:nup)),y,mi,NReps,nNR,ba,conv)
    if(conv.eq.0) then
        b(0,nstp)=ba(0)
        b(A(1:nup),nstp)=ba(1:nup)
        call eta_mk(n,nup,X(:,A(1:nup)),ba,eta)
        call mu_mk_bin(n,eta,mi,mu)
        call deviance_bin(n,y,mi,mu,dev(nstp))
        nnonzero(nstp)=nup+1
    else
        nup=0
    end if
    deallocate(ba,stat=conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
end if
phi(nstp)=1.d0
nav=nup
call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
call sqrt_i_b_mk(n,p,X2,dmu_dth,sqrt_ib)
call rao_c(n,p,X,y,w(1:p),mu,sqrt_ib,ruv_new)
ru(:,nstp)=ruv_new
ai=maxloc(abs(ruv_new(A((nup+1):p))),dim=1)
call shift_A(p,A,nav,ai,1)
g=abs(ruv_new(A(nup+1)))
g_seq(nstp)=g
if(g.le.g0.or.g_hat.eq.1.) then
    np=nstp
    return
end if
if(g_hat.ne.2.) g0=g0+g_hat*(g-g0)
nav = nav + 1
action = .true.
do k = 2, np
    test = .false.
    ruv_old = ruv_new
    if(action) then
        allocate(ba(0:nav),dba(0:nav),ba_crct(0:nav),dabsruac(1:(p-nav)),stat=conv)
        if(conv.ne.0) then
            conv=4
            np=nstp
            return
        end if
        action=.false.
    end if
    ba(0) = b(0,nstp)
    ba(1:nav) = b(A(1:nav),nstp)
    call prd_bin_c(mthd,g,g0,n,p,X,X2,A,nav,nup,ba,mi,mu,dmu_dth,sqrt_ib,w(1:p),ruv_old,dg_max,dba,dg,conv,ai,final)
    if(conv.ne.0) then
        np=nstp
        return
    end if
    do nsbstp=1,ncrct
        call crct_bin_c(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),y,nup,ba,dba,g,dg,w(A(1:nav)),ruv_old(A(1:nav)),NReps,&
                nNR,mi,mu,dmu_dth,ba_crct,conv)
        if(conv.ne.0) then
            conv = 0
            dg = dg * cf
            if(dg .le. dzero) then
                conv = 2
                np = nstp
                return
            end if
        else
            call sqrt_i_b_mk(n,p,X2,dmu_dth,sqrt_ib)
            call rao_c(n,p,X,y,w(1:p),mu,sqrt_ib,ruv_new)
            dg_tmp1 = dg
            test(1) = .false.
            if(.not.final) then
                do ii = nav + 1, p
                    dabsruac(ii - nav) = abs(ruv_new(A(ii))) - g + dg
                    if(dabsruac(ii - nav).gt.eps) then
                        test(1) = .true.
                        if(ruv_new(A(ii)).gt.0.d0) then
                            dg_tmp1=min(dg_tmp1,dg*(ruv_old(A(ii))-g)/(ruv_old(A(ii))-ruv_new(A(ii))-dg))
                        else
                            dg_tmp1=min(dg_tmp1,dg*(ruv_old(A(ii))+g)/(ruv_old(A(ii))-ruv_new(A(ii))+dg))
                        end if
                    end if
                end do
            end if
            dg_tmp2 = dg
            test(2) = .false.
            if(mthd.eq.1) then
                do ii = nup + 1, nav
                    if(ba_crct(ii)*ruv_new(A(ii)).lt.0.d0.and.abs(ba_crct(ii)).gt.zero.and.abs(ruv_new(A(ii))).gt.zero) then
                        test(2) = .true.
                        dg_tmp2 = min(dg_tmp2,dg*ba(ii)/(ba(ii)-ba_crct(ii)))
                    end if
                end do
            end if
            if(test(1).or.test(2)) then
                dg=min(dg_tmp1,dg_tmp2)
                if(dg.le.dzero) exit
                call mu_mk_bin(n,eta,mi,mu)
                call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
                call sqrt_i_b_mk(n,p,X2,dmu_dth,sqrt_ib)
                call rao_c(n,p,X,y,w(1:p),mu,sqrt_ib,ruv_old)
            else
                exit
            end if
        end if
    end do
    if(nsbstp .ge. ncrct) then
        conv = 2
        np = nstp
        return
    end if
    test(3) = abs(abs(ruv_new(A(nup + 1))) - g0).le.eps
    if(.not.test(1) .and. .not.test(2)) then
        nstp = nstp + 1
        b(0,nstp)=ba_crct(0)
        b(A(1:nav),nstp) = ba_crct(1:nav)
        phi(nstp) = 1.d0
        call eta_mk(n,nav,X(:,A(1:nav)),ba_crct,eta)
        call mu_mk_bin(n,eta,mi,mu)
        call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
        call sqrt_i_b_mk(n,p,X2,dmu_dth,sqrt_ib)
        call rao_c(n,p,X,y,w(1:p),mu,sqrt_ib,ruv_new)
        call deviance_bin(n,y,mi,mu,dev(nstp))
        if(dev(nstp).gt.dev(nstp - 1)) then
            conv = 6
            np = nstp
            return
        end if
        nnonzero(nstp) = nav + 1
        ru(:,nstp) = ruv_new
        g = abs(ruv_new(A(nup+1)))
        g_seq(nstp) = g

        if(any(abs(ba_crct((nup+1):nav)).le.zero).and.mthd.eq.1.and.ai.lt.0) then
            ai = abs(ai)
            ruv_new(A(ai)) = sign(g, ruv_new(A(ai)))
            b(A(ai), nstp) = 0.d0
            nnonzero(nstp) = nnonzero(nstp) - 1
            call shift_A(p, A, nav, ai, -1)
            nav = nav - 1
            action = .true.
            final = .false.
        end if
        if(minval(abs(dabsruac)).le.eps.and..not.final.and.ai.gt.0) then
            ai = minloc(abs(dabsruac), dim = 1)
            call shift_A(p,A,nav,ai,1)
            nav = nav + 1
            action = .true.
            if(nav .ge. nv) final = .true.
        end if
    else
        action = .true.
        j = p-nav
        i = nav
        if(test(1)) then
            do ai = 1, j
                if(dabsruac(ai) .ge. eps) then
                    call shift_A(p,A,nav,ai+i-nav,1)
                    nav = nav + 1
                    if(nav .ge. nv) then
                        final = .true.
                        exit
                    end if
                end if
            end do
        end if
        if(test(2)) then
            do ai = nup + 1, i
                if(ba_crct(ai)*ruv_new(A(ai)).lt.0.d0.and.abs(ba_crct(ai)).gt.zero.and.abs(ruv_new(A(ai))).gt.zero) then
                    call shift_A(p, A, nav, ai, -1)
                    nav = nav - 1
                    if(nav .eq. 0) then
                        conv = 6
                        np = nstp
                        return
                    end if
                    final = .false.
                end if
            end do
        end if
    end if
    if(action) then
        deallocate(ba,dba,ba_crct,dabsruac,stat=conv)
        if(conv.ne.0) then
            conv = 4
            np = nstp
            return
        end if
    end if
    if(test(3)) exit
end do
if(k.ge.np) then
    conv = 7
    np = nstp
    return
end if
np = nstp
end subroutine pc_bin_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine used in the starting step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bastart_bin_c(n,nav,Xa,X2a,y,mi,NReps,n_step,ba,conv)
integer :: n,nav,n_step,conv,i,h,k
double precision :: Xa(n,nav),X2a(n,nav),y(n),mi(n),NReps,ba(0:nav),eta(n),mu(n),r(n),dmu_dth(n),dba(0:nav)
double precision :: Hsn(0:nav,0:nav),sum_abs_db,work
integer :: ipiv(nav + 1)
Hsn = 0.d0
do  i = 1, n_step
    call eta_mk(n,nav,Xa,ba,eta)
    call mu_mk_bin(n,eta,mi,mu)
    r = y - mu
    dba(0) = sum(r)
    forall(h = 1:nav) dba(h) = dot_product(Xa(:,h), r)
    if(sum(abs(dba)).le.NReps) exit
    call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
    Hsn(0,0) = sum(dmu_dth)
    do k = 1, nav
        Hsn(0, k) = dot_product(dmu_dth, Xa(:,k))
        do h = 1, k - 1
            Hsn(h, k) = dot_product(Xa(:,h), dmu_dth*Xa(:,k))
        end do
        Hsn(k, k) = dot_product(dmu_dth, X2a(:,h))
    end do
    ipiv = 0
    call dsysv('U', nav + 1, 1, Hsn, nav + 1, ipiv, dba, nav + 1, work, 1, conv)
    if(conv.ne.0) then
        conv = 4
        return
    end if
    sum_abs_db = sum(abs(dba))
    if(sum_abs_db.ne.sum_abs_db) then
        conv = 4
        return
    end if
    ba = ba + dba
end do
if(i.eq.n_step) then
    conv = 3
    return
end if
end subroutine bastart_bin_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the prediction step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prd_bin_c(mthd,g,g0,n,p,X,X2,A,nav,nup,ba,mi,mu,dmu_dth,sqrt_ib,w,ruv,dg_max,dba,dg,conv,ai,final)
integer :: mthd,n,p,nav,nup,conv,A(p),ai,j
double precision :: g,g0,X(n,p),X2(n,p),ba(0:nav),mi(n),mu(n),dmu_dth(n),sqrt_ib(p),w(p),ruv(p),dg_max,dg,dg_out,dba(0:nav)
double precision :: d2mu_dth2(1:n),Drua(0:nav,0:nav)
logical :: final
dba = 0.d0
dba((nup+1):nav) = sign(1.d0,ruv(A((nup+1):nav)))
call d2mu_dth2_mk_bin(n,mi,mu,dmu_dth,d2mu_dth2)
call jacob_c(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),nup,dmu_dth,d2mu_dth2,sqrt_ib(A(1:nav)),w(A(1:nav)),ruv(A(1:nav)),Drua)
call solve(nav + 1,-Drua,dba,conv)
if(conv.ne.0) then
    conv = 1
    return
end if
if(final) then
    if(dg_max.gt.0.d0) then
        dg = min(dg_max, g - g0)
    else
        dg = g-g0
    end if
else
    call step_size_c(n,g,g0,p,nav,X(:,A(1:nav)),X(:,A((nav+1):p)),X2(:,A((nav+1):p)),dba,dmu_dth,d2mu_dth2,&
        sqrt_ib(A((nav+1):p)),w(A((nav+1):p)),ruv(A((nav+1):p)),dg_max,ai,dg)
end if
if(mthd.eq.1) then
    do j = nup + 1, nav
        if(abs(ba(j)).ne.0.d0) then
            dg_out = ba(j) / dba(j)
            if(dg_out .gt. 0.d0 .and. dg_out .le. dg) then
                dg = dg_out
                ai = -j
            end if
        end if
    end do
end if
end subroutine prd_bin_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the corrector step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crct_bin_c(n,nav,Xa,X2a,y,nup,ba,dba,g,dg,wa,rua,NReps,nNR,mi,mu,dmu_dth,ba_crct,conv)
integer ::  n,nav,nup,nNR,conv
double precision :: Xa(n,nav),X2a(n,nav),y(n),g,dg,wa(nav),rua(nav),NReps,mi(n),mu(n),dmu_dth(n),va(nav)
double precision, dimension(0:nav) :: ba,dba,ba_crct,ba_prd
va = 0.d0
va((nup+1):nav) = (g - dg) * sign(1.d0,rua((nup+1):nav))
ba_prd = ba - dg * dba
call newt_bin_c(n,nav,va,Xa,X2a,y,nup,wa,NReps,nNR,mi,mu,dmu_dth,ba_prd,conv)
if(conv.ne.0) return
ba_crct = ba_prd
end subroutine crct_bin_c
!
subroutine newt_bin_c(n,nav,va,Xa,X2a,y,nup,wa,NReps,n_step,mi,mu,dmu_dth,ba_crct,conv)
integer :: n,nav,nup,n_step,conv,i,j
double precision :: wa(nav),ba_crct(0:nav),va(nav),Xa(n,nav),X2a(n,nav),y(n),NReps,mi(n),mu(n),dmu_dth(n),sum_abs_f,sum_abs_db
double precision :: dba(0:nav),eta(n),d2mu_dth2(n),sqrt_i_ba(nav),rua(nav),Drua(0:nav,0:nav),r(n)
do  i=1,n_step
    call eta_mk(n,nav,Xa,ba_crct,eta)
    call mu_mk_bin(n,eta,mi,mu)
    call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
    call sqrt_i_b_mk(n,nav,X2a,dmu_dth,sqrt_i_ba)
    call rao_c(n,nav,Xa,y,wa,mu,sqrt_i_ba,rua)
    r = y - mu
    dba(0) = sum(r)
    forall(j=1:nup) dba(j)=dot_product(Xa(:,j),r)
    forall(j=(nup+1):nav) dba(j)=rua(j)-va(j)
    sum_abs_f=sum(abs(dba))
    if(sum_abs_f.le.NReps) exit
    call d2mu_dth2_mk_bin(n,mi,mu,dmu_dth,d2mu_dth2)
    call jacob_c(n,nav,Xa,X2a,nup,dmu_dth,d2mu_dth2,sqrt_i_ba,wa,rua,Drua)
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
if(i.eq.n_step) then
    conv=2
    return
end if
end subroutine newt_bin_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines related to the binomial family
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mu_mk_bin(n,eta,mi,mu)
integer :: i,n
double precision :: eta(n),mi(n),mu(n)
double precision, parameter :: eps=2.d0**(-52)
forall(i=1:n) mu(i) = mi(i)*max(min(1.d0/(1.d0+exp(-eta(i))),1.d0-eps),eps)
end subroutine mu_mk_bin
!
subroutine dmu_dth_mk_bin(n,mi,mu,dmu_dth)
integer :: n
double precision :: mi(n),mu(n),dmu_dth(n)
dmu_dth=mu*(1.d0-mu/mi)
end subroutine dmu_dth_mk_bin
!
subroutine  d2mu_dth2_mk_bin(n,mi,mu,dmu_dth,d2mu_dth2)
integer :: n
double precision ::mi(n),mu(n),dmu_dth(n),d2mu_dth2(n)
d2mu_dth2=dmu_dth*(1.d0-2.d0*mu/mi)
end subroutine  d2mu_dth2_mk_bin
!
subroutine  d2th_dmu2_mk_bin(n,mi,mu,d2th_dmu2)
integer :: n
double precision :: mi(n),mu(n),d2th_dmu2(n)
d2th_dmu2=1.d0/(mi-mu)**2-1.d0/mu**2
end subroutine  d2th_dmu2_mk_bin
!
subroutine deviance_bin(n,y,mi,mu,dev)
integer :: n
double precision :: y(n),mi(n),mu(n),dev
dev=2.d0*(sum(y*log(y/mu),mask=y.ne.0.d0) + sum((mi-y)*log((mi-y)/(mi-mu)),mask=y.ne.mi))
end subroutine deviance_bin
!
subroutine w_mk_bin_c(n,p,mi,X,X2,w)
integer :: n,p,i
double precision :: mi(n),X(n,p),X2(n,p),w(0:p),eta(n),mu(n),dmu_dth(n)
if(w(1).ne.0.d0) then
    call eta_mk(n,p,X,w,eta)
    call mu_mk_bin(n,eta,mi,mu)
    call dmu_dth_mk_bin(n,mi,mu,dmu_dth)
    w(0) = 1.d0
    forall(i=1:p) w(i)=0.5d0*sum(dmu_dth*X2(:,i))*w(i)**2
else
    w = 1.d0
end if
end subroutine w_mk_bin_c
