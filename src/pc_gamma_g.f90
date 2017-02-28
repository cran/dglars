subroutine pc_gamma_g(linkid,n,p,X,y,nup,w,b,phi,ru,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,ncrct,cf,NReps,nNR,conv)
double precision, parameter :: dzero=1.0e-8,zero=1.0e-5
integer :: linkid, n,p,nup,nv,nav,np,ncrct,nNR,conv,A(p),nnonzero(np),mthd,ai,ii,i,j,k,nstp,nsbstp
double precision :: X(n,p),y(n),w(0:p),dev(np),g_seq(np),eps,b(0:p,np),phi(np),ru(p,np),dg_max,g0,g_hat,NReps
double precision :: cf,dg,dg_tmp1,dg_tmp2,g
double precision :: eta(n),mu(n),ruv_old(p),ruv_new(p),dmu_dth(n),dmu_de(n),dth_de(n),sqrt_ib(p),X2(n,p),mi(n)
double precision, dimension(:), allocatable :: ba,dba,ba_crct,dabsruac
logical :: test(3),action,final
mi = 1.d0
final = .false.
dg = 0.d0
X2 = X*X
nstp = 1
call w_mk_gamma_g(linkid,n,p,mi,X,X2,w,conv)
if(conv.eq.5) then
    np=nstp
    return
end if
mu=sum(y)/n
call linkfun(linkid,mu,b(0,nstp))
eta=b(0,nstp)
dev(nstp)=2.d0*(n*log(mu(1))-sum(log(y)))
if(nup.ne.0) then
    allocate(ba(0:nup),stat=conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
    ba=0.d0
    ba(0)=b(0,nstp)
    call bastart_gamma_g(linkid,n,nup,X(:,A(1:nup)),X2(:,A(1:nup)),y,mi,NReps,nNR,ba,conv)
    if(conv.eq.0) then
        b(0,nstp)=ba(0)
        b(A(1:nup),nstp)=ba(1:nup)
        call eta_mk(n,nup,X(:,A(1:nup)),ba,eta)
        call mu_mk(linkid,n,eta,mi,mu)
        if(any(mu.le.0.d0)) then
            conv=5
            np=nstp
            return
        end if
        call deviance_gamma(n,y,mu,dev(nstp))
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
nnonzero(nstp)=nup+1
nav=nup
call dmu_dth_mk_gamma(n,mu,dmu_dth)
call phi_hat(n,mu,y,dmu_dth,nnonzero(nstp),phi(nstp))
call dmu_de_mk(linkid,n,mi,eta,dmu_de)
dth_de=dmu_de/dmu_dth
call sqrt_i_b_mk(n,p,X2,dth_de*dmu_de,sqrt_ib)
call rao_g(n,p,X,y,w(1:p),mu,dth_de,sqrt_ib,ruv_new)
ru(:,nstp)=ruv_new
ai=maxloc(abs(ruv_new(A((nup+1):p))),dim=1)
call shift_A(p,A,nav,ai,1)
g=abs(ruv_new(A(nup+1)))
g_seq(nstp)=g
if(g.le.g0.or.g_hat.eq.1.d0) then
    np=nstp
    return
end if
if(g_hat.ne.2.d0) g0=g0+g_hat*(g-g0)
nav = nav + 1
action = .true.
do k= 2, np
    test=.false.
    ruv_old=ruv_new
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
    call prd_gamma_g(linkid,mthd,g,g0,n,p,X,X2,y,A,nav,nup,ba,mi,eta,mu,dth_de,dmu_de,dmu_dth,sqrt_ib,&
                    w(1:p),ruv_old,dg_max,dba,dg,conv,ai,final)
    if(conv.ne.0) then
        np=nstp
        return
    end if
    do nsbstp=1,ncrct
        call crct_gamma_g(linkid,n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),y,nup,ba,dba,g,dg,w(A(1:nav)),ruv_old(A(1:nav)),&
                        NReps,nNR,mi,eta,mu,dth_de,dmu_de,ba_crct,conv)
        if(conv.ne.0) then
            if(conv.eq.5) then
                np=nstp
                return
            end if
            conv=0
            dg=dg*cf
            if(dg.le.dzero) then
                conv=2
                np=nstp
                return
            end if
        else
            call sqrt_i_b_mk(n,p,X2,dth_de*dmu_de,sqrt_ib)
            call rao_g(n,p,X,y,w(1:p),mu,dth_de,sqrt_ib,ruv_new)
            dg_tmp1=dg
            test(1)=.false.
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
            dg_tmp2=dg
            test(2)=.false.
            if(mthd.eq.1) then
                do ii = nup + 1, nav
                    if(ba_crct(ii)*ruv_new(A(ii)).lt.0.d0 .and. abs(ba_crct(ii)).gt.zero .and. abs(ruv_new(A(ii))).gt.zero) then
                        test(2)=.true.
                        dg_tmp2=min(dg_tmp2,dg*ba(ii)/(ba(ii)-ba_crct(ii)))
                    end if
                end do
            end if
            if(test(1).or.test(2)) then
                dg=min(dg_tmp1,dg_tmp2)
                if(dg.le.dzero) exit
                call mu_mk(linkid,n,eta,mi,mu)
                if(any(mu.le.0.d0)) then
                    conv=5
                    np=nstp
                    return
                end if
                call dmu_dth_mk_gamma(n,mu,dmu_dth)
                call dmu_de_mk(linkid,n,mi,eta,dmu_de)
                dth_de=dmu_de/dmu_dth
                call sqrt_i_b_mk(n,p,X2,dth_de*dmu_de,sqrt_ib)
                call rao_g(n,p,X,y,w(1:p),mu,dth_de,sqrt_ib,ruv_new)
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
    test(3) = abs(abs(ruv_new(A(nup+1))) - g0) .le. eps
    if(.not.test(1) .and. .not.test(2)) then
        nstp = nstp+1
        b(0,nstp) = ba_crct(0)
        b(A(1:nav), nstp) = ba_crct(1:nav)
        call eta_mk(n,nav,X(:,A(1:nav)),ba_crct,eta)
        call mu_mk(linkid,n,eta,mi,mu)
        if(any(mu .le. 0.d0)) then
            conv = 5
            np = nstp
            return
        end if
        call dmu_dth_mk_gamma(n,mu,dmu_dth)
        call dmu_de_mk(linkid,n,mi,eta,dmu_de)
        dth_de = dmu_de / dmu_dth
        call sqrt_i_b_mk(n,p,X2,dth_de*dmu_de,sqrt_ib)
        call rao_g(n,p,X,y,w(1:p),mu,dth_de,sqrt_ib,ruv_new)
        call deviance_gamma(n,y,mu,dev(nstp))
        if(dev(nstp).gt.dev(nstp - 1)) then
            conv = 6
            np = nstp
            return
        end if
        nnonzero(nstp) = nav + 1
        call phi_hat(n,mu,y,dmu_dth,nnonzero(nstp),phi(nstp))
        ru(:, nstp) = ruv_new
        g = abs(ruv_new(A(nup + 1)))
        g_seq(nstp) = g
        if(any(abs(ba_crct((nup+1):nav)).le.zero) .and. mthd.eq.1 .and. ai.lt.0) then
            ai = abs(ai)
            ruv_new(A(ai)) = sign(g, ruv_new(A(ai)))
            b(A(ai), nstp) = 0.d0
            nnonzero(nstp) = nnonzero(nstp) - 1
            call phi_hat(n,mu,y,dmu_dth,nnonzero(nstp),phi(nstp))
            call shift_A(p, A, nav, ai, -1)
            nav = nav - 1
            action = .true.
            final = .false.
        end if
        if(minval(abs(dabsruac)) .le. eps .and. .not. final .and. ai .gt.0) then
            ai = minloc(abs(dabsruac),dim=1)
            call shift_A(p, A, nav, ai, 1)
            nav = nav + 1
            action = .true.
            if(nav .ge. nv) final = .true.
        end if
    else
        action = .true.
        j = p - nav
        i = nav
        if(test(1)) then
            do ai = 1, j
                if(dabsruac(ai).ge.eps) then
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
if(k .ge. np) then
    conv = 7
    np = nstp
    return
end if
np = nstp
end subroutine pc_gamma_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine w_mk_gamma_g(linkid,n,p,mi,X,X2,w,check)
integer :: linkid,n,p,i,check
double precision :: mi(n),X(n,p),X2(n,p),w(0:p),eta(n),mu(n),dmu_de(n),dmu_dth(n),wght(n)
if(w(1).ne.0.) then
    call eta_mk(n,p,X,w,eta)
    call mu_mk(linkid,n,eta,mi,mu)
    if(any(mu.le.0.)) then
        check=5
        return
    end if
    call dmu_dth_mk_gamma(n,mu,dmu_dth)
    call dmu_de_mk(linkid,n,mi,eta,dmu_de)
    wght=dmu_de**2/dmu_dth
    w(0) = 1.d0
    forall(i=1:p) w(i)=0.5d0*sum(wght*X2(:,i))*w(i)**2
else
    w = 1.d0
end if
end subroutine w_mk_gamma_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the starting step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bastart_gamma_g(linkid,n,nav,Xa,X2a,y,mi,NReps,n_step,ba,conv)
integer :: linkid,n,nav,n_step,conv,i,h,k
double precision :: Xa(n,nav),X2a(n,nav),y(n),mi(n),NReps,ba(0:nav),eta(n),mu(n),dl_de(n),dmu_de(n),dmu_dth(n),dba(0:nav)
double precision :: Hsn(0:nav,0:nav),sum_abs_db,dth_de(n),wght(n), work
integer :: ipiv(nav + 1)
Hsn = 0.d0
do    i = 1, n_step
    call eta_mk(n,nav,Xa,ba,eta)
    call mu_mk(linkid,n,eta,mi,mu)
    if(any(mu.le.0.d0)) then
        conv = 5
        return
    end if
    call dmu_dth_mk_gamma(n,mu,dmu_dth)
    call dmu_de_mk(linkid,n,mi,eta,dmu_de)
    dth_de = dmu_de / dmu_dth
    dl_de = dth_de * (y - mu)
    dba(0) = sum(dl_de)
    forall(h = 1:nav) dba(h) = dot_product(Xa(:,h), dl_de)
    if(sum(abs(dba)).le.NReps) exit
    wght = dth_de*dmu_de
    Hsn(0,0) = sum(wght)
    do k = 1, nav
        Hsn(0, k) = dot_product(wght, Xa(:,k))
        do h = 1, k - 1
            Hsn(h, k) = dot_product(Xa(:, h), wght*Xa(:, k))
        end do
        Hsn(k, k) = dot_product(wght, X2a(:,h))
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
end subroutine bastart_gamma_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the prediction step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prd_gamma_g(linkid,mthd,g,g0,n,p,X,X2,y,A,nav,nup,ba,mi,eta,mu,dth_de,dmu_de,dmu_dth,sqrt_ib,w,ruv,&
                        dg_max,dba,dg,conv,ai,final)
integer :: linkid,mthd,n,p,nav,nup,conv,A(p),ai,j
double precision :: g,g0,X(n,p),X2(n,p),y(n),ba(0:nav),mi(n),eta(n),mu(n),sqrt_ib(p),w(p),ruv(p),dg_max,dg,dg_out,dba(0:nav)
double precision :: dth_de(n),d2th_de2(n),dmu_de(n),d2mu_de2(n),dmu_dth(n),d2th_dmu2(n),r(n),wght1(n),wght2(n),Drua(0:nav,0:nav)
logical    :: final
dba = 0.d0
dba((nup+1):nav) = sign(1.d0, ruv(A((nup+1):nav)))
r = y - mu
call d2mu_de2_mk(linkid,n,mi,eta,d2mu_de2)
call d2th_dmu2_mk_gamma(n,mu,d2th_dmu2)
d2th_de2 = dmu_de**2 * d2th_dmu2 + d2mu_de2 / dmu_dth
wght1 = dth_de * dmu_de - d2th_de2 * r
wght2 = 2.d0 * dth_de * d2mu_de2 + dmu_de**3 * d2th_dmu2
call jacob_g(n,nav,X(:,A(1:nav)),X2(:,A(1:nav)),nup,wght1,wght2,sqrt_ib(A(1:nav)),w(A(1:nav)),ruv(A(1:nav)),Drua)
call solve(nav + 1, -Drua, dba, conv)
if(conv.ne.0)    then
    conv = 1
    return
end if
if(final) then
    if(dg_max .gt. 0.d0) then
        dg = min(dg_max, g - g0)
    else
        dg = g - g0
    end if
else
    call step_size_g(n,g,g0,p,nav,X(:,A(1:nav)),X(:,A((nav+1):p)),X2(:,A((nav+1):p)),dba,wght1,wght2,&
                    sqrt_ib(A((nav+1):p)),w(A((nav+1):p)),ruv(A((nav+1):p)),dg_max,ai,dg)
end if
if(mthd.eq.1) then
    do j = nup + 1, nav
        if(abs(ba(j)) .ne. 0.d0) then
            dg_out = ba(j) / dba(j)
            if(dg_out .gt. 0.d0 .and. dg_out .le. dg) then
                dg = dg_out
                ai = -j
            end if
        end if
    end do
end if
end subroutine prd_gamma_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used in the corrector step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crct_gamma_g(linkid,n,nav,Xa,X2a,y,nup,ba,dba,g,dg,wa,rua,NReps,nNR,mi,eta,mu,dth_de,dmu_de,ba_crct,conv)
integer :: linkid,n,nav,nup,nNR,conv
double precision :: Xa(n,nav),X2a(n,nav),y(n),g,dg,wa(nav),rua(nav),NReps,mi(n),mu(n),dth_de(n),dmu_de(n),va(nav)
double precision, dimension(0:nav) :: eta(n),ba,dba,ba_crct,ba_prd
va = 0.d0
va((nup+1):nav) = (g - dg)*sign(1.d0,rua((nup+1):nav))
ba_prd=ba-dg*dba
call newt_gamma_g(linkid,n,nav,va,Xa,X2a,y,nup,wa,NReps,nNR,mi,eta,mu,dth_de,dmu_de,ba_prd,conv)
if(conv.ne.0) return
ba_crct=ba_prd
end subroutine crct_gamma_g
!
subroutine newt_gamma_g(linkid,n,nav,va,Xa,X2a,y,nup,wa,NReps,n_step,mi,eta,mu,dth_de,dmu_de,ba_crct,conv)
integer    :: linkid,n,nav,nup,n_step,conv,i,j
double precision :: wa(nav),ba_crct(0:nav),va(nav),Xa(n,nav),X2a(n,nav),y(n),NReps,mi(n),mu(n),sum_abs_f,sum_abs_db
double precision :: dba(0:nav),eta(n),r(n),sqrt_i_ba(nav),rua(nav),Drua(0:nav,0:nav)
double precision :: dmu_dth(n),d2th_dmu2(n),dmu_de(n),d2mu_de2(n),dth_de(n),d2th_de2(n),dl_de(n),wght1(n),wght2(n)
do    i=1, n_step
    call eta_mk(n,nav,Xa,ba_crct,eta)
    call mu_mk(linkid,n,eta,mi,mu)
    if(any(mu.le.0.d0)) then
        conv=5
        return
    end if
    call dmu_dth_mk_gamma(n,mu,dmu_dth)
    call dmu_de_mk(linkid,n,mi,eta,dmu_de)
    dth_de=dmu_de/dmu_dth
    call sqrt_i_b_mk(n,nav,X2a,dth_de*dmu_de,sqrt_i_ba)
    call rao_g(n,nav,Xa,y,wa,mu,dth_de,sqrt_i_ba,rua)
    r=y-mu
    dl_de=dth_de*r
    dba(0)=sum(dl_de)
    forall(j=1:nup) dba(j)=dot_product(Xa(:,j),dl_de)
    forall(j=(nup+1):nav) dba(j)=rua(j)-va(j)
    sum_abs_f=sum(abs(dba))
    if(sum_abs_f.le.NReps) exit
    call d2mu_de2_mk(linkid,n,mi,eta,d2mu_de2)
    call d2th_dmu2_mk_gamma(n,mu,d2th_dmu2)
    d2th_de2=dmu_de**2*d2th_dmu2+d2mu_de2/dmu_dth
    wght1=dth_de*dmu_de-d2th_de2*r
    wght2=2.d0*dth_de*d2mu_de2+dmu_de**3*d2th_dmu2
    call jacob_g(n,nav,Xa,X2a,nup,wght1,wght2,sqrt_i_ba,wa,rua,Drua)
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
end subroutine newt_gamma_g







