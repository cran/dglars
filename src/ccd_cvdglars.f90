!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used for the cross-validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ccd_cvdglars(familyid,linkid,n,p,X,y,mi,nup,A,w,foldid,nfold,np,g_seq,g0,ng,g_cv,nstp,eps,NReps,nNR,g_hat,mthd,b,phi,&
dev_m,dev_v,conv)
integer :: familyid,linkid,n,p,nup,A(0:p),foldid(n),nfold,np,ng,nstp,nNR,mthd,conv(nfold+1)
double precision :: X(n,p),y(n),mi(n),w(0:p),g_seq(np),g0,g_cv(ng),eps,NReps,g_hat,b(0:p),phi,dev_m(ng),dev_v(ng)
! internal variables
logical :: g_flag
integer :: i,lfold,np_cv,nstp_cv,A_cv(0:p),nav,g_id,nnonzero(np),denom
double precision :: b_cv(0:p,np),w_cv(0:p),phi_cv(np),ru_cv(p,np),dev(np),dev_cv(ng)
if(g_seq(1) .eq. 0.d0) then
    g_flag = .true.
else
    g_flag = .false.
end if
denom = 0
dev_m = 0.d0
dev_v = 0.d0
lfold = n / nfold
do i = 1, nfold
    w_cv = w
    b_cv = 0.d0
    phi_cv = 0.d0
    ru_cv = 0.d0
    dev = 0.d0
    A_cv = A
    nav = 0
    nnonzero = 0
    if(g_flag) g_seq = 0.d0
    np_cv = np
    nstp_cv = nstp
    select case (familyid)
        case (2) !binomial
            call ccd_bin_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),mi(foldid(lfold+1:n)),&
                            nup,w_cv,np_cv,g0,2.d0,nstp_cv,eps,NReps,nNR,mthd,b_cv,phi_cv,ru_cv,dev,g_seq,&
                            A_cv,nnonzero,nav,conv(i))
        case (3) !poisson
            call ccd_pois_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,np_cv,g0,2.d0,&
                            nstp_cv,eps,NReps,nNR,mthd,b_cv,phi_cv,ru_cv,dev,g_seq,A_cv,nnonzero,nav,conv(i))
    end select
    if(conv(i) .eq. 0) then
        call linterpol(familyid,linkid,lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),mi(foldid(1:lfold)),np_cv,b_cv,g_seq,&
                        ng,g_cv,dev_cv)
        denom = denom + 1
        dev_m = dev_m + dev_cv
        dev_v = dev_v + dev_cv**2
    end if
    foldid = cshift(foldid, lfold)
end do
if(denom.le.1) return
dev_m = dev_m / denom
dev_v = (dev_v / denom - dev_m**2) * denom / (denom - 1)
b_cv = 0.d0
phi_cv = 0.d0
ru_cv = 0.d0
dev = 0.d0
A_cv = A
nav = 0
nnonzero = 0
g_seq = 0.d0
g_id = minloc(dev_m, dim = 1)
g_hat = g_cv(g_id)
select case (familyid)
    case (2) !binomial
        call ccd_bin_c(n,p,X,y,mi,nup,w,np,g0,g_hat,nstp,eps,NReps,nNR,mthd,b_cv,phi_cv,ru_cv,dev,g_seq,A_cv,nnonzero,nav,&
                        conv(nfold + 1))
    case (3) !poisson
        call ccd_pois_c(n,p,X,y,nup,w,np,g0,g_hat,nstp,eps,NReps,nNR,mthd,b_cv,phi_cv,ru_cv,dev,g_seq,A_cv,nnonzero,nav,&
                        conv(nfold + 1))
end select
b = b_cv(:,np)
phi = phi_cv(np)
g_hat = g0
g_cv(1) = g_seq(1)
end subroutine ccd_cvdglars
