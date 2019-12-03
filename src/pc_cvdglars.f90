!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used for the cross-validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pc_cvdglars(familyid,linkid,n,p,X,y,mi,nup,A,w,foldid,nfold,ng,g,b,phi,dev_m,dev_v,g_hat,nv,mthd,g0,dg_max,&
eps,np,ncrct,cf,NReps,nNR,conv)
    integer :: familyid,linkid,n,p,nup,A(p),foldid(n),nfold,ng,nv,mthd,np,ncrct,nNR,conv(nfold+1)
    double precision :: X(n,p),y(n),mi(n),w(0:p),g(ng),b(0:p),phi,dev_m(ng),dev_v(ng),g_hat,g0,dg_max,eps,cf,NReps
    ! internal variables
    integer :: i,lfold,np_cv,A_cv(p),nav,g_id,nnonzero(np),denom
    double precision :: b_cv(0:p,np),w_cv(0:p),phi_cv(np),ru_cv(p,np),dev(np),g_seq(np),dev_cv(ng)
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
        g_seq = 0.d0
        np_cv = np
        select case (familyid)
            case (1) !gaussian
                if(linkid.eq.1) then
                    call pc_gaussian_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,phi_cv,&
                                        ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,eps,np_cv,conv(i))
                else
                    call pc_gaussian_g(linkid,n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,&
                                        phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,ncrct,&
                                        cf,NReps,nNR,conv(i))
                end if
            case (2) !binomial
                if(linkid.eq.8) then
                    call pc_bin_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),mi(foldid(lfold+1:n)),&
                                        nup,w_cv,b_cv,phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,&
                                        eps,np_cv,ncrct,cf,NReps,nNR,conv(i))
                else
                    call pc_bin_g(linkid,n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),mi(foldid(lfold+1:n)),&
                                        nup,w_cv,b_cv,phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,&
                                        eps,np_cv,ncrct,cf,NReps,nNR,conv(i))
                end if
            case (3) !poisson
                if(linkid.eq.2) then
                    call pc_pois_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,phi_cv,&
                                        ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,ncrct,cf,&
                                        NReps,nNR,conv(i))
                else
                    call pc_pois_g(linkid,n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,&
                                        phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,&
                                        ncrct,cf,NReps,nNR,conv(i))
                end if
            case (4) !Gamma
                if(linkid.eq.3) then
                    call pc_gamma_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,phi_cv,&
                                        ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,ncrct,cf,&
                                        NReps,nNR,conv(i))
                else
                    call pc_gamma_g(linkid,n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,&
                                        phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,&
                                        ncrct,cf,NReps,nNR,conv(i))
                end if
            case (5) !inverse gaussian
                if(linkid.eq.9) then
                    call pc_invgaus_c(n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,phi_cv,&
                                        ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,ncrct,cf,&
                                        NReps,nNR,conv(i))
                else
                    call pc_invgaus_g(linkid,n-lfold,p,X(foldid(lfold+1:n),:),y(foldid(lfold+1:n)),nup,w_cv,b_cv,&
                                        phi_cv,ru_cv,dev,A_cv,nv,nav,nnonzero,g_seq,mthd,g0,2.d0,dg_max,eps,np_cv,&
                                        ncrct,cf,NReps,nNR,conv(i))
                end if
        end select
        if(conv(i).eq.0) then
            call linterpol(familyid,linkid,lfold,p,X(foldid(1:lfold),:),y(foldid(1:lfold)),mi(foldid(1:lfold)),np_cv,b_cv,g_seq,&
                        ng,g,dev_cv)
            denom = denom + 1
            dev_m = dev_m + dev_cv
            dev_v = dev_v + dev_cv**2
        end if
        call shift_vec(n,foldid,lfold)
        !foldid = cshift(foldid, lfold)
    end do
    if(denom.le.1) return
    dev_m = dev_m / denom
    dev_v = (dev_v / denom - dev_m**2) * denom / (denom - 1)
    b_cv = 0.d0
    phi_cv = 0.d0
    ru_cv = 0.d0
    dev = 0.d0
    nav = 0
    nnonzero = 0
    g_seq = 0.d0
    g_id = minloc(dev_m, dim = 1)
    g_hat = g(g_id)
    select case (familyid)
        case (1) !gaussian
            if(linkid.eq.1) then
                call pc_gaussian_c(n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,eps,np,&
                                conv(nfold+1))
            else
                call pc_gaussian_g(linkid,n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,&
                                eps,np,ncrct,cf,NReps,nNR,conv(nfold + 1))
            end if
        case (2) !binomial
            if(linkid.eq.8) then
                call pc_bin_c(n,p,X,y,mi,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,&
                                ncrct,cf,NReps,nNR,conv(nfold + 1))
            else
                call pc_bin_g(linkid,n,p,X,y,mi,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,&
                                np,ncrct,cf,NReps,nNR,conv(nfold + 1))
            end if
        case (3) !poisson
            if(linkid.eq.2) then
                call pc_pois_c(n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,&
                                ncrct,cf,NReps,nNR,conv(nfold + 1))
            else
                call pc_pois_g(linkid,n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,&
                                ncrct,cf,NReps,nNR,conv(nfold + 1))
            end if
        case (4) !Gamma
            if(linkid.eq.3) then
                call pc_gamma_c(n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,&
                                ncrct,cf,NReps,nNR,conv(nfold + 1))
            else
                call pc_gamma_g(linkid,n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,&
                                np,ncrct,cf,NReps,nNR,conv(nfold + 1))
            end if
        case (5) !inverse gaussian
            if(linkid.eq.9) then
                call pc_invgaus_c(n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,np,&
                                ncrct,cf,NReps,nNR,conv(nfold + 1))
            else
                call pc_invgaus_g(linkid,n,p,X,y,nup,w,b_cv,phi_cv,ru_cv,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,dg_max,eps,&
                                np,ncrct,cf,NReps,nNR,conv(nfold + 1))
            end if
    end select
    b = b_cv(:,np)
    phi = phi_cv(np)
    g_hat = g0
    g(1) = g_seq(1)
end subroutine pc_cvdglars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linterpol(familyid,linkid,n,p,X,y,mi,np,b,g_seq,ng,g,dev)
    integer :: familyid,linkid,n,p,np,ng
    double precision :: X(n,p),y(n),mi(n),b(0:p,np),g_seq(np),g(ng),dev(ng)
    integer :: i,m,right,left,model,check
    double precision :: g_min,g_max,g_frac(np),b_prd(0:p),dlt,eta(n),mu(n)
    model = 10 * familyid + linkid
    g_min = minval(g_seq)
    g_max = maxval(g_seq)
    g_frac = (g_seq - g_min) / (g_max - g_min)
    do i = 1, ng
        if(any(g_frac .eq.g(i))) then
            right = count(g_frac .ge. g(i))
            b_prd = b(:, right)
        else
            right = count(g_frac.gt.g(i))
            left = right + 1
            dlt = (g(i) - g_frac(right)) / (g_frac(left) - g_frac(right))
            b_prd = b(:,right) + dlt * (b(:, left) - b(:, right))
        end if
        eta = b_prd(0)
        do m = 1, p
            if(abs(b_prd(m)).gt.0.d0) eta = eta + X(:, m) * b_prd(m)
        end do
        select case (model)
            case (11)
                mu = eta
            case (28)
                call mu_mk_bin(n, eta, mi, mu)
            case (32)
                call mu_mk_pois(n, eta, mu)
            case (43)
                call mu_mk_gamma(n, eta, mu, check)
            case (59)
                call mu_mk_invgaus(n, eta, mu, check)
            case default
                call mu_mk(linkid, n, eta, mi, mu)
        end select
        select case (familyid)
            case (1)
                call deviance_gaussian(n, y, mu, dev(i))
            case (2)
                call deviance_bin(n, y, mi, mu, dev(i))
            case (3)
                call deviance_pois(n, y, mu, dev(i))
            case (4)
                call deviance_gamma(n, y, mu, dev(i))
            case (5)
                call deviance_invgaus(n, y, mu, dev(i))
            end select
    end do
end subroutine linterpol
