subroutine ccd_bin_c(n,p,X,y,mi,nup,w,np,g0,g_hat,nstp,eps,NReps,nNR,mthd,b,phi,ru,dev,g_seq,A,nnonzero,nav,conv)
double precision, parameter :: big = 9.9e30, zero = 1.0e-5
integer :: n,p,nup,np,nstp,nNR,mthd,A(0:p),nnonzero(np),nav,conv
double precision :: X(n,p),y(n),mi(n),w(0:p),g0,g_hat,eps,NReps,b(0:p,np),phi(np),ru(p,np),dev(np),g_seq(np)
integer :: step,m,nav_old,nstp_ccd,ai
double precision :: X2(n,p),ob(0:p),nb(0:p),eta(n),mu(n),r(n),V(n),sV,g,Xai(n),dlai,dl(p)
double precision :: Iai,I(p),cf,dbai_max,zai,adlai,dbai,adlai_max
double precision, dimension(:), allocatable :: bstart
logical :: new_step
X2 = X**2
dl = 0.d0
ob = 0.d0
nb = 0.d0
call w_mk_bin_c(n, p, mi, X, X2, w)
mu = sum(y) / sum(mi)
ob(0) = log(mu(1) / (1.d0 - mu(1)))
eta = ob(0)
mu = mi * mu
if(nup .ne. 0) then
    allocate(bstart(0:nup), stat = conv)
    if(conv.ne.0) then
        conv = 4
        np = 1
        return
    end if
    call bastart_bin_c(n, nup, X(:, A(1:nup)), X2(:, A(1:nup)), y, mi, NReps, nNR, bstart, conv)
    if(conv .eq. 0) then
        ob(A(0:nup)) = bstart
        call eta_mk(n, nup, X(:, A(1:nup)), ob(A(0:nup)), eta)
        call mu_mk_bin(n, eta, mi, mu)
    else
        conv = 0
        nup = 0
    end if
    deallocate(bstart, stat = conv)
    if(conv.ne.0) then
        conv = 4
        np = 1
        return
    end if
end if
nav = nup
r = y - mu
V = mu * (1.d0 - mu / mi)
sV = sum(V)
do m = 1, p
    I(m) = dot_product(V, X2(:, m))
end do
if(g_seq(1) .eq. 0.d0) then
    if(np .gt. 1) then
        g = 0.d0
        do m = nav + 1, p
            ai = A(m)
            dlai = dot_product(X(:, ai), r)
            Iai = I(ai)
            g = max(g, w(ai) * abs(dlai) / sqrt(Iai))
        end do
        if(g_hat .ne. 2.d0) g0 = g0 + g_hat * (g - g0)
        cf = exp((log(g0) - log(g))/(np - 1))
        g_seq(1) = g
        do step = 2, np
            g_seq(step) = g_seq(step - 1) * cf
        end do
    else
        g_seq(1) = g0
    end if
end if
nstp_ccd = 0
do step = 1, np
    g = g_seq(step)
    do
        nstp_ccd = nstp_ccd + 1
        if(nstp_ccd .gt. nstp) then
            np = step - 1
            conv = 3
            return
        end if
        dbai_max = 0.d0
        do m = 1, nup
            ai = A(m)
            Xai = X(:, ai)
            Iai = I(ai)
            dlai = dot_product(Xai, r)
            dbai = dlai / Iai
            nb(ai) = ob(ai) + dbai
            dbai_max = max(dbai_max, abs(dbai))
            r = r - V * Xai * dbai
            ob(ai) = nb(ai)
        end do
        do m = nup + 1, nav
            ai = A(m)
            Xai = X(:, ai)
            Iai = I(ai)
            dlai = dot_product(Xai, r)
            zai = w(ai) * (dlai + Iai * ob(ai))
            adlai = abs(zai) - sqrt(Iai) * g
            nb(ai) = sign(1.d0, zai) * adlai / (w(ai) * Iai)
            dbai = nb(ai) - ob(ai)
            dbai_max = max(dbai_max, abs(dbai))
            r = r - dbai * V * Xai
            ob(ai) = nb(ai)
        end do
        dbai = sum(r) / sV
        nb(0) = ob(0) + dbai
        r = r - V * dbai
        ob(0) = nb(0)
        dbai_max = max(dbai_max, abs(dbai))
        if(dbai_max .le. eps) then
            call eta_mk(n, nav, X(:, A(1:nav)), nb(A(0:nav)), eta)
            call mu_mk_bin(n, eta, mi, mu)
            r = y - mu
            V = mu * (1.d0 - mu / mi)
            sV = sum(V)
            do m = 1, nav
                ai = A(m)
                dl(ai) = dot_product(X(:, ai), r)
                I(ai) = dot_product(V, X2(:, ai))
            end do
            adlai_max = max(0.d0, abs(sum(r)))
            do m = 1, nup
                ai = A(m)
                dlai = dl(ai)
                adlai_max = max(adlai_max, abs(dlai))
            end do
            do m = nup + 1, nav
                ai = A(m)
                dlai = dl(ai)
                Iai = I(ai)
                adlai_max = max(adlai_max, abs(abs(w(ai) * dlai) - sqrt(Iai) * g))
            end do
            if(adlai_max .le. eps) then
                new_step = .true.
                nav_old = nav
                do m = nup + 1, nav_old
                    ai = A(m)
                    dlai = dl(ai)
                    if(mthd .eq. 1 .and. abs(dlai) .gt. zero .and. abs(nb(ai)) .gt. zero) then
                        if(dlai * nb(ai) .le. 0.d0) then
                            new_step = .false.
                            nb(ai) = 0.d0
                            ob(ai) = 0.d0
                            A(m) = A(nav)
                            A(nav) = ai
                            nav = nav - 1
                            exit
                        end if
                    end if
                end do
                if(new_step) then
                    nav_old = nav
                    do m = nav_old + 1, p
                        ai = A(m)
                        Xai = X(:, ai)
                        dl(ai) = dot_product(Xai, r)
                        I(ai) = dot_product(V, X2(:, ai))
                        dlai = dl(ai)
                        Iai = I(ai)
                        if(abs(w(ai) * dlai) .gt. (sqrt(Iai) * g)) then
                            new_step = .false.
                            nav = nav + 1
                            A(m) = A(nav)
                            A(nav) = ai
                        end if
                    end do
                end if
                if(new_step) exit
            end if
        end if
    end do
    b(A(0:nav), step) = nb(A(0:nav))
    phi(step) = 1.d0
    ru(:, step) = w(1:p) * dl / sqrt(I)
    nnonzero(step) = nav + 1
    call deviance_bin(n, y, mi, mu, dev(step))
end do
end subroutine ccd_bin_c
