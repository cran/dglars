subroutine pc_gaussian_c(n,p,X,y,nup,w,b,phi,ru,dev,A,nv,nav,nnonzero,g_seq,mthd,g0,g_hat,eps,np,conv)
integer :: n,p,nup,A(p),nv,nav,nnonzero(np),mthd,np,conv
double precision :: X(n,p),y(n),w(0:p),b(0:p,np),phi(np),ru(p,np),dev(np),g_seq(np),g0,g_hat,eps
integer :: i,j,k,nstp,ai,ipiv(p),As(p)
double precision :: ym,xm(p),sqrt_Im(p),Z(n,p),yc(n),mu(n),dmu(n),ZtZ(p,p),Ztyc(p),ruv(p),g,work,v(p)
double precision :: druv(p),dg,dg_in,dg_out,ba(p),dba(p)
logical :: final
double precision, dimension(:), allocatable :: vec
double precision, dimension(:,:), allocatable :: ZatZa
final=.false.
As = 0
ba = 0.d0
ZtZ = 0.d0
mu = 0.d0
ym = sum(y)/n
yc = y - ym
do i = 1, p
    xm(i) = sum(X(:,i))/n
    sqrt_Im(i) = sqrt(dot_product(X(:,i), X(:,i)))
    Z(:,i) = (X(:,i) - xm(i)) / sqrt_Im(i)
    Ztyc(i) = dot_product(Z(:,i), yc)
end do
if(w(1).ne.0.d0) then
    w(0) = 1.d0
    forall(j=1:p) w(j) = 0.5d0 * (sqrt_Im(j) * w(j))**2
else
    w = 1.d0
end if
v = 1.d0 / w(1:p)
nstp = 1
b(0, nstp) = ym
call deviance_gaussian(n, yc, mu, dev(nstp))
if(nup.ne.0) then
    allocate(vec(nup), ZatZa(nup, nup), stat = conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
    v(A(1:nup)) = 0.d0
    As(1:nup) = A(1:nup)
    do j = 1, nup
        do i = 1, j
            ZtZ(As(i),As(j)) = dot_product(Z(:,As(i)), Z(:,As(j)))
        end do
    end do
    vec = Ztyc(As(1:nup))
    ZatZa = ZtZ(As(1:nup), As(1:nup))
    ipiv = 0
    call dsysv('U', nup, 1, ZatZa, nup, ipiv(1:nup), vec, nup, work, 1, conv)
    if(conv.ne.0) then
        conv = 5
        np=nstp
        return
    end if
    ba(As(1:nup)) = vec
    b(As(1:nup),nstp) = vec / sqrt_Im(As(1:nup))
    b(0,nstp) = ym - dot_product(xm(As(1:nup)), b(As(1:nup),nstp))
    do i = 1, nup
        mu = mu + Z(:, As(i)) * ba(As(i))
    end do
    call deviance_gaussian(n, yc, mu, dev(nstp))
    deallocate(vec, ZatZa, stat = conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
end if
nnonzero(nstp) = nup + 1
nav = nup
phi(nstp) = dev(nstp) / (n - nnonzero(nstp))
forall(i=1:p) ruv(i) = w(i) * (Ztyc(i) - dot_product(Z(:,i), mu))
ru(:,nstp) = ruv
ai = maxloc(abs(ruv(A((nup+1):p))), dim = 1)
call shift_A(p,A,nav,ai,1)
nav = nav + 1
if(nup.eq.0) then
    As(1) = A(1)
else
    i = count(As(1:(nav-1)) .lt. A(nav))
    if(i.lt.nav-1) As((i+2):nav) = As((i+1):(nav-1))
    As(i+1) = A(nav)
end if
g = abs(ruv(A(nup+1)))
g_seq(nstp) = g
if(g.le.g0 .or. g_hat .eq. 1.d0) then
    np=nstp
    return
end if
if(g_hat.ne.2.d0) g0 = g0 + g_hat * (g - g0)
do k = 1, np
    allocate(vec(nav), ZatZa(nav,nav), stat = conv)
    if(conv.ne.0) then
        conv=4
        np=nstp
        return
    end if
    if(ZtZ(A(nav), A(nav)).eq.0.d0) then
        do i = 1, nav
            if(A(i) .le. A(nav)) then
                ZtZ(A(i),A(nav)) = dot_product(Z(:,A(i)), Z(:,A(nav)))
            else
                ZtZ(A(nav),A(i)) = dot_product(Z(:,A(i)), Z(:,A(nav)))
            end if
        end do
    end if
    vec = - v(As(1:nav)) * sign(1.d0, ruv(As(1:nav)))
    ZatZa = ZtZ(As(1:nav), As(1:nav))
    ipiv = 0
    call dsysv('U', nav, 1, ZatZa, nav, ipiv(1:nav), vec, nav, work, 1, conv)
    if(conv.ne.0) then
        conv = 5
        np=nstp
        return
    end if
    dba(As(1:nav)) = vec
    dmu = 0.d0
    do i = 1, nav
        dmu = dmu + Z(:, As(i)) * vec(i)
    end do
    dg = g - g0
    if(.not.final) then
        ai = 1
        do i = nav + 1, p
            druv(A(i)) = w(A(i)) * dot_product(Z(:, A(i)), dmu)
            dg_in = (g - ruv(A(i))) / (1.d0 + druv(A(i)))
            if(dg_in .le. 0.d0 .or. dg_in .gt. g) dg_in = (g + ruv(A(i))) / (1.d0 - druv(A(i)))
            if(dg_in .lt. dg .and. dg_in .gt. 0.d0) then
                dg = dg_in
                ai = i - nav
            end if
        end do
    end if
    if(mthd .eq. 1) then
        do i = nup + 1, nav
            if(abs(ba(A(i))) .ne. 0.d0) then
                dg_out = ba(A(i)) / dba(A(i))
                if(dg_out .gt. 0.d0 .and. dg_out .le. dg) then
                    dg = dg_out
                    ai = -i
                end if
            end if
        end do
    end if
    g = g - dg
    ba(As(1:nav)) = ba(As(1:nav)) - dg * dba(As(1:nav))
    mu = mu - dg * dmu
    if(nup.gt.0) ruv(A(1:nup)) = 0.d0
    forall(i = (nup + 1):nav) ruv(A(i)) = sign(g, ruv(A(i)))
    forall(i = (nav + 1):p) ruv(A(i)) = ruv(A(i)) + dg * druv(A(i))
    nstp = nstp + 1
    b(As(1:nav), nstp) = ba(As(1:nav)) / sqrt_Im(As(1:nav))
    b(0, nstp) = ym - dot_product(xm(As(1:nav)), b(As(1:nav), nstp))
    ru(:,nstp) = ruv
    call deviance_gaussian(n, yc, mu, dev(nstp))
    nnonzero(nstp) = nav + 1
    phi(nstp) = dev(nstp) / (n - nnonzero(nstp))
    g_seq(nstp) = g
    if(mthd .eq. 1 .and. ai .lt. 0) then
        i = abs(ai)
        ba(A(i)) = 0.d0
        b(A(i), nstp) = 0.d0
        nnonzero(nstp) = nnonzero(nstp) - 1
        phi(nstp) = dev(nstp) / (n - nnonzero(nstp))
        call shift_A(p, A, nav, i, -1)
        nav = nav - 1
        i = count(As(1:(nav + 1)) .lt. A(nav + 1))
        As((i + 1):nav) = As((i + 2):(nav + 1))
        As(nav + 1) = 0
        final = .false.
    end if
    if(.not.final .and. ai .gt. 0) then
        call shift_A(p, A, nav, ai, 1)
        nav = nav + 1
        i = count(As(1:(nav-1)) .lt. A(nav))
        if(i .lt. nav - 1) As((i+2):nav) = As((i+1):(nav-1))
        As(i + 1) = A(nav)
        if(nav .ge. nv) final = .true.
    end if
    deallocate(vec, ZatZa, stat = conv)
    if(conv.ne.0) then
        conv = 4
        np=nstp
        return
    end if
    if(abs(g - g0) .le. eps) exit
end do
np = nstp
end subroutine pc_gaussian_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines related to the gaussian family
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mu_mk_gaussian(n,eta,mu)
integer :: n
double precision :: eta(n),mu(n)
mu=eta
end subroutine mu_mk_gaussian

subroutine deviance_gaussian(n, y, mu, dev)
integer :: n
double precision :: y(n), mu(n), dev, r(n)
r = y - mu
dev = sum(r**2)
end subroutine deviance_gaussian
