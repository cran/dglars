subroutine predict(familyid,linkid,n,p,X,y,mi,np,b,g_seq,ng,g,coef,phi)
integer :: familyid,linkid,n,p,np,ng,nnonzero
double precision :: X(n,p),y(n),mi(n),b(0:p,np),g_seq(np),g(ng)
double precision :: coef(0:p,ng),phi(ng)
integer :: i,m,right,left
double precision :: dlt,b_prd(0:p),etah(n),muh(n),dmu_dth(n)
do i = 1, ng
    if(any(g_seq .eq. g(i))) then
        right = count(g_seq .ge. g(i))
        b_prd = b(:, right)
    else
        right = count(g_seq .gt. g(i))
        left = right + 1
        dlt = (g(i) - g_seq(right)) / (g_seq(left) - g_seq(right))
        b_prd = b(:, right) + dlt * (b(:,left) - b(:, right))
    end if
    coef(:, i) = b_prd
    if(familyid.ne.2 .and. familyid.ne.3) then
        nnonzero = 1
        etah = b_prd(0)
        do m = 1, p
            if(abs(b_prd(m)) .gt. 0.d0) then
                nnonzero = nnonzero + 1
                etah = etah + X(:, m) * b_prd(m)
            end if
        end do
        call mu_mk(linkid, n, etah, mi, muh)
        if(familyid.eq.1) then
            phi(i) = sum((y - muh)**2) / (n - nnonzero)
        else
            select case (familyid)
                case (4)
                    call dmu_dth_mk_gamma(n, muh, dmu_dth)
                case (5)
                    call dmu_dth_mk_invgaus(n, muh, dmu_dth)
            end select
            call phi_hat(n, muh, y, dmu_dth, nnonzero, phi(i))
        end if
    end if
end do
end subroutine predict
