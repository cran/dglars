function qnorm ( p )
!*****************************************************************************80
!
!! qnorm inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10**16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.  If P is outside this range,
!    an "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) D_NORMAL_01_CDF_INVERSE, the normal deviate 
!    value with the property that the probability of a standard normal 
!    deviate being less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real   ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) qnorm 
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00

  if ( p <= 0.0D+00 ) then
    qnorm = - huge ( p )
    return
  end if

  if ( 1.0D+00 <= p ) then
    qnorm = huge ( p )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    qnorm = q * r8poly_value ( 8, a, r ) &
                                 / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

!    if ( r <= 0.0D+00 ) then
!      qnorm = - 1.0D+00
!      stop
!    end if

    r = sqrt ( -log ( r ) )

    if ( r <= split2 ) then

      r = r - const2
      qnorm = r8poly_value ( 8, c, r ) &
                               / r8poly_value ( 8, d, r )

    else

      r = r - split2
      qnorm = r8poly_value ( 8, e, r ) &
                               / r8poly_value ( 8, f, r )
   
    end if

    if ( q < 0.0D+00 ) then
      qnorm = - qnorm
    end if

  end if

  return
end

function r8poly_value ( n, a, x )
!*****************************************************************************80
!
!! R8POLY_VALUE evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value
  real ( kind = 8 ) x

  r8poly_value = 0.0D+00
  do i = n, 1, -1
    r8poly_value = r8poly_value * x + a(i)
  end do

  return
end
!
function dnorm( eta )
!*****************************************************************************80
! this function computes the density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  Author:
!
!    Luigi Augugliaro
!
!  Parameters:
!
!   Input,  double precision :: eta
!   Output, double precision :: dnorm, the value of the standard normal distribution evaluaed at eta
implicit none
real (kind = 8) eta,x
real (kind = 8) dnorm
real (kind = 8), parameter :: const = 1.d0/sqrt(8.d0*atan(1.d0))
if(eta.eq.0.d0) then
    dnorm = const
    return
end if
x = abs(eta)
dnorm = const * exp(-0.5d0 * x**2)
return
end
!
function qcauchy( p )
!*****************************************************************************80
! this function computes the quantiles function of the Cauchy istribution
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  Author:
!
!    Luigi Augugliaro
!
!  Parameters:
!
!   Input,  double precision :: p
!   Output, double precision :: qcauchy, the quantile of the Cauchy distribution corresponding to p
implicit none
real ( kind = 8 ) p
real ( kind = 8 ) qcauchy
real (kind = 8), parameter :: pi = 4.d0*atan(1.d0)
if(p .eq. 0.5d0) then
    qcauchy = 0.d0
else
    qcauchy = tan(pi * (p - 0.5d0))
end if
return
end
function pcauchy( eta )
!*****************************************************************************80
! this function computes the cumulative distribution function of the Cauchy distribution
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  Author:
!
!    Luigi Augugliaro
!
!  Parameters:
!
!   Input,  double precision :: eta
!   Output, double precision :: pcauchy, the value of the cumulative distribution function evaluate at eta
implicit none
real ( kind = 8 ) eta
real ( kind = 8 ) pcauchy
real (kind = 8), parameter :: pi = 4.d0*atan(1.d0)
if(eta .eq. 0.d0) then
    pcauchy=0.5d0
else
    pcauchy=0.5d0+atan(eta)/pi
end if
return
end
!
function dcauchy( eta )
!*****************************************************************************80
! this function computes the denisty function of the Cauchy didstribution
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  Author:
!
!    Luigi Augugliaro
!
!  Parameters:
!
!   Input,  double precision :: eta
!   Output, double precision :: dcauchy, the value of the denisty function of the Cauchy didstribution
implicit none
real ( kind = 8 ) eta
real ( kind = 8 ) dcauchy
real (kind = 8), parameter :: pi = 4.d0*atan(1.d0)
dcauchy = 1.d0/(pi * (1.d0 + eta**2))
return
end








