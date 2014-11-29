! Adapted from
!   Armstrong JQSRT 7, 85 (1967)
! I also looked at
!   http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library43.html?VOIGT
!
! Note that sometimes the Voigt function is defined as H(a, u), which is equivalent to
!       voigt_scalar(u, a)
! here.
!
! Fujun Du  2014-03-27 Thu 18:34:57
!
module voigt

implicit none

private

double precision, parameter, dimension(10) :: &
  w=(/0.462243670D0,   0.286675505D0,    0.109017206D0, 0.0248105209D0, &
      0.00324377334D0, 0.000228338636D0, 7.80255648D-6, 1.08606937D-7, &
      4.39934099D-10,  2.22939365D-13/)

double precision, parameter, dimension(10) :: &
  t=(/0.245340708D0,   0.737473729D0,    1.23407622D0,  1.73853771D0, &
      2.25497400D0,    2.78880606D0,     3.34785457D0,  3.94476404D0, &
      4.60368245D0,    5.38748089D0/)

double precision, parameter, dimension(34) :: &
  c=(/ 0.1999999999972224D0,  -0.1840000000029998D0,   0.1558399999965025D0,  &
      -0.1216640000043988D0,   0.0877081599940391D0,  -0.0585141248086907D0,  &
       0.0362157301623914D0,  -0.0208497654398036D0,   0.0111960116346270D0,  &
      -0.56231896167109D-2,    0.26487634172265D-2,   -0.11732670757704D-2,   &
       0.4899519978088D-3,    -0.1933630801528D-3,     0.722877446788D-4,     &
      -0.256555124979D-4,      0.86620736841D-5,      -0.27876379719D-5,      &
       0.8566873627D-6,       -0.2518433784D-6,        0.709360221D-7,        &
      -0.191732257D-7,         0.49801256D-8,         -0.12447734D-8,         &
       0.2997777D-9,          -0.696450D-10,           0.156262D-10,          &
      -0.33897D-11,            0.7116D-12,            -0.1447D-12,            &
       0.285D-13,             -0.55D-14,               0.10D-14,              &
      -0.2D-15/)


public :: voigt_scalar
public :: voigt_vector


contains


pure function voigt_vector(x, y, n)
  integer, intent(in) :: n
  double precision, dimension(n) :: voigt_vector
  double precision, dimension(n), intent(in) :: x
  double precision, intent(in) :: y
  integer i
  do i=1, n
    voigt_vector(i) = voigt_scalar(x(i), y)
  end do
end function voigt_vector



pure function voigt_scalar(x, y)
  double precision voigt_scalar
  double precision, intent(in) :: x, y
  !
  if (((y .lt. 1D0) .and. (x .lt. 4D0)) .or. &
      ((y .lt. 1.8D0/(1D0+x)) .and. (x .gt. 4D0))) then
    !
    voigt_scalar = k1(x, y)
    !
  else if ((y .gt. 1D0) .and. (y .le. 2.5D0) .and. (x .le. 4D0)) then
    !
    voigt_scalar = k2(x, y)
    !
  else
    !
    voigt_scalar = k3(x, y)
    !
  end if
end function voigt_scalar


pure function k1(x, y)
  double precision k1
  double precision, intent(in) :: x, y
  integer i
  double precision :: u1, bno1, bno2, bn, dno1, dno2, dn, coef, f, x1, funct, g
  double precision :: q, yn, xx, tmp
  !
  tmp = (x+y)*(x-y)
  if (tmp .ge. 222D0) then
    u1 = 0D0
  else
    u1 = exp(-tmp) * cos(2D0*y*x)
  end if
  !
  if (x .le. 5D0) then
    !Clenshaw's Algoriuhm
    bno1 = 0.0D0
    bno2 = 0.0D0
    !
    x1 = x/5D0
    !
    coef = 4D0*x1*x1 - 2D0
    !
    do i=34, 1, -1
      bn   = coef*bno1 - bno2 + c(i)
      bno2 = bno1
      bno1 = bn
    end do
    !
    f = x1 * (bn - bno2)
    dno1 = 1D0 - 2D0 * x * f
    dno2 = f
  else
    xx = 1D0/(x*x)
    dno1 = &
      -xx*(0.5D0 + xx*(0.75D0 + xx*(1.875D0 + xx*(6.5625D0 + &
           xx*(29.53125D0 + xx*(162.4218D0 + xx*1055.7421D0))))))
    dno2 = (1D0 - dno1) / (2D0 * x)
  end if
  !
  funct = y * dno1
  if (y .gt. 1D-8) then
    q = 1D0
    yn = y
    do i=2, 100
      dn = -2D0*(x*dno1 + dno2) / dble(i)
      dno2 = dno1
      dno1 = dn
      if (mod(i, 2) .eq. 1) then
         q  = -q
         yn = yn * (y*y)
         g  = dn * yn
         funct = funct + q * g
         if (abs(g/funct) .le. 1D-10) then
           exit
         end if
      end if
    end do
  end if
  k1 = u1 - 1.12837917D0 * funct
end function k1


pure function k2(x, y)
  double precision k2
  double precision, intent(in) :: x, y
  double precision r, s, g
  integer i
  g = 0D0
  do i=1, 10
    r = t(i) - x
    s = t(i) + x
    g = g + w(i) * (4D0*t(i)**2 - 2D0) * &
            (r*atan(r/y) + s * atan(s/y) - &
             0.5D0 * y * (log(y*y+r*r) + log(y*y+s*s)))
  end do
  k2 = 0.318309886 * g
end function k2


pure function k3(x, y)
  double precision k3
  double precision, intent(in) :: x, y
  double precision g
  integer i
  g = 0D0
  do i=1, 10
    g = g + w(i) * (1D0 / ((x-t(i))**2 + y*y) + &
                    1D0 / ((x+t(i))**2 + y*y))
  end do
  k3 = 0.318309886 * y * g
end function k3


!pure function Voigt_za(a, x)
!  ! 2007MNRAS_375_1043Zaghloul
!  double precision Voigt_za
!  double precision, intent(in) :: x, a
!  double precision v1, v2, u, du, t, t1
!  double precision, parameter :: tmp = (2D0 / sqrt(phy_Pi))
!  !
!  if (a .le. 26.6D0) then
!    v1 = exp(a*a - x*x) * erfc(a) * cos(2D0 * a * x)
!  else
!    t = a*a ! a**2
!    t1 = t*t ! a**4
!    v1 = (1D0 - 1D0 / (2D0 * t) &
!      + 3D0 / (4D0 * t1) - 15D0 / (8D0 * t*t1) &
!      + 105D0 / (16D0 * t1*t1) - 945D0 / (32D0 * t*t1*t1) &
!      + 10395D0 / (64D0 * t1*t1*t1)) &
!      / (a * sqrt(phy_Pi)) &
!      * exp(-x*x) * cos(2D0 * a * x)
!  end if
!  v2 = 0D0
!  if (x .gt. 0D0) then
!    u = 0D0
!    du = min(phy_Pi/(1D5*a), x*1D-5)
!    do
!      if (u .gt. x) then
!        exit
!      end if
!      v2 = v2 + ( &
!        exp((u-x)*(u+x)) * sin(2D0*a*(x-u)) + &
!        exp((u+du*0.5D0-x)*(u+du*0.5D0+x)) * &
!          sin(2D0*a*(x-u-du*0.5D0)) * 4D0 + &
!        exp((u+du-x)*(u+du+x)) * sin(2D0*a*(x-u-du)) &
!        )
!      u = u + du
!    end do
!  end if
!  Voigt_za = v1 + v2 * tmp * du/6D0
!end function Voigt_za



end module voigt
