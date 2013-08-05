module analytical_H2O_OH

use phy_const
implicit none

private

double precision, parameter :: &
  sigma_2 = 1.8D-18, & ! sigma_OH
  sigma_4 = 1.2D-17    ! sigma_H2O

double precision, parameter :: &
  F0 = 1D9, &
  n0 = 1D12, &
  T0 = 1000D0, &
  X_O = 3D-4, &
  z_scale = 1D0, &
  z_max = 5D0

public dF_dL, ana_solve, Density, Temperature, n_H2O, n_OH


contains


function dF_dL(F, L)
  double precision dF_dL, F, L
  dF_dL = -F * kappa(L, F) * phy_AU2cm
end function dF_dL


function n_OH(L, F)
  double precision n_OH, L, F
  double precision f_modify, t_scale
  t_scale = 1D0 / (k1(L, F) * 0.5D0 * Density(L)) / phy_SecondsPerYear
  f_modify = exp(-t_scale/1D7)
  n_OH = Density(L) * X_O / &
    (1D0 + &
     k3(L, F) * 0.5D0 * Density(L) / k4(L, F) + &
     k2(L, F) / (k1(L, F) * 0.5D0 * Density(L))) * &
    f_modify
end function n_OH


function n_H2O(L, F)
  double precision n_H2O, L, F
  double precision f_modify, t_scale
  t_scale = 1D0 / (k1(L, F) * 0.5D0 * Density(L)) / phy_SecondsPerYear
  f_modify = exp(-t_scale/1D7)
  n_H2O = Density(L) * X_O / &
    (1D0 + &
     k4(L, F) / (k3(L, F) * 0.5D0 * Density(L)) * (1D0 + k2(L, F) / (k1(L, F) * 0.5D0 * Density(L)))) * &
    f_modify
end function n_H2O


function kappa(L, F)
  double precision kappa, L, F
  kappa = n_OH(L, F) * sigma_2 + n_H2O(L, F) * sigma_4
end function kappa


function Temperature(L, F)
  double precision Temperature, L, F
  Temperature = T0 * exp(-L/(z_scale*1.5D0)) + 20D0
end function Temperature


function Density(L)
  double precision Density, L
  Density = n0 * exp(-((z_max - L) / z_scale)**2)
end function Density


function k1(L, F)
  double precision k1, L, F, T
  T = Temperature(L, F)
  k1 = 3.14D-13 * (T/300D0)**2.70D0 * exp(-3150D0/T)
end function k1


function k3(L, F)
  double precision k3, L, F, T
  T = Temperature(L, F)
  k3 = 2.05D-12 * (T/300D0)**1.52D0 * exp(-1736D0/T)
end function k3


function k2(L, F)
  double precision k2, L, F
  k2 = F * sigma_2
end function k2


function k4(L, F)
  double precision k4, L, F
  k4 = F * sigma_4
end function k4


subroutine ana_solve(ndim, L_s, F_s)
  external ana_f, ana_jac
  integer ndim
  integer i
  double precision dL
  double precision, dimension(ndim) :: L_s, F_s
  double precision, dimension(1) :: y
  double precision t, tout
  double precision, dimension(128) :: RWORK
  integer, dimension(128) :: IWORK
  double precision :: &
    ATOL = 1D-30, &
    RTOL = 1D-10
  integer :: &
    NEQ = 1, &
    ITOL = 1, &
    ITASK = 1, &
    ISTATE = 1, &
    IOPT = 0, &
    LRW = 128, &
    LIW = 128, &
    MF = 22
  !
  dL = z_max/dble(ndim-1)
  L_s(1) = 0D0
  do i=2, ndim
    L_s(i) = L_s(i-1) + dL
  end do
  F_s(1) = F0
  do i=2, ndim
    y(1) = F_s(i-1)
    t = L_s(i-1)
    tout = L_s(i)
    call DLSODE( &
       ana_f, &
       !
       NEQ, &
       y, &
       !
       t, &
       tout, &
       !
       ITOL, &
       RTOL, &
       ATOL, &
       ITASK, &
       ISTATE, &
       IOPT, &
       RWORK, &
       LRW, &
       IWORK, &
       LIW, &
       !
       ana_jac, &
       !
       MF)
    if (ISTATE .LT. 0) then
      write(*,*) 'Error!'
    end if
    L_s(i) = tout
    F_s(i) = y(1)
  end do
end subroutine ana_solve


end module analytical_H2O_OH


subroutine ana_f(NEQ, t, y, ydot)
  use analytical_H2O_OH
  implicit none
  integer NEQ
  double precision t, y(NEQ), ydot(NEQ)
  ! t: L
  ! y: F
  ! ydot: dF/dL = -F * kappa
  ydot(1) = dF_dL(y(1), t)
end subroutine ana_f


subroutine ana_jac(NEQ, t, y, j, ian, jan, pdj)
  implicit none
  double precision t
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, j
end subroutine ana_jac


program main
use analytical_H2O_OH
implicit none

integer i

integer, parameter :: ndim = 128

double precision, dimension(ndim) :: L_s, F_s

double precision L, F

call ana_solve(ndim, L_s, F_s)

do i=1, ndim
  L = L_s(i)
  F = F_s(i)
  write(*,*) i, L, F, Temperature(L, F), Density(L), n_H2O(L, F), n_OH(L, F)
end do

end program main
