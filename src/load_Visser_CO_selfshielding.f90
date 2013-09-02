module load_Visser_CO_selfshielding

use trivials

! Visser 2009

implicit none

integer, parameter, private :: nrow = 6, ncol = 8

double precision, dimension(nrow), parameter, private :: &
  logN_H2 = (/0, 19, 20, 21, 22, 23/)

double precision, dimension(ncol), parameter, private :: &
  logN_12CO = (/0, 13, 14, 15, 16, 17, 18, 19/)

! For Tex(CO) = 5 K
!double precision, dimension(ncol, nrow), parameter, private :: &
!  f_12CO = reshape((/ &
!    1.000D00, 8.080D-1, 5.250D-1, 2.434D-1, 5.467D-2, 1.362D-2, 3.378D-3, 5.240D-4, &
!    8.176D-1, 6.347D-1, 3.891D-1, 1.787D-1, 4.297D-2, 1.152D-2, 2.922D-3, 4.662D-4, &
!    7.223D-1, 5.624D-1, 3.434D-1, 1.540D-1, 3.515D-2, 9.231D-3, 2.388D-3, 3.899D-4, &
!    3.260D-1, 2.810D-1, 1.953D-1, 8.726D-2, 1.907D-2, 4.768D-3, 1.150D-3, 1.941D-4, &
!    1.108D-2, 1.081D-2, 9.033D-3, 4.441D-3, 1.102D-3, 2.644D-4, 7.329D-5, 1.437D-5, &
!    3.938D-7, 3.938D-7, 3.936D-7, 3.923D-7, 3.901D-7, 3.893D-7, 3.890D-7, 3.875D-7 /), &
!    (/ncol, nrow/))

! For Tex(CO) = 50 K
double precision, dimension(ncol, nrow), parameter, private :: &
  f_12CO = reshape((/ &
    1.000E+0, 9.405E-1, 7.046E-1, 4.015E-1, 9.964E-2, 1.567E-2, 3.162E-3, 4.839E-4, &
    7.546E-1, 6.979E-1, 4.817E-1, 2.577E-1, 6.505E-2, 1.135E-2, 2.369E-3, 3.924E-4, &
    5.752E-1, 5.228E-1, 3.279E-1, 1.559E-1, 3.559E-2, 6.443E-3, 1.526E-3, 2.751E-4, &
    2.493E-1, 2.196E-1, 1.135E-1, 4.062E-2, 7.864E-3, 1.516E-3, 4.448E-4, 9.367E-5, &
    1.550E-3, 1.370E-3, 6.801E-4, 2.127E-4, 5.051E-5, 1.198E-5, 6.553E-6, 3.937E-6, &
    8.492E-8, 8.492E-8, 8.492E-8, 8.492E-8, 8.492E-8, 8.492E-8, 8.488E-8, 8.453E-8 /), &
    (/ncol, nrow/))


contains


function get_12CO_shielding(N_H2, N_12CO)
  double precision  get_12CO_shielding, N_H2, N_12CO
  double precision logN_H2_, logN_12CO_
  integer i, i1, j1
  logN_H2_   = log10(max(N_H2, 1D0))
  logN_12CO_ = log10(max(N_12CO, 1D0))
  if (logN_H2_ .gt. logN_H2(nrow)) then
    i1 = nrow - 1
  else if (logN_H2_ .lt. logN_H2(1)) then
    i1 = 1
  else
    do i=1, nrow-1
      if ((logN_H2(i) .le. logN_H2_) .and. (logN_H2(i+1) .ge. logN_H2_)) then
        i1 = i
        exit
      end if
    end do
  end if
  if (logN_12CO_ .gt. logN_12CO(ncol)) then
    j1 = ncol - 1
  else if (logN_12CO_ .lt. logN_12CO(1)) then
    j1 = 1
  else
    do i=1, ncol-1
      if ((logN_12CO(i) .le. logN_12CO_) .and. (logN_12CO(i+1) .ge. logN_12CO_)) then
        j1 = i
        exit
      end if
    end do
  end if
  get_12CO_shielding = calc_four_point_linear_interpol( &
    logN_12CO_,    logN_H2_, &
    logN_12CO(j1), logN_12CO(j1+1), &
    logN_H2(i1),   logN_H2(i1+1), &
    log(f_12CO(j1,i1)),   log(f_12CO(j1,i1+1)), &
    log(f_12CO(j1+1,i1)), log(f_12CO(j1+1,i1+1)))
  get_12CO_shielding = exp(get_12CO_shielding)
end function get_12CO_shielding


end module load_Visser_CO_selfshielding
