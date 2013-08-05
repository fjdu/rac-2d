module load_Visser_CO_selfshielding

use trivials

! Visser 2009

implicit none

integer, parameter, private :: nrow = 6, ncol = 8

double precision, dimension(nrow), parameter, private :: &
  logN_H2 = (/0, 19, 20, 21, 22, 23/)

double precision, dimension(ncol), parameter, private :: &
  logN_12CO = (/0, 13, 14, 15, 16, 17, 18, 19/)

double precision, dimension(ncol, nrow), parameter, private :: &
  f_12CO = reshape((/ &
    1.000D00, 8.080D-1, 5.250D-1, 2.434D-1, 5.467D-2, 1.362D-2, 3.378D-3, 5.240D-4, &
    8.176D-1, 6.347D-1, 3.891D-1, 1.787D-1, 4.297D-2, 1.152D-2, 2.922D-3, 4.662D-4, &
    7.223D-1, 5.624D-1, 3.434D-1, 1.540D-1, 3.515D-2, 9.231D-3, 2.388D-3, 3.899D-4, &
    3.260D-1, 2.810D-1, 1.953D-1, 8.726D-2, 1.907D-2, 4.768D-3, 1.150D-3, 1.941D-4, &
    1.108D-2, 1.081D-2, 9.033D-3, 4.441D-3, 1.102D-3, 2.644D-4, 7.329D-5, 1.437D-5, &
    3.938D-7, 3.938D-7, 3.936D-7, 3.923D-7, 3.901D-7, 3.893D-7, 3.890D-7, 3.875D-7 /), &
    (/ncol, nrow/))


contains


function get_12CO_shielding(N_H2, N_12CO)
  double precision  get_12CO_shielding, N_H2, N_12CO
  double precision logN_H2_, logN_12CO_
  integer i, i1, j1
  logN_H2_ = log10(N_H2)
  logN_12CO_ = log10(N_12CO)
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
    logN_12CO_, logN_H2_, &
    logN_12CO(j1), logN_12CO(j1+1), &
    logN_H2(i1), logN_H2(i1+1), &
    f_12CO(j1,i1), f_12CO(j1,i1+1), &
    f_12CO(j1+1,i1), f_12CO(j1+1,i1+1))
end function get_12CO_shielding


end module load_Visser_CO_selfshielding
