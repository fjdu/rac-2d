program main

use configure

use disk
use dust
use spline_1d_2d
use barycentric_1d_2d

use my_timer

implicit none

integer i, j

type(atimer) timer

call get_command_argument(0, disk_params_ini%filename_exe, i, j)
call get_command_argument(1, filename_config, i, j)
if (i .EQ. 0) then
  filename_config = 'configure.dat'
end if

call config_do

call timer%init('Main')

call disk_iteration

if (a_disk_iter_params%flag_converged) then
  write(*, '(A/)') "Iteration has converged!"
else
  write(*, '(A/)') "Iteration hasn't converged. :("
end if

if (FileUnitOpened(a_book_keeping%fU)) then
  close(a_book_keeping%fU)
end if


call timer%elapse

end program main
