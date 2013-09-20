program main

use configure
use disk
use dust
use my_timer

implicit none

integer i, j

type(atimer) timer
type(date_time) a_date_time

call get_command_argument(0, disk_params_ini%filename_exe, i, j)
call get_command_argument(1, filename_config, i, j)
if (i .EQ. 0) then
  filename_config = 'configure.dat'
end if

call config_do

call timer%init('Main')

!call a_test_case
!call b_test_case
!stop


call disk_iteration

if (FileUnitOpened(a_book_keeping%fU)) then
  write(a_book_keeping%fU, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
  close(a_book_keeping%fU)
end if

call timer%elapse

end program main
