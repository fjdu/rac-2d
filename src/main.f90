program main

use configure
use disk
use my_timer
use trivials
use phy_const
use ray_tracing
use cdms

implicit none

integer i, j

type(atimer) timer
type(date_time) a_date_time

!!! TEST REGION
!!! END OF TEST REGION

call get_command_argument(0, a_disk_iter_params%filename_exe, i, j)
call get_command_argument(1, filename_config, i, j)
if (i .EQ. 0) then
  filename_config = 'configure.dat'
end if

call config_do

call timer%init('Main')

call disk_iteration

if (a_disk_iter_params%rerun_single_points) then
  call do_rerun_single_points
  stop
end if

call post_disk_iteration
!
if (a_disk_iter_params%do_continuum_transfer) then
  call continuum_tran_prep
  !
  if (len_trim(raytracing_conf%dir_save_image) .eq. 0) then
    raytracing_conf%dir_save_image = a_disk_iter_params%iter_files_dir
  end if
  write(*, '(/A)') 'Doing continuum radiative transfer...'
  call make_cubes_continuum
end if
!
if (a_disk_iter_params%do_line_transfer) then
  call line_tran_prep
  !
  raytracing_conf%VeloHalfWidth = 1.2D0 * sqrt( &
    phy_GravitationConst_SI * a_disk%star_mass_in_Msun * phy_Msun_SI / &
      (root%xmin * phy_AU2m))
  !
  if (len_trim(raytracing_conf%dir_save_image) .eq. 0) then
    raytracing_conf%dir_save_image = a_disk_iter_params%iter_files_dir
  end if
  !
  call line_excitation_do
  !
  write(*, '(/A)') 'Doing line radiative transfer...'
  call make_cubes_line
end if
!
if (FileUnitOpened(a_book_keeping%fU)) then
  write(a_book_keeping%fU, '(A)') '! Current time: ' // &
    trim(a_date_time%date_time_str())
  close(a_book_keeping%fU)
end if

write(*, '(A)') 'Current time: ' // trim(a_date_time%date_time_str())
call timer%elapse

end program main
