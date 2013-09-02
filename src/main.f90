program main

use configure

use disk
use dust
use spline_1d_2d
use barycentric_1d_2d
use statistic_equilibrium

use my_timer

implicit none

integer i, j

type(atimer) timer


!! Test
!type(type_molecule_energy_set), target :: a_molecule_set
!double precision x, a
!!
!a_molecule_using => a_molecule_set
!call load_moldata_LAMBDA('../inp/12C16O_H2.dat')
!!
!a_molecule_using%Tkin = 20D0
!a_molecule_using%density_mol = 1D2
!a_molecule_using%dv = 1D5
!a_molecule_using%length_scale = 1D13
!a_molecule_using%f_occupation = 0D0
!a_molecule_using%f_occupation(1) = 1D0
!!a_molecule_using%f_occupation = a_molecule_using%level_list%weight * &
!!    exp(-a_molecule_using%level_list%energy / a_molecule_using%Tkin)
!!a_molecule_using%f_occupation = a_molecule_using%f_occupation / sum(a_molecule_using%f_occupation)
!do i=1, a_molecule_using%colli_data%n_partner
!  a_molecule_using%colli_data%list(i)%dens_partner = 1D1
!end do
!!
!call timer%init('Main')
!call statistic_equil_solve
!call timer%elapse
!!
!call calc_cooling_rate
!!
!write(*,*) a_molecule_using%level_list%energy
!write(*,*) a_molecule_using%f_occupation
!write(*,*) a_molecule_using%f_occupation(2)/a_molecule_using%f_occupation(1)
!write(*,*) exp(-a_molecule_using%level_list(2)%energy / a_molecule_using%Tkin) * &
!    a_molecule_using%level_list(2)%weight / a_molecule_using%level_list(1)%weight
!write(*,*) a_molecule_using%cooling_rate_total
!!
!x = precision(0D0)
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = tiny(0D0)
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-5
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-6
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-7
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-8
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-9
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-10
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-15
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-20
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-30
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-50
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-150
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!x = 1D-250
!write(*,*) x, sin(x)/x, (1D0 - exp(-x)) / x
!!
!call timer%init('Main')
!x = 1D0
!do i=1, 10000000
!  x = x + 1D-16
!end do
!call timer%elapse
!write(*,*) x
!!
!call timer%init()
!x = 0D0
!do i=1, 10000000
!  x = x + 1D0 / dble(i)
!end do
!call timer%elapse
!write(*,*) x
!!
!call timer%init()
!x = 0D0
!do i=10000000, 1, -1
!  x = x + 1D0 / dble(i)
!end do
!call timer%elapse
!write(*,*) x
!!
!call timer%init()
!do i=1, 20000000
!  x = sqrt(sqrt(dble(i)))
!end do
!call timer%elapse
!!
!call timer%init()
!do i=1, 20000000
!  x = dble(i)**0.25D0
!end do
!call timer%elapse
!!
!a = -10.73D0
!call timer%init()
!do i=1, 20000000
!  x = dble(i)**a
!end do
!call timer%elapse
!!
!call timer%init()
!do i=1, 20000000
!  x = exp(log(dble(i))*a)
!end do
!call timer%elapse
!stop
!!end test



call get_command_argument(0, disk_params_ini%filename_exe, i, j)
call get_command_argument(1, filename_config, i, j)
if (i .EQ. 0) then
  filename_config = 'configure.dat'
end if

call config_do

!call a_test_case
!stop

call timer%init('Main')

call disk_iteration

if (FileUnitOpened(a_book_keeping%fU)) then
  close(a_book_keeping%fU)
end if


call timer%elapse

end program main
