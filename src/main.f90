program main

use configure
use disk
!use dust
use my_timer
use montecarlo
use load_Draine_dusts
use trivials

implicit none

integer i, j

type(atimer) timer
type(date_time) a_date_time


!double precision, dimension(4) :: x = (/1D0, 2D0, 3D0, 4D0/)
!double precision, dimension(4) :: y = (/1D0, 5D0, 1D0, 1D0/)
!write(*,*) discrete_integral(4, x, y, 1.0D0, 3.0D0)

!stop


!dust_config%dir = '/Users/fdu/Dropbox/work/protoplanetary_disk/inp/'
!dust_config%ndust = 3
!dust_config%filenames(1) = 'suvSil_81'
!dust_config%weights(1) = 0.7D0
!dust_config%filenames(2) = 'Gra_81'
!dust_config%weights(2) = 0.2D0
!dust_config%filenames(3) = 'SiC_81'
!dust_config%weights(3) = 0.1D0
!call load_dusts
!call mix_dusts(dust_config%ndust, dust_data, dust_config%weights, &
!     dust_mixed)
!
!do i=1, dust_config%ndust
!  write(*,*) i, dust_data(i)%dustname_full
!  write(*,*) dust_data(i)%r(75:81)
!  write(*,*) dust_data(i)%w(235:241)
!  write(*,*) dust_data(i)%ab(235:241, 81)
!  write(*,*) dust_data(i)%sc(235:241, 81)
!  write(*,*) dust_data(i)%g(235:241, 81)
!end do
!write(*,*) 'Mixed:'
!write(*,*) dust_mixed%r(75:81)
!write(*,*) dust_mixed%w(235:241)
!write(*,*) dust_mixed%ab(235:241, 81)
!write(*,*) dust_mixed%sc(235:241, 81)
!write(*,*) dust_mixed%g(235:241, 81)
!stop
!Test

!type(type_dust_optical_property) d
!type(type_LUT_Tdust) lut
!integer n1, n2, n3, n
!type(type_stellar_spectrum) star
!type(type_montecarlo_dust) mc
!double precision lam_this, s
!integer idx
!call load_stellar_spectrum('tw_hya_spec_combined.dat', star)
!do i=star%n-20, star%n
!  write(*,*) i, star%lam(i), star%vals(i)
!end do
!mc%nph = 4000
!mc%eph = get_stellar_luminosity(star) / dble(mc%nph)
!lam_this = star%lam(1)
!idx = 1
!s = 0
!!do i=1, mc%nph
!!  write(*,*) lam_this
!!  lam_this = get_next_lam(lam_this, idx, star, mc%eph)
!!  s = s + mc%eph
!!  if (idx .ge. star%n) then
!!    exit
!!  end if
!!end do
!!write(*,*) s, get_stellar_luminosity(star)
!call load_dust_data('silicate.Kappa', d)
!write(*,*) get_idx_for_kappa(0.01D0, d)
!write(*,*) get_idx_for_kappa(0.1D0, d)
!write(*,*) get_idx_for_kappa(0.3D0, d)
!write(*,*) get_idx_for_kappa(3.3D0, d)
!write(*,*) get_idx_for_kappa(33D0, d)
!write(*,*) get_idx_for_kappa(33D2, d)
!write(*,*) get_idx_for_kappa(33D3, d)
!stop
!call make_LUT_Tdust(d, lut)
!!write(*,*) d%n, d%g(1:10)
!write(*,*) lut%n
!do i=1,20
!  write(*,*) lut%Tds(i), lut%vals(i)
!end do
!write(*,*) lut%Tds(990:1000)
!write(*,*) lut%vals(990:1000)
!write(*,*) get_Tdust_from_LUT(0.0D0, lut)
!write(*,*) get_Tdust_from_LUT(0.1D0, lut)
!write(*,*) get_Tdust_from_LUT(0.5D0, lut)
!write(*,*) get_Tdust_from_LUT(1.5D0, lut)
!write(*,*) get_Tdust_from_LUT(2.5D0, lut)
!write(*,*) get_Tdust_from_LUT(2.5D3, lut)
!write(*,*) get_Tdust_from_LUT(2.5D6, lut)
!write(*,*) get_Tdust_from_LUT(2.5D9, lut)
!write(*,*) get_Tdust_from_LUT(2.5D13, lut)
!write(*,*) get_Tdust_from_LUT(1.5D15, lut)
!write(*,*) get_Tdust_from_LUT(2.5D23, lut)
!write(*,*) get_Tdust_from_LUT(2.5D33, lut)
!
!type(type_spectrum_generic) sp
!type(type_distribution_table) di
!sp%n = 3
!allocate(sp%intervals(2, sp%n), sp%vals(sp%n))
!sp%intervals(1, 1) = 1D0
!sp%intervals(2, 1) = 2D0
!sp%intervals(1, 2) = 6D0
!sp%intervals(2, 2) = 9D0
!sp%intervals(1, 3) = 16D0
!sp%intervals(2, 3) = 19D0
!sp%vals(1) = 3D0
!sp%vals(2) = 1D0
!sp%vals(3) = 1D1
!call convert_spec_to_distri(sp, di, 0)
!write(*,*) di%pvals
!call init_random_seed
!call srand(86456)
!n = 20000000
!n1 = 0
!n2 = 0
!n3 = 0
!call timer%init('Main')
!do i=1, n
!  if (get_a_sample(di) .eq. 1) then
!    n1 = n1 + 1
!  end if
!  if (get_a_sample(di) .eq. 2) then
!    n2 = n2 + 1
!  end if
!  if (get_a_sample(di) .eq. 3) then
!    n3 = n3 + 1
!  end if
!end do
!call timer%elapse
!write(*,*) real(n1)/real(n), di%pvals(1)
!write(*,*) real(n2)/real(n), di%pvals(2) - di%pvals(1)
!write(*,*) real(n3)/real(n), di%pvals(3) - di%pvals(2)
!
!!write(*,*) planck_B_nu(5D3, 1D13)
!!write(*,*) planck_B_nu(5D3, 1D14)
!!write(*,*) planck_B_nu(5D3, 1D15)
!!write(*,*) planck_B_nu(5D3, 1D16)
!!!
!!write(*,*) planck_B_lambda(5D3, 1D-3) * 1D-7 * 1D-3 * 1D4 * 1D-7
!!write(*,*) planck_B_lambda(5D3, 1D-4) * 1D-7 * 1D-3 * 1D4 * 1D-7
!!write(*,*) planck_B_lambda(5D3, 0.5D-4) * 1D-7 * 1D-3 * 1D4 * 1D-7
!!write(*,*) planck_B_lambda(5D3, 1D-5) * 1D-7 * 1D-3 * 1D4 * 1D-7
!!write(*,*) planck_B_lambda(5D3, 1D-6) * 1D-7 * 1D-3 * 1D4 * 1D-7
!!stop
!write(*,*) Voigt(1D0, 1D-20)
!write(*,*) Voigt(10D0, 1D-4)
!write(*,*) Voigt(1D-3, 1D-3)
!write(*,*) Voigt(0D0, 0.25D0)
!write(*,*) Voigt(1D0, 0.5D0)
!write(*,*) Voigt(5D0, 5D0)
!write(*,*) Voigt(1D0, 10D0)
!call make_H_Lya(300D0, HI_0)
!stop
!End test

call get_command_argument(0, a_disk%filename_exe, i, j)
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
