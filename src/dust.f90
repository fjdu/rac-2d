module dust

use trivials
use phy_const

implicit none

integer, parameter :: const_len_name_dust_material = 16
integer, parameter :: cont_n_dust_materials = 4

type :: type_dust_optical_property_
  character(len=const_len_name_dust_material) :: material
  integer :: nlen
  double precision, dimension(:), allocatable :: &
    freq, &
    wavelen, &
    sig_abs, &
    sig_sca, &
    g_sca
end type type_dust_optical_property_


type :: type_dust_species
  character(len=const_len_name_dust_material) :: material
  double precision :: &
    radius_CGS, &
    material_density_CGS, &
    particle_mass_CGS, &
    num_fraction, &
    T, &
    n_CGS, &
    rho_CGS!, &
    !sig_abs, &
    !sig_sca, &
    !g_sca
end type type_dust_species


type :: type_dust_local_collection
  integer nItem
  type(type_dust_species), dimension(:), allocatable :: list
end type type_dust_local_collection


type :: type_dust_optical_property_p
  type(type_dust_optical_property_), pointer :: p
end type type_dust_optical_property_p


type :: type_dust_property_global_set
  integer nItem
  !type(type_dust_optical_property_), dimension(4) :: list
  type(type_dust_optical_property_p), dimension(:), allocatable :: list
end type type_dust_property_global_set


character(len=const_len_name_dust_material), dimension(cont_n_dust_materials) :: &
  dust_materials = (/ &
    'silicate        ', &
    'graphite        ', &
    'ice             ', &
    'forsterite      '/)


type(type_dust_property_global_set) :: a_dust_set

type(type_dust_local_collection) :: a_cell_dust

contains


subroutine load_dust_optical_property(filename, iItem)
  character(len=*) filename
  character(len=128) strtmp
  character, parameter :: commentchar = '!'
  integer fU, nrow, ncol, ios, strlen, i
  integer, intent(in) :: iItem
  nrow = GetFileLen_comment_blank(filename, commentchar)
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a free file unit!  In load_dust_optical_property.'
    stop
  end if
  allocate(a_dust_set%list(iItem)%p)
  a_dust_set%list(iItem)%p%nlen = nrow - 2
  allocate( &
    a_dust_set%list(iItem)%p%freq(nrow-2), &
    a_dust_set%list(iItem)%p%wavelen(nrow-2), &
    a_dust_set%list(iItem)%p%sig_abs(nrow-2), &
    a_dust_set%list(iItem)%p%sig_sca(nrow-2) &
    )
  call openFileSequentialRead(fU, filename, 99999)
  i = 0
  do
    read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
    if (ios .LT. 0) exit
    strtmp = adjustl(strtmp)
    strlen = len_trim(strtmp)
    if ((strtmp(1:1) .NE. commentchar) .AND. &
        (strlen .NE. 0)) then
      if (.not. IsDigitChar(strtmp(1:1))) then
        a_dust_set%list(iItem)%p%material = trim(strtmp)
      else if (strlen .le. 6) then
        read(strtmp(1:strlen), '(I6)') ncol
        if (ncol .gt. 3) then
          if (.not. allocated(a_dust_set%list(iItem)%p%g_sca)) then
            allocate(a_dust_set%list(iItem)%p%g_sca(nrow-2))
          end if
        end if
      else
        i = i + 1
        if (ncol .eq. 3) then
          read(strtmp, '(3F14.3)') &
                    a_dust_set%list(iItem)%p%wavelen(i), &
                    a_dust_set%list(iItem)%p%sig_abs(i), &
                    a_dust_set%list(iItem)%p%sig_sca(i)
        else if (ncol .eq. 4) then
          read(strtmp, '(4F14.3)') &
                    a_dust_set%list(iItem)%p%wavelen(i), &
                    a_dust_set%list(iItem)%p%sig_abs(i), &
                    a_dust_set%list(iItem)%p%sig_sca(i), &
                    a_dust_set%list(iItem)%p%g_sca(i)
        else
        end if
        a_dust_set%list(iItem)%p%freq(i) = &
          phy_SpeedOfLight_CGS / (a_dust_set%list(iItem)%p%wavelen(i) * 1D-4)
      end if
    end if
  end do
  close(fU)
end subroutine load_dust_optical_property


subroutine get_dust_abs_sca_crosssec(iDust, freq, abs_sig, sca_sig, sca_g)
  integer, intent(in) :: iDust
  double precision, intent(in) :: freq ! in Hz
  double precision, intent(out) :: abs_sig, sca_sig
  double precision, intent(out), optional :: sca_g
  integer idx
  idx = binary_search(a_dust_set%list(iDust)%p%freq, a_dust_set%list(iDust)%p%nlen, freq, 1)
  abs_sig = a_dust_set%list(iDust)%p%sig_abs(idx)
  sca_sig = a_dust_set%list(iDust)%p%sig_sca(idx)
  if (present(sca_g)) then
    if (allocated(a_dust_set%list(iDust)%p%g_sca)) then
      sca_g = a_dust_set%list(iDust)%p%g_sca(idx)
    else
      sca_g = 0D0
    end if
  end if
end subroutine get_dust_abs_sca_crosssec


subroutine calc_dust_particle_params_basic(d)
  type(type_dust_species), intent(inout) :: d
  d%particle_mass_CGS = 4D0*phy_Pi/3D0 * d%radius_CGS**3 * d%material_density_CGS
end subroutine calc_dust_particle_params_basic


function stellar_flux_total(L, R)
  double precision stellar_flux_total
  double precision, intent(in) :: L, R ! L in Lsun, R in AU
  stellar_flux_total = L * phy_Lsun_CGS / (4D0 * phy_Pi * (R*phy_AU2cm)**2)
end function stellar_flux_total


!function stellar_flux_blackbody_total(T, R, freq)
!  double precision stellar_flux_blackbody_total
!  double precision, intent(in) :: T, R, freq
!  stellar_flux_blackbody_total = phy_StefanBoltzmann_CGS * T**4 
!end function stellar_flux_blackbody_total


end module dust
