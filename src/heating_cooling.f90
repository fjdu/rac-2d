module heating_cooling

use phy_const
use data_struct
use trivials
use lamda
use statistic_equilibrium
use chemistry
use grid
use cont_lookuptable
use binary_array_io
use spline_1d_2d

implicit none

type :: type_heating_cooling_config
  logical :: use_analytical_CII_OI = .true.
  logical :: use_mygasgraincooling = .true.
  logical :: use_chemicalheatingcooling = .true.
  logical :: use_Xray_heating = .true.
  logical :: use_phdheating_H2 = .true.
  logical :: use_phdheating_H2OOH = .true.
  logical :: dust_gas_linear_couple = .false.
  integer :: solve_method = 1 ! 1: ODE; 2: Newton
  double precision :: heating_eff_chem = 1D0
  double precision :: heating_eff_H2form = 0.1D0
  double precision :: heating_eff_phd_H2 = 1D0
  double precision :: heating_eff_phd_H2O = 0.1D0
  double precision :: heating_eff_phd_OH = 0.1D0
  double precision :: heating_Xray_en = 18.7D0
  double precision :: cooling_gg_coeff = 0.3D0
  character(len=128) :: dir_transition_rates = './transitions/'
  character(len=128) :: filename_CII = 'C+.dat'
  character(len=128) :: filename_NII = 'N+.dat'
  character(len=128) :: filename_OI = 'Oatom.dat'
  character(len=128) :: filename_FeII = 'Fe+.dat'
  character(len=128) :: filename_SiII = 'Si+.dat'
  logical :: IonCoolingWithLut = .true.
end type type_heating_cooling_config

double precision, public :: hc_Tgas, hc_Tdust

type(type_cell_rz_phy_basic), pointer :: hc_params => null()

type(type_heating_cooling_rates_list) :: heating_cooling_rates

type(type_molecule_energy_set), target :: &
  molecule_CII, molecule_NII, molecule_OI, molecule_FeII, molecule_SiII

type(type_spline_2D) :: spl_NII, spl_SiII, spl_FeII

type(type_heating_cooling_config) :: heating_cooling_config

namelist /heating_cooling_configure/ &
  heating_cooling_config

double precision, private, parameter :: very_small_num = 1D-100
double precision, public, parameter :: frac_dust_lose_en = 0.8D0

! To be used by the DLSODE solver.
! Defined here so that this common block can be preserved.
double precision RLS(218)
integer ILS(37)
COMMON /DLS001/ RLS,ILS

contains


subroutine heating_cooling_prepare
  ! Load the molecular data
  integer n_level_max
  if (.not. heating_cooling_config%use_analytical_CII_OI) then
    call load_a_mol_data(heating_cooling_config%filename_CII, molecule_CII)
    call load_a_mol_data(heating_cooling_config%filename_OI, molecule_OI)
  end if
  !
  if (heating_cooling_config%IonCoolingWithLut) then
    call create_spline2d_from_file( &
      combine_dir_filename(heating_cooling_config%dir_transition_rates, &
                           heating_cooling_config%filename_NII), spl_NII)
    call create_spline2d_from_file( &
      combine_dir_filename(heating_cooling_config%dir_transition_rates, &
                           heating_cooling_config%filename_FeII), spl_FeII)
    call create_spline2d_from_file( &
      combine_dir_filename(heating_cooling_config%dir_transition_rates, &
                           heating_cooling_config%filename_SiII), spl_SiII)
  else
    call load_a_mol_data(heating_cooling_config%filename_NII, molecule_NII)
    call load_a_mol_data(heating_cooling_config%filename_FeII, molecule_FeII)
    call load_a_mol_data(heating_cooling_config%filename_SiII, molecule_SiII)
  end if
  !
  if ((.not. heating_cooling_config%use_analytical_CII_OI) .or. &
      (.not. heating_cooling_config%IonCoolingWithLut)) then
    ! Initialize the storage to be used by the solver
    n_level_max = max(molecule_CII%n_level, molecule_NII%n_level, &
      molecule_OI%n_level, molecule_FeII%n_level, molecule_SiII%n_level)
    !
    call init_statistic_sol(n_level_max, heating_cooling_config%solve_method)
  end if
  !
end subroutine heating_cooling_prepare


subroutine prepare_cont_lut_for_line_cooling
  type(type_cell), pointer :: c
  integer i
  write(*, '(A/)') 'Preparing continuum lookup table for line cooling.'
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    call allocate_local_cont_lut(c)
    call make_local_cont_lut(c)
  end do
end subroutine prepare_cont_lut_for_line_cooling


subroutine load_a_mol_data(fname, mol)
  character(len=*), intent(in) :: fname
  type(type_molecule_energy_set), target, intent(inout) :: mol
  integer i
  write(*, '(A, A)') 'Working with file: ', combine_dir_filename( &
      heating_cooling_config%dir_transition_rates, fname)
  if (len_trim(fname) .gt. 0) then
    mol_sta_sol => mol
    call load_moldata_LAMDA(combine_dir_filename( &
      heating_cooling_config%dir_transition_rates, fname), mol_sta_sol)
    nullify(mol_sta_sol)
    !
    write(*, '(A, A)') mol%name_molecule, ' is loaded.'
    write(*, '(2X, I8, A)') mol%n_level, ' energy levels.'
    write(*, '(2X, I8, A)') mol%rad_data%n_transition, ' radiative transitions.'
    write(*, '(2X, I8, A)') mol%colli_data%n_partner, ' collisional partners.'
    do i=1, mol%colli_data%n_partner
      write(*, '(2X, A, I3, A, A)') 'Partner ', i, ': ',  mol%colli_data%list(i)%name_partner
      write(*, '(12X, I8, A)') mol%colli_data%list(i)%n_transition, ' collisional transitions.'
    end do
  else
    mol%n_level = 0
  end if
end subroutine load_a_mol_data


subroutine heating_cooling_prepare_molecule
  integer i
  mol_sta_sol%Tkin = hc_Tgas
  mol_sta_sol%dv = hc_params%velo_width_turb
  mol_sta_sol%length_scale = &
    min(hc_params%coherent_length, hc_params%Ncol_toISM / hc_params%n_gas)
  mol_sta_sol%f_occupation = mol_sta_sol%level_list%weight * &
      exp(-mol_sta_sol%level_list%energy / mol_sta_sol%Tkin)
  mol_sta_sol%f_occupation = mol_sta_sol%f_occupation &
                             / sum(mol_sta_sol%f_occupation)
  do i=1, mol_sta_sol%colli_data%n_partner
    if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'H2') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        hc_params%n_gas * hc_params%X_H2
    else if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'o-H2') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        0.75D0 * hc_params%n_gas * hc_params%X_H2
    else if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'p-H2') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        0.25D0 * hc_params%n_gas * hc_params%X_H2
    else if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'H') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        hc_params%n_gas * hc_params%X_HI
    else if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'H+') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        hc_params%n_gas * hc_params%X_Hplus
    else if (mol_sta_sol%colli_data%list(i)%name_partner .eq. 'e') then
      mol_sta_sol%colli_data%list(i)%dens_partner = &
        hc_params%n_gas * hc_params%X_E
    else
      mol_sta_sol%colli_data%list(i)%dens_partner = 0D0
    end if
  end do
end subroutine heating_cooling_prepare_molecule


function heating_chemical()
  integer i, i0
  double precision heating_chemical, tmp
  if ((.not. heating_cooling_config%use_chemicalheatingcooling) .or. &
      (hc_Tgas .le. 0D0)) then
    heating_chemical = 0D0
    return
  end if
  tmp = chem_params%Tgas
  if (chem_params%Tgas .ne. hc_Tgas) then
    chem_params%Tgas = hc_Tgas
    call chem_cal_rates
  end if
  heating_chemical = 0D0
  do i=1, chem_net%nReacWithHeat
    i0 = chem_net%iReacWithHeat(i)
    heating_chemical = heating_chemical + &
      chem_net%rates(i0) * &
      chemsol_stor%y(chem_net%reac(1, i0)) * &
      chemsol_stor%y(chem_net%reac(2, i0)) * &
      chem_net%heat(i)
  end do
  heating_chemical = heating_chemical * hc_params%n_gas / phy_SecondsPerYear * &
    heating_cooling_config%heating_eff_chem
  chem_params%Tgas = tmp
end function heating_chemical


subroutine heating_chemical_termbyterm(Tgas, nr, vecr)
  integer i, i0
  double precision, intent(in) :: Tgas
  integer, intent(in) :: nr
  double precision, intent(out), dimension(nr) :: vecr
  double precision tmp
  tmp = chem_params%Tgas
  chem_params%Tgas = Tgas
  call chem_cal_rates
  vecr = 0D0
  do i=1, chem_net%nReacWithHeat
    i0 = chem_net%iReacWithHeat(i)
    vecr(i0) = &
      chem_net%rates(i0) * &
      chemsol_stor%y(chem_net%reac(1, i0)) * &
      chemsol_stor%y(chem_net%reac(2, i0)) * &
      chem_net%heat(i) * heating_cooling_config%heating_eff_chem
  end do
  chem_params%Tgas = tmp
end subroutine heating_chemical_termbyterm



function heating_photoelectric_small_grain()
  ! Wolfire 1995, which is actually taken from
  ! Bakes 1994, equation 42 and 43.
  ! It is for very small graphitic grains and PAHs.
  ! Applicable for T up to 1D4 K.
  ! Note that in Bakes 1994 they made important changes from their
  ! section 4 to section 5.
  !
  ! Output unit = erg s-1 cm-3
  !
  double precision heating_photoelectric_small_grain
  double precision chi, n_e, t1, t2, tmp
  if ((hc_params%X_E .le. 0D0) .or. (hc_Tgas .le. 0D0)) then
    heating_photoelectric_small_grain = 0D0
    return
  end if
  chi = hc_params%G0_UV_toISM * &
            exp(-phy_UVext2Av * hc_params%Av_toISM) + &
        hc_params%G0_UV_toStar * &
            exp(-phy_UVext2Av * hc_params%Av_toStar)
  n_e = hc_params%X_E * hc_params%n_gas
  tmp = chi*sqrt(hc_Tgas) / (n_e + very_small_num)
  if ((tmp .le. 0D0) .or. isnan(tmp)) then
    t1 = 0D0
  else
    t1 = exp(0.73D0 * log(tmp))
  end if
  if (hc_Tgas .le. 0D0) then
    t2 = 0D0
  else
    t2 = exp(0.70D0 * log(1D-4 * hc_Tgas))
  end if
  heating_photoelectric_small_grain = &
    1D-24 * chi * hc_params%n_gas * &
    hc_params%PAH_abundance/const_PAH_abundance_0 * &
    ( &
    4.87D-2 / (1D0 + 4D-3 * t1) &
    + &
    3.65D-2 * t2 / (1D0 + 2D-4*tmp))
end function heating_photoelectric_small_grain


!function heating_photoelectric_small_grain()
!  ! Weingartner 2001, equation 44, table 2
!  double precision heating_photoelectric_small_grain
!  double precision, parameter :: & ! Rv = 5.5, bc = 3.0D-5
!    C0 = 1.84D0, &
!    C1 = 3.81D0, &
!    C2 = 0.08348D0, &
!    C3 = 0.00391D0, &
!    C4 = 0.089D0, &
!    C5 = 0.328D0, &
!    C6 = 0.778D0
!  associate( &
!    chi => hc_params%UV_G0_factor * exp(-phy_UVext2Av*hc_params%Av), &
!    n_gas   => hc_params%n_gas, &
!    n_e   => hc_params%X_E * hc_params%n_gas, &
!    Tgas  => hc_Tgas)
!    associate(tmp => chi*sqrt(Tgas)/(n_e))
!      heating_photoelectric_small_grain = &
!        1D-26 * chi * n_gas * hc_params%dust_depletion * &
!        (C0 + C1 * Tgas**C4) / &
!        (1D0 + C2 * tmp**C5 * (1D0 + C3 * tmp**C6))
!    end associate
!  end associate
!end function heating_photoelectric_small_grain


function heating_formation_H2()
  ! Rollig 2006
  ! Sternberg 1989, equation D5
  ! 2.4D-12 erg = 1/3 * 4.5 eV
  !
  ! 2014-06-27 Fri 02:21:05
  ! See http://adsabs.harvard.edu/abs/1978ApJ...226..477H
  double precision heating_formation_H2
  double precision, parameter :: energyPerEvent = 2.4D-12
  heating_formation_H2 = &
    energyPerEvent * hc_params%R_H2_form_rate * &
    heating_cooling_config%heating_eff_H2form
end function heating_formation_H2


function heating_cosmic_ray_Woitke()
  ! Woitke 2009, equation 100
  ! Would not be very different from Bruderer 2009
  double precision heating_cosmic_ray_Woitke
  associate( &
        n_HI             => hc_params%n_gas * hc_params%X_HI, &
        n_H2             => hc_params%n_gas * hc_params%X_H2, &
        zeta_cosmicray_H2 => hc_params%zeta_cosmicray_H2)
    heating_cosmic_ray_Woitke = &
      zeta_cosmicray_H2 * (n_HI * 5.5D-12 + n_H2 * 2.5D-11) &
      * exp(-hc_params%Ncol_toISM / const_cosmicray_attenuate_N)
  end associate
end function heating_cosmic_ray_Woitke


function heating_cosmic_ray()
  ! Bruderer 2009
  !Return unit = erg s-1 cm-3
  ! 1.5D-11 = 9 * 1.6D-19 / 1D-7
  double precision heating_cosmic_ray
  heating_cosmic_ray = &
    1.5D-11 * hc_params%zeta_cosmicray_H2 * hc_params%n_gas &
    * exp(-hc_params%Ncol_toISM / const_cosmicray_attenuate_N)
end function heating_cosmic_ray


function heating_vibrational_H2()
  ! Rollig 2006, equation C.2 - C.3
  ! It causes problems: T cannot be solved if it is included.
  ! 9.4D-22 =  2.9e-10*23500.0*1.38e-16
  ! It seems this function duplicates another one.
  !   Yes, the other one is redundant.  That one is a
  !   special case of this one; there the f factor in
  !   Rollig 2006 is taken to be 2.8D-3.
  double precision heating_vibrational_H2
  if (hc_Tgas .le. 0D0) then
    heating_vibrational_H2 = 0D0
    return
  end if
  associate( &
        n_gas   => hc_params%n_gas, &
        X_H2    => hc_params%X_H2, &
        chi => hc_params%G0_UV_toISM * &
                 exp(-phy_UVext2Av*hc_params%Av_toISM) * &
                 hc_params%f_selfshielding_toISM_H2 + &
               hc_params%G0_UV_H2phd * &
                 hc_params%f_selfshielding_toStar_H2, &
        gamma_10 => 5.4D-13 * sqrt(hc_Tgas))
    heating_vibrational_H2 = &
      (n_gas * X_H2) * chi * 9.4D-22 / &
      (1D0 + (1.9D-6 + chi*4.7D-10) / (n_gas * gamma_10))
  end associate
end function heating_vibrational_H2


function heating_photodissociation_H2()
  ! Tielens 2005, P72, equation 3.18, 3.19
  double precision heating_photodissociation_H2
  if (.not. heating_cooling_config%use_phdheating_H2) then
    heating_photodissociation_H2 = 0D0
    return
  end if
  associate( &
        chi => hc_params%G0_UV_toISM * &
                 exp(-phy_UVext2Av*hc_params%Av_toISM) * &
                 hc_params%f_selfshielding_toISM_H2 + &
               hc_params%G0_UV_H2phd * &
                 hc_params%f_selfshielding_toStar_H2, &
        n_gas => hc_params%n_gas, &
        X_H2  => hc_params%X_H2)
    heating_photodissociation_H2 = &
      4D-14 * (n_gas*X_H2) * 3.4D-10 * chi * &
      heating_cooling_config%heating_eff_phd_H2
end associate
end function heating_photodissociation_H2


function heating_photodissociation_H2O()
  ! LyAlpha = 1.63e-11 erg = 10.2 eV
  ! H(H2O) - H(OH) - H(H) = 498.826e3 J mol-1 = 8.282e-12 erg = 5.18 eV.
  double precision heating_photodissociation_H2O
  if (.not. heating_cooling_config%use_phdheating_H2OOH) then
    heating_photodissociation_H2O = 0D0
    return
  end if
  associate( &
        ! 2014-06-19 Thu 00:12:57 ! self-shielding factor added
        chi => hc_params%phflux_Lya*hc_params%f_selfshielding_toStar_H2O, &
        LyAlpha_cross_H2O => const_LyAlpha_cross_H2O, &
        ph_disso_en_H2O   => 8.07D-12 * &
            heating_cooling_config%heating_eff_phd_H2O, &
        n_H2O => hc_params%n_gas * hc_params%X_H2O)
    heating_photodissociation_H2O = &
      ph_disso_en_H2O * n_H2O * LyAlpha_cross_H2O * chi
  end associate
end function heating_photodissociation_H2O


function heating_photodissociation_OH()
  ! LyAlpha = 1.63e-11 erg = 10.2 eV
  ! H(OH) - H(O) - H(H) = 428.188 J mol-1 = 7.11e-12 erg = 4.44 eV.
  double precision heating_photodissociation_OH
  if (.not. heating_cooling_config%use_phdheating_H2OOH) then
    heating_photodissociation_OH = 0D0
    return
  end if
  associate( &
        ! 2014-06-19 Thu 00:16:18 ! self-shielding factor added
        chi => hc_params%phflux_Lya * hc_params%f_selfshielding_toStar_OH, &
        LyAlpha_cross_OH => const_LyAlpha_cross_OH, &
        ph_disso_en_OH   => 9.19D-12 * &
            heating_cooling_config%heating_eff_phd_OH, &
        n_OH => hc_params%n_gas * hc_params%X_OH)
    heating_photodissociation_OH = &
      ph_disso_en_OH * n_OH * LyAlpha_cross_OH * chi
  end associate
end function heating_photodissociation_OH


function heating_ionization_CI()
  ! Tielens 2005, P66, equation 3.8
  double precision heating_ionization_CI
  associate( &
        chi => hc_params%G0_UV_toISM * &
                 exp(-phy_UVext2Av * hc_params%Av_toISM) + &
               hc_params%G0_UV_toStar * &
                 exp(-phy_UVext2Av * hc_params%Av_toStar), &
        n_gas => hc_params%n_gas, &
        X_CI  => hc_params%X_CI)
    heating_ionization_CI = &
      2.2D-22 * X_CI * n_gas * chi
  end associate
end function heating_ionization_CI


function heating_Xray_Bethell()
  use load_Bethell_Xray_cross
  double precision heating_Xray_Bethell
  double precision en_heating_per_pair, gam1, gam2, R, n_gas_crit
  double precision tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, Q_el_rot, Q_diss, Q_dirvib, Q_BCvib, Q_vib
  double precision eta_H_e, eta_H2_e, eps1, eps2, epsB, epsC
  !double precision, parameter :: en_X = 1D0! keV
  !double precision, parameter :: en_deposit = 18D0 * phy_eV2erg ! 18 eV; AGN paper
  !double precision sigma
  !sigma = sigma_Xray_Bethell(en_X, &
  !  hc_params%dust_depletion, &
  !  hc_params%ratioDust2HnucNum, &
  !  hc_params%GrainRadius_CGS)
  !sigma = hc_params%sigma_Xray
  !heating_Xray_Bethell = sigma * hc_params%n_gas * en_deposit * &
  !  hc_params%Xray_flux_0 * exp(-sigma*hc_params%Ncol_toStar) !ISM)
  !
  ! Glassgold 2012, table 4
  if (.not. heating_cooling_config%use_Xray_heating) then
    heating_Xray_Bethell = 0D0
    return
  end if
  !heating_Xray_Bethell = hc_params%zeta_Xray_H2 * hc_params%n_gas * &
  !  heating_cooling_config%heating_Xray_en * phy_eV2erg
  !
  ! Glassgold 2012; Tielens 1985 A13 and A14
  if (hc_Tgas .gt. 0D0) then
    gam1 = 1D-12 * sqrt(hc_Tgas) * exp(-1000D0/hc_Tgas)
    gam2 = 1.4D-12 * sqrt(hc_Tgas) * exp(-18100D0/(hc_Tgas + 1200D0))
  else
    gam1 = 0D0
    gam2 = 0D0
  end if
  R = 2D-7
  !
  tmp1 = hc_params%X_H2 / (hc_params%X_H2 + hc_params%X_HI)
  !
  ! Glassgold 2012, eq 7-10
  if (hc_params%X_E .gt. 0D0) then
    tmp2 = 7.95D0 * hc_params%X_E**0.678D0
    tmp3 = 2.17D0 * hc_params%X_E**0.366D0
    tmp4 = 22D0 * hc_params%X_E**0.574D0
    tmp5 = 23500D0 * hc_params%X_E**0.955D0
    tmp6 = 10700D0 * hc_params%X_E**0.907D0
    tmp7 = 7.09D0 * hc_params%X_E**0.779D0
    tmp8 = 6.88D0 * hc_params%X_E**0.802D0
  else
    tmp2 = 0D0
    tmp3 = 0D0
    tmp4 = 0D0
    tmp5 = 0D0
    tmp6 = 0D0
    tmp7 = 0D0
    tmp8 = 0D0
  end if
  eta_H_e  = 1D0 - (1D0 - 0.117D0) / (1D0 + tmp2)
  eta_H2_e = 1D0 - (1D0 - 0.055D0) / (1D0 + tmp3)
  Q_el_rot = 37D0 * (hc_params%X_HI * eta_H_e + hc_params%X_H2 * eta_H2_e) &
         / (hc_params%X_HI + hc_params%X_H2)
  !
  ! ibid, eq 11
  Q_diss = 2.14D0 * tmp1 / (1D0 + tmp4)
  !
  ! ibid eq 12,13
  eps1 = 7.81D0 * (1D0 + tmp5)
  eps2 = 109D0  * (1D0 + tmp6)
  Q_dirvib = 19D0 * tmp1 * (1D0/eps1 + 2D0/eps2)
  !
  ! ibid eq 14,15
  epsB = 117D0 * (1D0 + tmp7)
  epsC = 132D0 * (1D0 + tmp8)
  Q_BCvib = 147D0 * tmp1 * (1D0/epsB + 1D0/epsC)
  !
  if ((gam1 + gam2) .gt. 0D0) then
    n_gas_crit = R / (gam1 * hc_params%X_HI + gam2 * hc_params%X_H2)
    Q_vib = hc_params%n_gas / (hc_params%n_gas + n_gas_crit) * (Q_dirvib + Q_BCvib)
  else
    Q_vib = 0D0
  end if
  !
  en_heating_per_pair = Q_el_rot + Q_diss + Q_vib
  !
  heating_Xray_Bethell = hc_params%zeta_Xray_H2 * hc_params%n_gas &
    * phy_eV2erg &
    * en_heating_per_pair
end function heating_Xray_Bethell


function heating_viscosity()
  ! From the AGN paper.
  double precision heating_viscosity
  double precision f_cutoff
  double precision :: heating_viscosity_stop_T = 2D4
  if (hc_Tgas .le. 0D0) then
    heating_viscosity = 0D0
    return
  end if
  associate( &
        rho => hc_params%n_gas * phy_mProton_CGS * &
               hc_params%MeanMolWeight, &
        c2  => phy_kBoltzmann_CGS * hc_Tgas &
               / (phy_mProton_CGS * hc_params%MeanMolWeight))
    f_cutoff = max(1D0 - hc_Tgas/heating_viscosity_stop_T, 0D0)
    heating_viscosity = 2.25D0 * hc_params%alpha_viscosity * rho * c2 &
      * hc_params%omega_Kepler * f_cutoff
  end associate
end function heating_viscosity


function cooling_photoelectric_small_grain()
  ! Wolfire 1995, which is actually taken from
  ! Bakes 1994, equation 44.
  !
  ! Note that the coefficient in Wolfire 1995 is 4.65D-30.
  ! Note sure why.
  !
  ! Output unit = erg s-1 cm-3
  double precision cooling_photoelectric_small_grain
  double precision chi, n_e, t0, t1, t2, t3, tmp
  if ((hc_params%X_E .le. 0D0) .or. (hc_Tgas .le. 0D0) .or. (hc_params%PAH_abundance .le. 0D0)) then
    cooling_photoelectric_small_grain = 0D0
    return
  end if
  chi = hc_params%G0_UV_toISM  * exp(-phy_UVext2Av * hc_params%Av_toISM) + &
        hc_params%G0_UV_toStar * exp(-phy_UVext2Av * hc_params%Av_toStar)
  n_e = hc_params%X_E * hc_params%n_gas
  tmp = chi * sqrt(hc_Tgas) / (n_e + very_small_num)
  if (tmp .le. 0D0) then
    cooling_photoelectric_small_grain = 0D0
    return
  end if
  t0 = log(hc_Tgas)
  t1 = exp(0.944D0 * t0)
  t2 = 0.735D0 * exp(-0.068D0 * t0)
  t3 = exp(t2 * log(tmp))
  cooling_photoelectric_small_grain = &
    hc_params%PAH_abundance / const_PAH_abundance_0 * &
    3.49D-30 * t1 * t3 * n_e * hc_params%n_gas
end function cooling_photoelectric_small_grain


!function cooling_photoelectric_small_grain()
!  ! Weingartner 2001, equation 45, table 3
!  double precision cooling_photoelectric_small_grain
!  double precision, parameter :: &
!    D0 = 0.4440D0, &
!    D1 = 2.067D0, &
!    D2 = -7.806D0, &
!    D3 = 1.687D0, &
!    D4 = 0.06251D0
!  associate( &
!    chi => hc_params%UV_G0_factor * exp(-phy_UVext2Av*hc_params%Av), &
!    n_gas => hc_params%n_gas, &
!    n_e   => hc_params%X_E * hc_params%n_gas, &
!    Tgas  => hc_Tgas)
!    associate(x => log(chi*sqrt(Tgas)/n_e))
!      cooling_photoelectric_small_grain = &
!        hc_params%dust_depletion * &
!        1D-28 * n_e * n_gas * Tgas**(D0 + D1/x) * &
!        exp(D2 + D3 * x - D4 * x * x)
!    end associate
!  end associate
!end function cooling_photoelectric_small_grain


function cooling_LymanAlpha()
  ! Tielens 2005, page 53, equation 2.62
  double precision cooling_LymanAlpha
  if (hc_Tgas .le. 0D0) then
    cooling_LymanAlpha = 0D0
    return
  end if
  cooling_LymanAlpha = 7.3D-19 * hc_params%n_gas**2 * &
    hc_params%X_HI * hc_params%X_E * exp(-118400D0 / hc_Tgas)
end function cooling_LymanAlpha


function cooling_free_bound()
  ! Draine 2011, page 139, equation 14.5 (optically thin case); and page 320, equation 27.22, 27.23
  double precision cooling_free_bound
  double precision t1, t2
  if (hc_Tgas .le. 0D0) then
    cooling_free_bound = 0D0
    return
  end if
  associate( &
        T   => hc_Tgas, &
        n_p => hc_params%n_gas * hc_params%X_Hplus, &
        n_E => hc_params%n_gas * hc_params%X_E, &
        T4  => hc_Tgas/1D4, &
        ZZ  => 1D0)
    t1 = log(T4 / ZZ)
    t2 = exp(t1 * (-0.7131D0-0.0115D0*t1))
    cooling_free_bound = n_E * n_p * &
      4.13D-13 * ZZ * t2 * & ! alpha_A
      (0.787D0 - 0.0230D0*t1) * phy_kBoltzmann_CGS * T ! Err_A
  end associate
end function cooling_free_bound


function cooling_free_free()
  ! Tielens 2005, page 53, equation 2.65
  ! Essentially a simpler version of Draine equation (10.12).
  ! Osterbrock 2006, p49
  double precision cooling_free_free
  if (hc_Tgas .le. 0D0) then
    cooling_free_free = 0D0
    return
  end if
  associate( &
        T   => hc_Tgas, &
        n_p => hc_params%n_gas * (hc_params%X_Hplus + hc_params%X_Heplus), &
        n_E => hc_params%n_gas * hc_params%X_E, &
        ZZ  => 1D0)
    cooling_free_free = &
      1.4D-27 * ZZ * sqrt(T) * 1.3D0 * n_E * n_p
  end associate
end function cooling_free_free


function cooling_vibrational_H2()
  ! Rollig 2006, equation C.1
  ! 8.26D-13 = 5988 * 1.38D-16
  double precision cooling_vibrational_H2
  if (hc_Tgas .le. 0D0) then
    cooling_vibrational_H2 = 0D0
    return
  end if
  associate( &
        n_gas   => hc_params%n_gas, &
        Tgas    => hc_Tgas, &
        X_H2    => hc_params%X_H2, &
        chi => hc_params%G0_UV_toISM * &
                 exp(-phy_UVext2Av*hc_params%Av_toISM) * &
                 hc_params%f_selfshielding_toISM_H2 + &
               hc_params%G0_UV_H2phd * &
                 hc_params%f_selfshielding_toStar_H2, &
        gamma_10 => 5.4D-13 * sqrt(hc_Tgas), &
        A_10 => 8.6D-7, &
        D_1 => 2.6D-11)
    cooling_vibrational_H2 = &
      8.26D-13 * gamma_10 * exp(-5988D0/Tgas) * &
        (n_gas*n_gas*X_H2) * &
        (A_10 + chi * D_1) / (gamma_10*n_gas + A_10 + chi*D_1)
  end associate
end function cooling_vibrational_H2


function cooling_Neufeld_H2_rot()
  use load_Neufeld_cooling_H2
  double precision cooling_Neufeld_H2_rot
  double precision t1
  if ((hc_Tgas .le. 0D0) .or. (hc_params%X_H2 .le. 0D0)) then
    cooling_Neufeld_H2_rot = 0D0
    return
  end if
  !
  a_Neufeld_cooling_H2_params%T = hc_Tgas
  call get_H2_rot_cool_params
  associate( &
      n_H2  => hc_params%n_gas * hc_params%X_H2, &
      L0    => a_Neufeld_cooling_H2_params%L0, &
      L_LTE => a_Neufeld_cooling_H2_params%L_LTE, &
      n_12  => a_Neufeld_cooling_H2_params%n_12, &
      alpha => a_Neufeld_cooling_H2_params%alpha)
    L0 = L0 + very_small_num
    L_LTE = L_LTE + very_small_num
    if (alpha .GT. 0D0) then
      t1 = exp(alpha * log(n_H2/n_12))
      cooling_Neufeld_H2_rot = n_H2 * n_H2 / &
        (1D0/L0 + n_H2/L_LTE + &
         1D0/L0 * t1 * (1D0 - n_12*L0/L_LTE))
    else
      cooling_Neufeld_H2_rot = n_H2 * n_H2 / (1D0/L0 + n_H2/L_LTE)
    end if
  end associate
end function cooling_Neufeld_H2_rot


function cooling_gas_grain_collision()
  ! Rollig 2006, equation A.7
  ! Z should be the abundance of grains relative to a certain value
  double precision cooling_gas_grain_collision
  double precision tmp, f_a, cs_H, cs_H2, coeff, Tg_minus_Td
  integer i
  !
  cooling_gas_grain_collision = 0D0
  !
  if (hc_Tgas .le. 0D0) then
    return
  end if
  !
  if (.not. heating_cooling_config%use_mygasgraincooling) then
    associate( &
          n_gas => hc_params%n_gas, &
          Tgas  => hc_Tgas, &
          Tdust => hc_Tdust, &
          r_g   => hc_params%GrainRadius_CGS, &
          Z     => hc_params%dust_depletion) ! Todo
      !cooling_gas_grain_collision = &
      !  3.5D-34 * sqrt(Tgas) * (Tgas - Tdust) * n_gas * n_gas * Z
      ! Hollenbach 1989, eq 2.15
      !cooling_gas_grain_collision = &
      !  1.2D-31 * n_gas * n_gas * sqrt(Tgas/1D3 * 1D-6/r_g) * &
      !  (1D0 - 0.8D0*exp(-75D0/Tgas)) * (Tgas - Tdust) * Z
      ! The AGN paper, equation (15)
      cooling_gas_grain_collision = &
        4.76D-33 * (1D0 - 0.8D0*exp(-75D0/Tgas)) * n_gas * n_gas * sqrt(Tgas) * &
        (Tgas - Tdust) * Z * (0.05D-4 / r_g)
    end associate
  else
    ! My own formula
    !
    !f_a = 1D0 - 0.8D0*exp(-75D0/hc_Tgas)
    f_a = heating_cooling_config%cooling_gg_coeff
    cs_H = sqrt((8D0/phy_Pi*phy_kBoltzmann_CGS/phy_mProton_CGS) * hc_Tgas)
    cs_H2 = cs_H / sqrt(2D0)
    tmp = 2D0 * phy_kBoltzmann_CGS * f_a * hc_params%n_gas * &
          (cs_H * (hc_params%X_HI + hc_params%X_Hplus) + cs_H2 * hc_params%X_H2)
    !
    do i=1, hc_params%ndustcompo
      coeff = tmp * &
              hc_params%sig_dusts(i) * &
              hc_params%n_dusts(i)
      if (heating_cooling_config%dust_gas_linear_couple) then
        Tg_minus_Td = (hc_Tgas - hc_params%Tdusts(i)) *  hc_params%dEmit_dTd(i) &
                                                      / (hc_params%dEmit_dTd(i) + coeff)
      else
        Tg_minus_Td =  hc_Tgas - hc_params%Tdusts(i)
      end if
      !
      hc_params%en_exchange_per_vol(i) = &
        max(coeff * Tg_minus_Td, &
            -frac_dust_lose_en * hc_params%en_gains(i) / hc_params%volume)
      !
      cooling_gas_grain_collision = &
        cooling_gas_grain_collision + hc_params%en_exchange_per_vol(i)
    end do
  end if
end function cooling_gas_grain_collision


function cooling_OI()
  double precision cooling_OI
  if (heating_cooling_config%use_analytical_CII_OI) then
    cooling_OI = cooling_OI_analytical()
  else
    cooling_OI = cooling_OI_my()
  end if
end function cooling_OI


function cooling_CII()
  double precision cooling_CII
  if (heating_cooling_config%use_analytical_CII_OI) then
    cooling_CII = cooling_CII_analytical()
  else
    cooling_CII = cooling_CII_my()
  end if
end function cooling_CII


function cooling_OI_my() result(val)
  double precision val
  val = calc_line_cooling_rate(molecule_OI, hc_params%X_OI)
end function cooling_OI_my


function cooling_CII_my() result(val)
  double precision val
  val = calc_line_cooling_rate(molecule_CII, hc_params%X_CII)
end function cooling_CII_my


function cooling_NII() result(val)
  double precision val
  double precision, parameter :: min_X_NII_for_cooling = 1D-15
  if ((hc_params%X_NII .le. min_X_NII_for_cooling) .or. &
      (hc_params%X_E .le. 0D0) .or. &
      (hc_params%n_gas .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    val = 0D0
    return
  end if
  if (heating_cooling_config%IonCoolingWithLut) then
    val = hc_params%X_NII * hc_params%n_gas * &
          exp(phy_ln10 * &  ! = 10**()
          spline2d_interpol(log10(hc_params%X_E*hc_params%n_gas), &
                            log10(hc_Tgas), spl_NII, extrapolate=.false.))
  else
    val = calc_line_cooling_rate(molecule_NII, hc_params%X_NII)
  end if
end function cooling_NII


function cooling_FeII() result(val)
  double precision val
  double precision, parameter :: min_X_FeII_for_cooling = 1D-15
  if ((hc_params%X_FeII .le. min_X_FeII_for_cooling) .or. &
      (hc_params%X_E .le. 0D0) .or. &
      (hc_params%n_gas .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    val = 0D0
    return
  end if
  if (heating_cooling_config%IonCoolingWithLut) then
    val = hc_params%X_FeII * hc_params%n_gas * &
          exp(phy_ln10 * &  ! = 10**()
          spline2d_interpol(log10(hc_params%X_E*hc_params%n_gas), &
                            log10(hc_Tgas), spl_FeII, extrapolate=.false.))
  else
    val = calc_line_cooling_rate(molecule_FeII, hc_params%X_FeII)
  end if
end function cooling_FeII


function cooling_SiII() result(val)
  double precision val
  double precision, parameter :: min_X_SiII_for_cooling = 1D-15
  if ((hc_params%X_SiII .le. min_X_SiII_for_cooling) .or. &
      (hc_params%X_E .le. 0D0) .or. &
      (hc_params%n_gas .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    val = 0D0
    return
  end if
  if (heating_cooling_config%IonCoolingWithLut) then
    val = hc_params%X_SiII * hc_params%n_gas * &
          exp(phy_ln10 * &  ! = 10**()
          spline2d_interpol(log10(hc_params%X_E*hc_params%n_gas), &
                            log10(hc_Tgas), spl_SiII, extrapolate=.false.))
  else
    val = calc_line_cooling_rate(molecule_SiII, hc_params%X_SiII)
  end if
end function cooling_SiII


function calc_line_cooling_rate(mol, X) result(val)
  double precision val
  type(type_molecule_energy_set), target, intent(inout) :: mol
  double precision, intent(in) :: X
  integer i
  double precision, dimension(:), allocatable :: RSAV
  integer, dimension(:), allocatable :: ISAV
  !
  if ((hc_Tgas .le. 0D0) .or. isnan(hc_Tgas) .or. (mol%n_level .le. 0)) then
    val = 0D0
    return
  end if
  !
  mol_sta_sol => mol
  mol_sta_sol%density_mol = hc_params%n_gas * X
  !
  call heating_cooling_prepare_molecule
  !
  if (heating_cooling_config%solve_method .eq. 1) then
    ! Save the internal common block used by both the DLSODES and DLSODE solver.
    !     Otherwise the content of the common block will be modified by the DLSODE
    !     solver inside statistic_equil_solve, and the DLSODES in the outside
    !     will not work properly.
    !
    allocate(RSAV(256), ISAV(128))
    CALL DSRCMS(RSAV,ISAV, 1)  ! Save
    call statistic_equil_solve
    CALL DSRCMS(RSAV,ISAV, 2)  ! Restore
    deallocate(RSAV, ISAV)
    !
  else if (heating_cooling_config%solve_method .eq. 2) then
    !
    call statistic_equil_solve_Newton
    !
  !else ! LTE
  end if
  !write(*, '(A, A)') mol_sta_sol%name_molecule, ' is solved.'
  !do i=1, mol_sta_sol%n_level
  !  write(*, '(2X, I4, ES12.3)') i, mol_sta_sol%f_occupation(i)
  !end do
  call calc_cooling_rate
  !
  val = mol_sta_sol%cooling_rate_total
  !
  !write(*, '(A, ES12.3)') 'Cooling rate = ', val
  !write(*, '(A, ES12.3)') 'X = ', X
  !write(*,*)
  nullify(mol_sta_sol)
end function calc_line_cooling_rate


function cooling_OI_analytical()
  ! Rollig 2006, equation A.5 and A.6
  ! Z: abundance of O in 3D-4
  ! Need to use the real abundance of OI.
  ! Hollenbach 1989
  ! O 6300: Tielens 2005, equation 2.69
  double precision cooling_OI_analytical
  double precision cooling_OI_63, cooling_OI_146, cooling_OI_6300
  double precision beta_63, beta_146, tau_63, tau_146
  double precision, parameter :: beta_63_N0 = 4.9D20, beta_146_N0 = 3.7D20
  double precision, parameter :: nu_OI_6300 = 4.7D14, Aul_OI_6300 = 6.5D-3 + 2.1D-3
  double precision t1, t2, t3, n_cr_E, n_cr_HI
  !
  if (hc_Tgas .le. 0D0) then
    cooling_OI_analytical = 0D0
    return
  end if
  associate( &
        n_gas => hc_params%n_gas, &
        Tgas  => hc_Tgas, &
        Z     => hc_params%X_OI / 3.2D-4, &
        N     => min(hc_params%Ncol_toISM, &
                     hc_params%Ncol_toStar, &
                     hc_params%n_gas * &
                     hc_params%coherent_length))
    tau_63 = N * Z / beta_63_N0
    tau_146 = N * Z / beta_146_N0
    beta_63 = tau2beta(tau_63)
    beta_146 = tau2beta(tau_146)
    !
    t1 = log(Tgas)
    t2 = exp(0.45D0 * t1)
    t3 = exp(0.66D0 * t1)
    associate( &
      tmp1 => n_gas + beta_63  * 1.66D-5 / (1.35D-11 * t2), &
      tmp2 => n_gas + beta_146 * 8.46D-5 / (4.37D-12 * t3), &
      tmp3 => exp(98D0/Tgas), &
      tmp4 => exp(228D0/Tgas))
      associate(tmp5 => n_gas*n_gas + &
                        tmp3 * tmp1 * (3D0*n_gas + tmp4*5D0*tmp2))
        cooling_OI_63  = 3.15D-14 * 8.46D-5 * beta_63  * Z * &
          3.2D-4 * n_gas * tmp3 * 3D0 * n_gas * tmp1 / tmp5
        cooling_OI_146 = 1.35D-14 * 1.66D-5 * beta_146 * Z * &
          3.2D-4 * n_gas * n_gas * n_gas / tmp5
        !cooling_OI_6300 = 1.8D-24 * hc_params%X_OI &
        !  * n_gas * n_gas * exp(-22800D0/Tgas)
        !
        n_cr_E = 1.3D6 * (Tgas/1D4)**(-0.58)
        n_cr_HI = 6.6D9
        cooling_OI_6300 = &
            phy_hPlanck_SI * nu_OI_6300 * Aul_OI_6300 * hc_params%X_OI * &
            (hc_params%X_E/n_cr_E + hc_params%X_HI/n_cr_HI) * n_gas**2
        !
        cooling_OI_analytical = cooling_OI_63 + cooling_OI_146 + cooling_OI_6300
      end associate
    end associate
  end associate
end function cooling_OI_analytical


function cooling_CII_analytical()
  ! Rollig 2006, equation A.2
  ! Z: Abudance of C in 1.4D-4
  ! Need to use the real abundance of CII
  ! Note that the value used in Rollig 2006 is 1.4D-4.
  ! I guess this is because they normalize relative to n(H2).
  ! Hollenbach 1989
  double precision cooling_CII_analytical
  double precision beta, tau
  double precision, parameter :: beta_N0 = 6.5D20
  if (hc_Tgas .le. 0D0) then
    cooling_CII_analytical = 0D0
    return
  end if
  associate( &
        Tgas    => hc_Tgas, &
        n_gas   => hc_params%n_gas, &
        Z       => hc_params%X_CII/1.4D-4, &
        N       => min(hc_params%Ncol_toISM, &
                       hc_params%Ncol_toStar, &
                       hc_params%n_gas * &
                       hc_params%coherent_length))
    tau =  N * Z / beta_N0
    beta = tau2beta(tau)
    associate( &
        cooling_CII_158 => 4.04D-24 * n_gas * Z * beta / &
          (1D0 + 0.5D0 * exp(92D0/Tgas) * (1D0 + 2600D0*beta/n_gas)))
      cooling_CII_analytical = cooling_CII_158
    end associate
  end associate
end function cooling_CII_analytical


function cooling_Neufeld_H2O_rot()
  use load_Neufeld_cooling_H2O
  double precision cooling_Neufeld_H2O_rot
  double precision t1
  !
  if ((hc_params%X_H2O .le. 0D0) .or. &
      (hc_params%X_H2 .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    cooling_Neufeld_H2O_rot = 0D0
    return
  end if
  a_Neufeld_cooling_H2O_params%T = hc_Tgas
  associate( &
        L0    => a_Neufeld_cooling_H2O_params%L0, &
        L_LTE => a_Neufeld_cooling_H2O_params%L_LTE, &
        n_12  => a_Neufeld_cooling_H2O_params%n_12, &
        alpha => a_Neufeld_cooling_H2O_params%alpha, &
        log10N=> a_Neufeld_cooling_H2O_params%log10N, &
        G     => hc_params%Neufeld_G, &
        n_M   => hc_params%n_gas * hc_params%X_H2O, &
        n_H2  => hc_params%n_gas * hc_params%X_H2, &
        dv_dz => hc_params%Neufeld_dv_dz)
    log10N = log10(min(G * n_M / (dv_dz + very_small_num), &
      n_M*hc_params%Ncol_toISM/hc_params%n_gas / &
      (9D0*hc_params%velo_width_turb*1D-5)))
    L0    = get_L0() + very_small_num
    L_LTE = get_L_LTE() + very_small_num
    n_12  = get_n_12() + very_small_num
    alpha = get_alpha()
    !
    t1 = exp(alpha * log(n_H2/n_12))
    !
    cooling_Neufeld_H2O_rot = n_H2 * n_M / &
      (1D0/L0 + n_H2/L_LTE + &
       1D0/L0 * t1 * (1D0 - n_12*L0/L_LTE))
  end associate
end function cooling_Neufeld_H2O_rot


function cooling_Neufeld_H2O_vib()
  use load_Neufeld_cooling_H2O
  double precision cooling_Neufeld_H2O_vib
  !
  if ((hc_params%X_H2O .le. 0D0) .or. &
      (hc_params%X_H2 .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    cooling_Neufeld_H2O_vib = 0D0
    return
  end if
  a_Neufeld_cooling_H2O_params%T = hc_Tgas
  associate( &
    L0        => a_Neufeld_cooling_H2O_params%L0_vib, &
    L_LTE     => a_Neufeld_cooling_H2O_params%L_LTE_vib, &
    log10N    => a_Neufeld_cooling_H2O_params%log10N, &
    G         => hc_params%Neufeld_G, &
    n_M       => hc_params%n_gas * hc_params%X_H2O, &
    n_H2      => hc_params%n_gas * hc_params%X_H2, &
    dv_dz     => hc_params%Neufeld_dv_dz)
    !
    log10N = log10(min(G * n_M / (dv_dz + very_small_num), &
      n_M*hc_params%Ncol_toISM/hc_params%n_gas / &
      (9D0*hc_params%velo_width_turb*1D-5)))
    L0    = get_L0_vib() + very_small_num
    L_LTE = get_L_LTE_vib() + very_small_num
    !
    cooling_Neufeld_H2O_vib = n_H2 * n_M / (1D0/L0 + n_H2/L_LTE)
  end associate
end function cooling_Neufeld_H2O_vib


function cooling_Neufeld_CO_rot()
  use load_Neufeld_cooling_CO
  double precision cooling_Neufeld_CO_rot
  !
  if ((hc_params%X_CO .le. 0D0) .or. &
      (hc_params%X_H2 .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    cooling_Neufeld_CO_rot = 0D0
    return
  end if
  a_Neufeld_cooling_CO_params%T = hc_Tgas
  associate( &
    ! L is the cooling coefficient.
    ! L     => a_Neufeld_cooling_CO_params%L, &
    ! The following four are to be interpolated or extrapolated from the Neufeld tables.
    L0    => a_Neufeld_cooling_CO_params%L0, &
    L_LTE => a_Neufeld_cooling_CO_params%L_LTE, &
    n_12  => a_Neufeld_cooling_CO_params%n_12, &
    alpha => a_Neufeld_cooling_CO_params%alpha, &
    log10N=> a_Neufeld_cooling_CO_params%log10N, &
    G     => hc_params%Neufeld_G, &
    n_M   => hc_params%n_gas * hc_params%X_CO, &
    n_H2  => hc_params%n_gas * hc_params%X_H2, &
    dv_dz => hc_params%Neufeld_dv_dz)
    !
    log10N = log10(min(G * n_M / (dv_dz + very_small_num), &
      n_M*hc_params%Ncol_toISM/hc_params%n_gas / &
      (9D0*hc_params%velo_width_turb*1D-5)))
    L0    = get_L0() + very_small_num
    L_LTE = get_L_LTE() + very_small_num
    n_12  = get_n_12() + very_small_num
    alpha = get_alpha()
    !
    cooling_Neufeld_CO_rot = n_H2 * n_M / &
      (1D0/L0 + n_H2/L_LTE + &
       1D0/L0 * (n_H2/n_12)**alpha * (1D0 - n_12*L0/L_LTE))
  end associate
end function cooling_Neufeld_CO_rot


function cooling_Neufeld_CO_vib()
  use load_Neufeld_cooling_CO
  double precision cooling_Neufeld_CO_vib
  !
  if ((hc_params%X_CO .le. 0D0) .or. &
      (hc_params%X_H2 .le. 0D0) .or. &
      (hc_Tgas .le. 0D0)) then
    cooling_Neufeld_CO_vib = 0D0
    return
  end if
  a_Neufeld_cooling_CO_params%T = hc_Tgas
  associate( &
    L0     => a_Neufeld_cooling_CO_params%L0_vib, &
    L_LTE  => a_Neufeld_cooling_CO_params%L_LTE_vib, &
    log10N => a_Neufeld_cooling_CO_params%log10N, &
    G      => hc_params%Neufeld_G, &
    n_M    => hc_params%n_gas * hc_params%X_CO, &
    n_H2   => hc_params%n_gas * hc_params%X_H2, &
    dv_dz  => hc_params%Neufeld_dv_dz)
    !
    log10N = log10(min(G * n_M / (dv_dz + very_small_num), &
      n_M*hc_params%Ncol_toISM/hc_params%n_gas / &
      (9D0*hc_params%velo_width_turb*1D-5)))
    L0    = get_L0_vib() + very_small_num
    L_LTE = get_L_LTE_vib() + very_small_num
    !
    cooling_Neufeld_CO_vib = n_H2 * n_M / (1D0/L0 + n_H2/L_LTE)
  end associate
end function cooling_Neufeld_CO_vib


function cooling_OH_rot()
  double precision cooling_OH_rot
  !
  ! Hollenbach 1989, equation 6.21
  !
  ! Gorti 2004, appendix D
  double precision, parameter :: A0 = 7.6D-4 ! s-1
  double precision, parameter :: E0 = 5.4D0  ! K
  double precision, parameter :: sig = 8D-16 ! cm-2
  double precision, parameter :: eta = 10D0
  double precision N, N_tau, tau, ctau, v_T, tmp, tmp1, ym, L
  !
  if ((hc_params%X_OH .le. 0D0) .or. &
      (hc_params%X_H2 .lt. 0D0) .or. &
      (hc_params%X_H2 .ge. 1D0) .or. &
      (hc_Tgas .le. 0D0)) then
    cooling_OH_rot = 0D0
    return
  end if
  !
  N = hc_params%X_OH * hc_params%n_gas * hc_params%coherent_length
  N_tau = 1.18D7 * hc_params%velo_width_turb*1D-5 * E0**3 / A0
  tau = 4D0 * N / N_tau / (eta * hc_Tgas/E0)
  ctau = tau * sqrt(phy_2Pi * log(2.13D0 + (tau/exp(1D0))**2))
  !
  v_T = sqrt((8D0/phy_Pi*phy_kBoltzmann_CGS/phy_mProton_CGS) * hc_Tgas)
  tmp = 4D0 * (hc_Tgas / E0) * A0 / (hc_params%n_gas * (1D0 - hc_params%X_H2) * sig * v_T)
  ym = log(1D0 + ctau / (1D0 + 10D0 * tmp))
  tmp1 = (2D0 + ym + 0.6D0 * ym**2) / (1D0 + ctau + tmp + 1.5D0 * sqrt(tmp))
  L = 2D0 * phy_kBoltzmann_CGS*hc_Tgas**2 * A0 / E0 * tmp1
  cooling_OH_rot = L * hc_params%n_gas * hc_params%X_OH
end function cooling_OH_rot


function heating_minus_cooling()
  double precision heating_minus_cooling
  !
  associate(r => heating_cooling_rates)
    r%heating_photoelectric_small_grain_rate = heating_photoelectric_small_grain()
    r%heating_formation_H2_rate              = heating_formation_H2()
    r%heating_cosmic_ray_rate                = heating_cosmic_ray()
    r%heating_vibrational_H2_rate            = heating_vibrational_H2()
    r%heating_ionization_CI_rate             = heating_ionization_CI()
    r%heating_photodissociation_H2_rate      = heating_photodissociation_H2()
    r%heating_photodissociation_H2O_rate     = heating_photodissociation_H2O()
    r%heating_photodissociation_OH_rate      = heating_photodissociation_OH()
    r%heating_Xray_Bethell_rate              = heating_Xray_Bethell()
    r%heating_viscosity_rate                 = heating_viscosity()
    r%heating_chem                           = heating_chemical()
    r%cooling_photoelectric_small_grain_rate = cooling_photoelectric_small_grain()
    r%cooling_vibrational_H2_rate            = cooling_vibrational_H2()
    r%cooling_gas_grain_collision_rate       = cooling_gas_grain_collision()
    r%cooling_OI_rate                        = cooling_OI()
    r%cooling_CII_rate                       = cooling_CII()
    r%cooling_NII_rate                       = cooling_NII()
    r%cooling_SiII_rate                      = cooling_SiII()
    r%cooling_FeII_rate                      = cooling_FeII()
    r%cooling_OH_rot_rate                    = cooling_OH_rot()
    r%cooling_Neufeld_H2O_rate_rot           = cooling_Neufeld_H2O_rot()
    r%cooling_Neufeld_H2O_rate_vib           = cooling_Neufeld_H2O_vib()
    r%cooling_Neufeld_CO_rate_rot            = cooling_Neufeld_CO_rot()
    r%cooling_Neufeld_CO_rate_vib            = cooling_Neufeld_CO_vib()
    r%cooling_Neufeld_H2_rot_rate            = cooling_Neufeld_H2_rot()
    r%cooling_LymanAlpha_rate                = cooling_LymanAlpha()
    r%cooling_free_bound_rate                = cooling_free_bound()
    r%cooling_free_free_rate                 = cooling_free_free()
    !
    heating_minus_cooling = &
        r%heating_photoelectric_small_grain_rate &  ! 1
      + r%heating_formation_H2_rate &               ! 2
      + r%heating_cosmic_ray_rate &                 ! 3
      + r%heating_vibrational_H2_rate &             ! 4
      + r%heating_ionization_CI_rate &              ! 5
      + r%heating_photodissociation_H2_rate &       ! 6
      + r%heating_photodissociation_H2O_rate &      ! 7
      + r%heating_photodissociation_OH_rate &       ! 8
      + r%heating_Xray_Bethell_rate         &       ! 9
      + r%heating_viscosity_rate            &       ! 10
      + r%heating_chem                      &       ! 11
      !
      - r%cooling_photoelectric_small_grain_rate &  ! 1
      - r%cooling_vibrational_H2_rate &             ! 2
      - r%cooling_gas_grain_collision_rate &        ! 3
      - r%cooling_OI_rate &                         ! 4
      - r%cooling_CII_rate &                        ! 5
      - r%cooling_Neufeld_H2O_rate_rot &            ! 6
      - r%cooling_Neufeld_H2O_rate_vib &            ! 7
      - r%cooling_Neufeld_CO_rate_rot &             ! 8
      - r%cooling_Neufeld_CO_rate_vib &             ! 9
      - r%cooling_Neufeld_H2_rot_rate &             ! 10
      - r%cooling_LymanAlpha_rate &                 ! 11
      - r%cooling_free_bound_rate &                 ! 12
      - r%cooling_free_free_rate  &                 ! 13
      - r%cooling_NII_rate        &                 ! 14
      - r%cooling_SiII_rate       &                 ! 15
      - r%cooling_FeII_rate       &                 ! 16
      - r%cooling_OH_rot_rate                       ! 17
    r%hc_net_rate = heating_minus_cooling
  end associate
end function heating_minus_cooling



function solve_bisect_T(T0, n_iter, converged)
  implicit none
  double precision solve_bisect_T
  double precision, intent(in) :: T0
  integer, intent(out) :: n_iter
  logical, intent(out) :: converged
  double precision x1, x2, f1, f2, dx, xmid, fmid
  integer i, j
  integer nmax_expand
  integer nmax_shrink
  double precision expand_factor
  double precision expand_factor_tmp, x_tmp
  logical found_diff_sign
  double precision rtol
  double precision atol
  !
  nmax_expand = 2048
  nmax_shrink = 1024
  expand_factor = 0.5D0
  rtol = 1D-5
  atol = 1D-1
  !
  ! Local copy of the gas and dust temperature
  hc_Tgas = hc_params%Tgas
  hc_Tdust = hc_params%Tdust
  !
  x1 = T0 / 1.1D0
  x2 = T0 * 1.1D0
  hc_Tgas = x1
  f1 = heating_minus_cooling()
  hc_Tgas = x2
  f2 = heating_minus_cooling()
  found_diff_sign = .false.
  converged = .false.
  solve_bisect_T = -1D0
  do i=1, nmax_expand
    if (f1*f2 .le. 0D0) then
      found_diff_sign = .true.
      exit
    end if
    if (abs(f1) .lt. abs(f2)) then
      x1 = max(1D0, x1 + expand_factor * (x1 - x2))
      hc_Tgas = x1
      f1 = heating_minus_cooling()
    else if (abs(f1) .ge. abs(f2)) then
      x2 = max(1D0, x2 + expand_factor * (x2 - x1))
      hc_Tgas = x2
      f2 = heating_minus_cooling()
    else
      if (isnan(f1)) then
        x1 = max(1D0, (x1 + expand_factor*x2)/(1D0+expand_factor))
        expand_factor_tmp = expand_factor / 2D0
        do j=1, 30
          x_tmp = max(1D0, x1 + expand_factor_tmp * (x1 - x2))
          hc_Tgas = x_tmp
          f1 = heating_minus_cooling()
          if (.not. isnan(f1)) then
            x1 = x_tmp
            exit
          else
            x_tmp = x1
            expand_factor_tmp  = expand_factor_tmp / 2D0
          end if
        end do
      else
        x2 = max(1D0, (x2 + expand_factor*x1)/(1D0+expand_factor))
        expand_factor_tmp = expand_factor / 2D0
        do j=1, 30
          x_tmp = max(1D0, x2 + expand_factor_tmp * (x2 - x1))
          hc_Tgas = x_tmp
          f2 = heating_minus_cooling()
          if (.not. isnan(f2)) then
            x2 = x_tmp
            exit
          else
            x_tmp = x2
            expand_factor_tmp  = expand_factor_tmp / 2D0
          end if
        end do
      end if
    end if
  end do
  if (.not. found_diff_sign) then
    write(*,*) 'Fail to find an appropriate bracketing interval.'
    converged = .false.
    n_iter = 0
    return
  else
    if (f1 .eq. 0D0) then
      solve_bisect_T = x1
      converged = .true.
      n_iter = 0
      return
    else if (f2 .eq. 0D0) then
      solve_bisect_T = x2
      converged = .true.
      n_iter = 0
      return
    else
      do i=1, nmax_shrink
        dx = x2 - x1
        xmid = 0.5D0 * (x1+x2)
        if (abs(dx) .lt. (rtol*xmid + atol)) then
          solve_bisect_T = xmid
          converged = .true.
          n_iter = i
          return
        end if
        hc_Tgas = xmid
        fmid = heating_minus_cooling()
        if (isnan(fmid)) then
          xmid = x2 - dx * 0.1D0
          hc_Tgas = xmid
          fmid = heating_minus_cooling()
        end if
        if (fmid*f1 .lt. 0D0) then
          x2 = xmid
          f2 = fmid
        else if (fmid .eq. 0D0) then
          solve_bisect_T = xmid
          converged = .true.
          n_iter = i
          return
        else
          x1 = xmid
          f1 = fmid
        end if
      end do
    end if
  end if
end function solve_bisect_T


subroutine print_out_h_c_rates(Tmin, Tmax, dT0, dT_ratio)
  double precision, intent(in) :: Tmin, Tmax, dT0, dT_ratio
  double precision dT, hc
  dT = dT0
  hc_Tgas = Tmin
  do
    hc = heating_minus_cooling()
    write(*,'(25ES12.4)') hc_Tgas, hc, &
        heating_cooling_rates%heating_photoelectric_small_grain_rate, &
        heating_cooling_rates%heating_formation_H2_rate, &
        heating_cooling_rates%heating_cosmic_ray_rate, &
        heating_cooling_rates%heating_vibrational_H2_rate, &
        heating_cooling_rates%heating_ionization_CI_rate, &
        heating_cooling_rates%heating_photodissociation_H2_rate, &
        heating_cooling_rates%heating_photodissociation_H2O_rate, &
        heating_cooling_rates%heating_photodissociation_OH_rate, &
        heating_cooling_rates%heating_Xray_Bethell_rate, &
        heating_cooling_rates%heating_viscosity_rate, &
        heating_cooling_rates%cooling_photoelectric_small_grain_rate, &
        heating_cooling_rates%cooling_vibrational_H2_rate, &
        heating_cooling_rates%cooling_gas_grain_collision_rate, &
        heating_cooling_rates%cooling_OI_rate, &
        heating_cooling_rates%cooling_CII_rate, &
        heating_cooling_rates%cooling_Neufeld_H2O_rate_rot, &
        heating_cooling_rates%cooling_Neufeld_H2O_rate_vib, &
        heating_cooling_rates%cooling_Neufeld_CO_rate_rot, &
        heating_cooling_rates%cooling_Neufeld_CO_rate_vib, &
        heating_cooling_rates%cooling_Neufeld_H2_rot_rate, &
        heating_cooling_rates%cooling_LymanAlpha_rate, &
        heating_cooling_rates%cooling_free_bound_rate, &
        heating_cooling_rates%cooling_free_free_rate
    if (hc_Tgas .gt. Tmax) then
      exit
    end if
    hc_Tgas = hc_Tgas + dT
    dT = dT * dT_ratio
  end do
end subroutine print_out_h_c_rates


function max_heating_rate()
  double precision max_heating_rate
  associate(r => heating_cooling_rates)
    max_heating_rate = max( &
      r%heating_photoelectric_small_grain_rate , &
      r%heating_formation_H2_rate              , &
      r%heating_cosmic_ray_rate                , &
      r%heating_vibrational_H2_rate            , &
      r%heating_ionization_CI_rate             , &
      r%heating_photodissociation_H2_rate      , &
      r%heating_photodissociation_H2O_rate     , &
      r%heating_photodissociation_OH_rate      , &
      r%heating_Xray_Bethell_rate              , &
      r%heating_viscosity_rate                 , &
      r%heating_chem                           , &
      -r%cooling_photoelectric_small_grain_rate)
  end associate
end function max_heating_rate


function max_cooling_rate()
  double precision max_cooling_rate
  associate(r => heating_cooling_rates)
    max_cooling_rate = max( &
      r%cooling_photoelectric_small_grain_rate , &
      r%cooling_vibrational_H2_rate            , &
      r%cooling_gas_grain_collision_rate       , &
      r%cooling_OI_rate                        , &
      r%cooling_CII_rate                       , &
      r%cooling_Neufeld_H2O_rate_rot           , &
      r%cooling_Neufeld_H2O_rate_vib           , &
      r%cooling_Neufeld_CO_rate_rot            , &
      r%cooling_Neufeld_CO_rate_vib            , &
      r%cooling_Neufeld_H2_rot_rate            , &
      r%cooling_LymanAlpha_rate                , &
      r%cooling_free_bound_rate                , &
      r%cooling_free_free_rate                 , &
      r%cooling_NII_rate                       , &
      r%cooling_SiII_rate                      , &
      r%cooling_FeII_rate                      , &
      r%cooling_OH_rot_rate                    , &
      -r%heating_chem                          )
  end associate
end function max_cooling_rate


function heating_cooling_is_very_slow() result(is_slow)
  logical is_slow
  double precision, parameter :: thresh_hc_slow = 1D-2
  is_slow = &
    heating_cooling_rates%hc_net_rate .lt. &
    thresh_hc_slow * max(max_heating_rate(), max_cooling_rate())
end function heating_cooling_is_very_slow



subroutine disp_h_c_rates
  associate(r => heating_cooling_rates)
  write(*,*) '1 ', r%heating_photoelectric_small_grain_rate    ! 1
  write(*,*) '2 ', r%heating_formation_H2_rate                 ! 2
  write(*,*) '3 ', r%heating_cosmic_ray_rate                   ! 3
  write(*,*) '4 ', r%heating_vibrational_H2_rate               ! 4
  write(*,*) '5 ', r%heating_ionization_CI_rate                ! 5
  write(*,*) '6 ', r%heating_photodissociation_H2_rate         ! 6
  write(*,*) '7 ', r%heating_photodissociation_H2O_rate        ! 7
  write(*,*) '8 ', r%heating_photodissociation_OH_rate         ! 8
  write(*,*) '9 ', r%heating_Xray_Bethell_rate                 ! 9
  write(*,*) '10', r%heating_viscosity_rate                    ! 10
  write(*,*) '11', r%heating_chem                              ! 11
  write(*,*) '1 ', r%cooling_photoelectric_small_grain_rate    ! 1
  write(*,*) '2 ', r%cooling_vibrational_H2_rate               ! 2
  write(*,*) '3 ', r%cooling_gas_grain_collision_rate          ! 3
  write(*,*) '4 ', r%cooling_OI_rate                           ! 4
  write(*,*) '5 ', r%cooling_CII_rate                          ! 5
  write(*,*) '6 ', r%cooling_Neufeld_H2O_rate_rot              ! 6
  write(*,*) '7 ', r%cooling_Neufeld_H2O_rate_vib              ! 7
  write(*,*) '8 ', r%cooling_Neufeld_CO_rate_rot               ! 8
  write(*,*) '9 ', r%cooling_Neufeld_CO_rate_vib               ! 9
  write(*,*) '10', r%cooling_Neufeld_H2_rot_rate               ! 10
  write(*,*) '11', r%cooling_LymanAlpha_rate                   ! 11
  write(*,*) '12', r%cooling_free_bound_rate                   ! 12
  write(*,*) '13', r%cooling_free_free_rate                    ! 13
  write(*,*)
  end associate
end subroutine disp_h_c_rates


end module heating_cooling
