module data_struct


implicit none

integer, parameter :: LongInt = 8

integer, parameter :: MaxNumOfDustComponents = 4

type :: type_point
  double precision x, y, z
end type type_point


type :: type_direction_cartesian
  double precision u, v, w
end type type_direction_cartesian


type :: type_direction_sph
  double precision theta, phi
end type type_direction_sph


type :: type_sphere_coor_quat
  double precision costheta, sintheta, cosphi, sinphi
end type type_sphere_coor_quat


type :: type_spectrum_generic
  integer n
  double precision, dimension(:,:), allocatable :: intervals
  double precision, dimension(:), allocatable :: vals
end type type_spectrum_generic



type :: type_distribution_table
  integer n
  double precision, dimension(:), allocatable :: pvals
end type type_distribution_table


type :: type_optical_property
  integer n
  double precision, dimension(:), allocatable :: lam
  double precision, dimension(:), allocatable :: ab, sc, g
end type type_optical_property



type :: type_global_material_collection
  integer ntype
  type(type_optical_property), dimension(:), allocatable :: list
end type type_global_material_collection



type :: type_local_encounter_collection
  integer ntype, nlam
  integer cr_count
  double precision, dimension(:), allocatable :: X
  double precision, dimension(:), allocatable :: summed_ab, summed_sc
  double precision, dimension(:), allocatable :: summed
  double precision, dimension(:,:), allocatable :: acc
  double precision, dimension(:), allocatable :: flux
  type(type_direction_cartesian), dimension(:), allocatable :: dir_wei
end type type_local_encounter_collection


type :: type_LUT_Tdust
  ! A look-up table
  ! Given a val, find out the corresponding T
  integer n, m
  double precision, dimension(:), allocatable :: Tds, vals
  double precision, dimension(:,:), allocatable :: table
end type type_LUT_Tdust



type :: type_stellar_spectrum
  integer n
  double precision lumi, mass, radius, T, lumi_UV, lumi0, lumi_UV0
  double precision, dimension(:), allocatable :: lam, vals, vals0
end type type_stellar_spectrum



type :: type_montecarlo_config
  double precision eph
  integer(kind=LongInt) nph, icount, nmax_cross, nmax_encounter
  integer fU
  character(len=128) fname_photons, fname_dust, fname_water, fname_star, mc_dir_in, mc_dir_out
  double precision minw, maxw, min_ang, max_ang
  logical use_blackbody_star, savephoton
end type type_montecarlo_config



type :: type_ray
  double precision x, y, z, vx, vy, vz
end type type_ray


type :: type_photon_packet
  type(type_ray) :: ray
  double precision lam, en
  double precision f, Inu
  integer iKap, iSpec, iTran
  integer e_count
end type type_photon_packet


type :: type_mole_f_occ
  integer nlevels
  double precision, dimension(:), allocatable :: vals
end type type_mole_f_occ


type :: type_cell_rz_phy_basic
  double precision rmin, rmax, rcen, dr, zmin, zmax, zcen, dz, daz
  double precision volume, surf_area, area_T, area_B, area_I, area_O
  integer ndustcompo
  double precision :: &
    Tgas, &
    Tdust, &
    !
    n_gas, &
    !n_dust, &
    !
    mgas_cell, &
    !
    Tdusts(MaxNumOfDustComponents), &
    rho_dusts(MaxNumOfDustComponents), &
    n_dusts(MaxNumOfDustComponents), &
    mp_dusts(MaxNumOfDustComponents), &
    mdusts_cell(MaxNumOfDustComponents), &
    !
    abso_wei(MaxNumOfDustComponents), &
    !
    en_gains(MaxNumOfDustComponents), &
    en_gains_abso(MaxNumOfDustComponents), &
    en_prevs(MaxNumOfDustComponents), &
    kphs(MaxNumOfDustComponents), &
    !
    en_gain_tot, &
    en_gain_abso_tot, &
    !
    sigdust_ave, &
    ndust_tot, &
    mdust_tot, &
    !
    ! UV_G0_factor, &
    UV_G0_factor_background, &
    ! LymanAlpha_G0_factor, &
    !
    ! LymanAlpha_number_flux_0, &
    ! LymanAlpha_energy_flux_0, &
    Xray_flux_0, &
    Av, &
    Ncol, &
    dNcol, &
    !
    phflux_Lya, &
    !
    flux_UV_star_unatten, &
    !
    G0_UV_toISM, &
    G0_UV_toStar, &
    G0_Lya_atten, &
    !
    Av_toISM, &
    Av_toStar, &
    !
    Ncol_toISM, &
    Ncol_toStar, &
    !
    omega_albedo, &
    zeta_cosmicray_H2, &
    !
    R_H2_form_rate_coeff, &
    R_H2_form_rate, &
    !
    f_selfshielding_H2, &
    f_selfshielding_CO, &
    f_selfshielding_H2O, &
    f_selfshielding_OH, &
    !
    f_selfshielding_toISM_H2, &
    f_selfshielding_toISM_CO, &
    f_selfshielding_toISM_H2O, &
    f_selfshielding_toISM_OH, &
    !
    f_selfshielding_toStar_H2, &
    f_selfshielding_toStar_CO, &
    f_selfshielding_toStar_H2O, &
    f_selfshielding_toStar_OH, &
    !
    SitesPerGrain, &
    GrainMaterialDensity_CGS, &
    GrainRadius_CGS, &
    !
    ratioDust2GasMass, &
    ratioDust2HnucNum, &
    dust_depletion, &
    MeanMolWeight, &
    !
    omega_Kepler, &
    velo_Kepler, &
    velo_gradient, &
    velo_width_turb, &
    coherent_length, &
    sound_speed, &
    !
    alpha_viscosity, &
    !
    t_final
  double precision :: X_H2, X_HI, X_CI, X_Cplus, X_OI, X_CO, &
                      X_H2O, X_OH, X_E, X_Hplus, X_gH
  double precision :: flux_UV, flux_Lya, flux_NIR, flux_MIR, flux_FIR
  double precision :: dir_UV_r, dir_UV_z, dir_Lya_r, dir_Lya_z, &
                      dir_NIR_r, dir_NIR_z, dir_MIR_r, dir_MIR_z, &
                      dir_FIR_r, dir_FIR_z
  double precision :: aniso_UV, aniso_Lya, aniso_NIR, aniso_MIR, aniso_FIR
  double precision :: pressure_thermal, gravity_z, gravity_acc_z
  !
  integer ab_count_dust, ab_count_water
  integer sc_count_dust, sc_count_HI
  !
  double precision ab_en_water
  !
end type type_cell_rz_phy_basic


type :: type_dust_MRN
  double precision :: rmin, rmax, n
  double precision :: rav, r2av, r3av
end type type_dust_MRN


type :: type_Andrews_disk
  logical :: useNumDens = .true.
  double precision :: particlemass = 1.4D0 * 1.67262158D-24
  double precision :: Md=0.01D0 ! Disk mass in Msun
  double precision :: rc=200D0  ! Disk outer boundary
  double precision :: hc=50D0   ! Scale height at outer boundary
  double precision :: gam=1D0   ! Power index for surface density
  double precision :: psi=1D0   ! Power index for scale height
end type type_Andrews_disk


type :: type_heating_cooling_rates_list
  double precision :: &
    hc_net_rate, &
    heating_photoelectric_small_grain_rate, &
    heating_formation_H2_rate, &
    heating_cosmic_ray_rate, &
    heating_vibrational_H2_rate, &
    heating_ionization_CI_rate, &
    heating_photodissociation_H2_rate, &
    heating_photodissociation_H2O_rate, &
    heating_photodissociation_OH_rate, &
    heating_Xray_Bethell_rate, &
    heating_viscosity_rate, &
    heating_chem, &
    cooling_photoelectric_small_grain_rate, &
    cooling_vibrational_H2_rate, &
    cooling_gas_grain_collision_rate, &
    cooling_OI_rate, &
    cooling_CII_rate, &
    cooling_Neufeld_H2O_rate_rot, &
    cooling_Neufeld_H2O_rate_vib, &
    cooling_Neufeld_CO_rate_rot, &
    cooling_Neufeld_CO_rate_vib, &
    cooling_Neufeld_H2_rot_rate, &
    cooling_LymanAlpha_rate, &
    cooling_free_bound_rate, &
    cooling_free_free_rate
end type type_heating_cooling_rates_list


type, private :: type_child_tmp
  type(type_cell), pointer :: p
end type type_child_tmp


type :: type_neighbor
  integer :: n = 0
  integer, dimension(:), allocatable :: idx
  double precision, dimension(:), allocatable :: fra
  double precision :: fra_tot = 0D0
end type type_neighbor

type :: type_cell
  double precision :: xmin=0D0, xmax=0D0, ymin=0D0, ymax=0D0
  double precision, dimension(:), allocatable :: val
  logical :: using = .false., converged = .false.
  integer :: order=0, nChildren=0, nOffspring=0, nleaves=0
  type(type_cell), pointer :: parent => null()
  type(type_child_tmp), pointer, dimension(:) :: children
  type(type_neighbor), pointer :: inner => null(), outer => null(), &
        below => null(), above => null(), around => null()
  type(type_cell_rz_phy_basic), pointer :: par => null()
  type(type_heating_cooling_rates_list), allocatable :: h_c_rates
  double precision, dimension(:), allocatable :: abundances
  double precision, dimension(:), allocatable :: col_den, col_den_acc
  double precision, dimension(:), allocatable :: col_den_toStar, col_den_toISM
  integer :: iIter = 0
  integer :: quality = 0
  type(type_local_encounter_collection) :: optical
  type(type_dust_MRN), allocatable :: mrn
  type(type_mole_f_occ), allocatable :: focc
end type type_cell


end module data_struct

