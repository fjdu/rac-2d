module data_struct


implicit none


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
  integer itype
  integer n
  double precision, dimension(:), allocatable :: lam
  double precision, dimension(:), allocatable :: ab, sc, tot, g
end type type_optical_property



type :: type_global_material_collection
  integer ntype
  type(type_optical_property), dimension(:), allocatable :: list
end type type_global_material_collection



type :: type_local_encounter_collection
  integer ntype, nlam
  double precision en_gain, kph
  integer ph_count
  double precision, dimension(:), allocatable :: X
  double precision, dimension(:,:), allocatable :: acc
  double precision, dimension(:), allocatable :: phweight
  double precision, dimension(:), allocatable :: summed
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
  double precision lumi, mass, radius, T
  double precision, dimension(:), allocatable :: lam, vals
end type type_stellar_spectrum



type :: type_montecarlo_config
  integer nph
  double precision eph
  integer icount, nmax_cross, nmax_encounter
  integer fU
  character(len=128) fname_photons, fname_dust, fname_star, mc_dir_in, mc_dir_out
  double precision minw, maxw, min_ang, max_ang
  logical use_blackbody_star
  double precision star_mass, star_radius, star_temperature
end type type_montecarlo_config



type :: type_ray
  double precision x, y, z, vx, vy, vz
end type type_ray


type :: type_photon_packet
  type(type_ray) :: ray
  double precision lam, en
  integer iKap, iSpec
  integer e_count
end type type_photon_packet



type :: type_cell_rz_phy_basic
  double precision rmin, rmax, rcen, dr, zmin, zmax, zcen, dz, daz, volume
  double precision :: &
    Tgas, &
    Tdust, &
    Tdust1, &
    n_gas, &
    n_dust, &
    mdust, &
    mdust_cell, &
    !
    UV_G0_factor, &
    UV_G0_factor_background, &
    LymanAlpha_G0_factor, &
    LymanAlpha_number_flux_0, &
    LymanAlpha_energy_flux_0, &
    Xray_flux_0, &
    Av, &
    Ncol, &
    dNcol, &
    omega_albedo, &
    zeta_cosmicray_H2, &
    !
    stickCoeffH, &
    R_H2_form_rate_coeff, &
    R_H2_form_rate, &
    !
    f_selfshielding_H2, &
    f_selfshielding_CO, &
    f_selfshielding_H2O, &
    f_selfshielding_OH, &
    !
    GrainMaterialDensity_CGS, &
    GrainRadius_CGS, &
    aGrainMin_CGS, &
    aGrainMax_CGS, &
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
    !
    alpha_viscosity, &
    !
    t_final
end type type_cell_rz_phy_basic


!type :: type_heating_cooling_rate_one
!  character(len=12) h_c_name
!  double precision val
!end type type_heating_cooling_rate_one
!
!
!
!type :: type_heating_cooling_rates_list
!  integer nitem
!  type(type_heating_cooling_rate_one) :: &
!    heating_photoelectric_small_grain_rate, &
!    heating_formation_H2_rate, &
!    heating_cosmic_ray_rate, &
!    heating_vibrational_H2_rate, &
!    heating_ionization_CI_rate, &
!    heating_photodissociation_H2_rate, &
!    heating_photodissociation_H2O_rate, &
!    heating_photodissociation_OH_rate, &
!    heating_Xray_Bethell_rate, &
!    heating_viscosity_rate, &
!    cooling_photoelectric_small_grain_rate, &
!    cooling_vibrational_H2_rate, &
!    cooling_gas_grain_collision_rate, &
!    cooling_OI_rate, &
!    cooling_CII_rate, &
!    cooling_Neufeld_H2O_rate_rot, &
!    cooling_Neufeld_H2O_rate_vib, &
!    cooling_Neufeld_CO_rate_rot, &
!    cooling_Neufeld_CO_rate_vib, &
!    cooling_Neufeld_H2_rot_rate, &
!    cooling_LymanAlpha_rate, &
!    cooling_free_bound_rate, &
!    cooling_free_free_rate
!end type type_heating_cooling_rates_list


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


type, private :: type_neighbor
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
  integer :: iIter = 0
  integer :: quality = 0
  type(type_local_encounter_collection) :: optical
end type type_cell


end module data_struct

