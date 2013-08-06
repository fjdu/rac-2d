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


type :: type_cell_rz_phy_basic
  double precision rmin, rmax, rcen, dr, zmin, zmax, zcen, dz, daz
  double precision :: &
    Tgas, &
    Tdust, &
    n_gas, &
    n_dust, &
    !
    UV_G0_factor, &
    Av, &
    LymanAlpha_flux_0, &
    Xray_flux_0, &
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
    !
    alpha_viscosity
    !
end type type_cell_rz_phy_basic


type :: type_heating_cooling_rates_list
  double precision :: &
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
  logical :: using = .false.
  integer :: order=0, nChildren=0, nOffspring=0
  type(type_cell), pointer :: parent => null()
  type(type_child_tmp), pointer, dimension(:) :: children
  type(type_neighbor), pointer :: inner, outer, below, above, around
  type(type_cell_rz_phy_basic), pointer :: par => null()
  type(type_heating_cooling_rates_list), allocatable :: h_c_rates
  double precision, dimension(:), allocatable :: abundances
  double precision, dimension(:), allocatable :: col_den, col_den_acc
  integer, allocatable :: iIter
end type type_cell


end module data_struct

