module data_struct

implicit none

integer, parameter :: MaxNumOfGasComponents  = 1
integer, parameter :: MaxNumOfDustComponents = 8

integer, parameter, private :: const_len_energy_level = 12
integer, parameter, public :: const_len_qnum = 15
integer, parameter, private :: const_len_molecule = 12

integer, parameter :: MaxNumOfFreqWin = 16
integer, parameter :: MaxNumOfLamWin  = 16

integer, parameter :: NumOfBasicParam = 2

integer, parameter :: LongInt = 8

type :: type_ray
  double precision x, y, z, vx, vy, vz
end type type_ray


type :: type_photon_packet
  ! Mainly for Monte Carlo.
  type(type_ray) :: ray
  double precision lam, en
  double precision f, Inu
  integer iKap, iSpec
  integer e_count
end type type_photon_packet


type :: type_photon_ray_multi
  ! Mainly for ray-tracing (to make image and spectra)
  type(type_ray) :: ray
  integer iTran, nf
  double precision, dimension(:), allocatable :: lam, f, Inu
  integer, dimension(:), allocatable :: iKap
end type type_photon_ray_multi


type :: type_position_cartesian
  double precision x, y, z
end type type_position_cartesian



type :: type_direction_cartesian
  double precision u, v, w
end type type_direction_cartesian


type :: type_sphere_coor_quat
  double precision costheta, sintheta, cosphi, sinphi
end type type_sphere_coor_quat


type :: type_photon_collector
  ! Use mu = cos(theta) instead of theta to save some calculations.
  integer nlam, iKap0, iKap1, nmu, nr, nphi
  double precision, dimension(:), allocatable :: &
    mu_min, mu_max, r_min, r_max, phi_min, phi_max
  double precision, dimension(:,:,:,:), allocatable :: energy
  integer, dimension(:,:,:,:), allocatable :: counts
end type type_photon_collector


type :: type_spectrum_generic
  integer n
  double precision, dimension(:,:), allocatable :: intervals
  double precision, dimension(:), allocatable :: vals
end type type_spectrum_generic


type :: type_distribution_table
  integer n
  double precision, dimension(:), allocatable :: pvals
end type type_distribution_table


type :: type_simple_integer_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: vals
end type type_simple_integer_list


type :: type_optical_property
  integer n
  double precision, dimension(:), allocatable :: lam
  double precision, dimension(:), allocatable :: ab, sc, g
end type type_optical_property


type :: type_global_material_collection
  integer ntype
  type(type_optical_property), dimension(:), allocatable :: list
  double precision, dimension(:), allocatable :: &
                              Xray_gas_abs, &
                              Xray_gas_sca, &
                              Xray_dus_abs, &
                              Xray_dus_sca
end type type_global_material_collection


type :: type_local_encounter_collection
  integer ntype, nlam
  integer cr_count
  double precision, dimension(:), allocatable :: X
  double precision, dimension(:), allocatable :: ext_tot
  double precision, dimension(:), allocatable :: albedo
  double precision, dimension(:), allocatable :: flux
  integer, dimension(:), allocatable :: phc
#ifdef SAVE_PHOTON_FIELD_DIR
  type(type_direction_cartesian), dimension(:), allocatable :: dir_wei
#endif
end type type_local_encounter_collection


type :: type_accum_extinction
  integer ntype
  integer :: current_iKap=-1, current_cell_id=-1
  double precision, dimension(:), allocatable :: accum
end type type_accum_extinction



type :: type_LUT_Tdust
  ! A look-up table
  ! Given a val, find out the corresponding T
  integer n, m
  double precision, dimension(:), allocatable :: Tds, vals
  double precision, dimension(:,:), allocatable :: table
end type type_LUT_Tdust


type :: type_continuum_lut
  integer :: n=0
  double precision, dimension(:), allocatable :: J
end type type_continuum_lut


type :: type_stellar_params
  integer n
  double precision mass, radius, T, &
    lumi, lumi_Vis, lumi_UV, lumi_Lya, lumi0, lumi_UV0
  double precision :: T_Xray=1D7, E0_Xray=0.1D0, E1_Xray=10D0, lumi_Xray=1D30
  double precision, dimension(:), allocatable :: lam, vals, vals0
end type type_stellar_params



type :: type_montecarlo_config
  double precision eph
  integer(kind=LongInt) nph, icount, nmax_cross, nmax_encounter
  character(len=128) fname_photons, fname_dust, fname_water, fname_star, fname_save_collected
  character(len=128) :: fname_save_seeds = 'random_seeds.dat'
  character(len=128) mc_dir_in, mc_dir_out
  double precision minw, maxw
  double precision :: starpos_r=0D0, starpos_z = 0D0
  logical :: use_blackbody_star=.true. !, savephoton=.false.
  logical :: allow_Xray_scattering = .true.
  logical :: disallow_any_scattering = .false.
  logical :: ph_init_symmetric=.false.
  double precision :: stellar_spectr_UV_rescale_factor = 1D0
  double precision :: refine_UV = 0.01D0, refine_LyA = 0.001D0, refine_Xray = 1D-4
  logical :: collect_photon=.false.
  double precision dist
  double precision :: collect_lam_min=1D0, collect_lam_max=1D6
  double precision :: collect_dmu=0.1D0
  integer :: collect_nmu=3, collect_nr=50, collect_nphi=50
  double precision, dimension(8) :: collect_ang_mins, collect_ang_maxs
  integer :: nlen_lut=1024
  double precision :: TdustMin=1D0, TdustMax=2D3
  logical :: do_fill_blank = .false.
  integer :: fill_blank_threshold = 10
end type type_montecarlo_config


type :: type_energy_level
  character(len=const_len_energy_level) :: name_energy
  !character(len=const_len_qnum*2) :: qnum
  integer id
  double precision :: energy
  double precision :: weight
end type type_energy_level


type :: type_rad_transition
  double precision Eup, Elow, freq, lambda
  double precision Aul, Bul, Blu, beta, J_ave, cooling_rate
  integer iup, ilow
  character(len=const_len_qnum*4) qnum
end type type_rad_transition


type :: type_collisional_transition
  character(len=const_len_molecule) :: name_partner
  double precision dens_partner
  integer n_transition, n_T
  integer, dimension(:), allocatable :: iup, ilow
  double precision, dimension(:), allocatable :: T_coll
  double precision, dimension(:,:), allocatable :: Cul
end type type_collisional_transition


type :: type_rad_set
  integer n_transition
  type(type_rad_transition), dimension(:), allocatable :: list
end type type_rad_set


type :: type_colli_set
  integer n_partner
  type(type_collisional_transition), dimension(:), allocatable :: list
end type type_colli_set


type :: type_molecule_energy_set
  character(len=const_len_molecule) name_molecule, name_surrogate, name_disp
  integer iSpe, iType
  double precision Tkin, density_mol, dv, length_scale, cooling_rate_total
  integer :: n_level = 0
  type(type_energy_level), dimension(:), allocatable :: level_list
  double precision, dimension(:), allocatable :: f_occupation
  type(type_rad_set), allocatable :: rad_data
  type(type_colli_set), allocatable :: colli_data
  double precision :: abundance_factor = 1D0
end type type_molecule_energy_set


type :: type_statistic_equil_params
  integer nitem
  double precision :: RTOL = 1D-2, ATOL = 1D-12
  double precision :: t_max = 1D8, dt_first_step = 1D-6, ratio_tstep = 1.2D0
  real :: max_runtime_allowed = 3.0
  integer n_record
  integer :: &
        NERR, &
        ITOL = 1, &
        ITASK = 1, &
        ISTATE = 1, &
        IOPT = 1, &
        LIW, &
        LRW, &
        MF = 21
  double precision, dimension(:), allocatable :: RWORK
  integer, dimension(:), allocatable :: IWORK
  logical is_good
end type type_statistic_equil_params

type :: type_mole_f_occ
  integer nlevels
  double precision, dimension(:), allocatable :: vals
end type type_mole_f_occ


type :: type_mole_exc_conf
  character(len=128) :: dirname_mol_data=''
  character(len=128) :: fname_mol_data=''
  character(len=128) :: fname_parti_data=''
  character(len=16) :: mole_name='', mole_name_surrogate='', mole_name_disp = ''
  character(len=128) :: dir_save_image=''
  character(len=8) :: line_database='lamda'
  logical :: save_contributions = .false.
  character(len=128) :: fname_contributions
  integer :: fU_save_contri, inu_sav=0
  double precision :: dlength_contri = 0.1D0
  integer nfreq_window
  double precision, dimension(MaxNumOfFreqWin) :: freq_mins, freq_maxs
  integer nlam_window
  double precision, dimension(MaxNumOfLamWin)  :: lam_mins, lam_maxs
  double precision abundance_factor
  double precision :: E_min = 0D0, E_max = 5D3
  double precision :: Aul_min = 0D0, Aul_max = 1D99
  double precision :: min_flux=0D0
  double precision :: VeloKepler = 0D0, VeloTurb = 0D0, VeloWidth = 1D0
  logical :: useLTE = .true.
  logical :: save_spectrum_only = .false.
  logical :: subtract_cont_tau = .true.
  logical :: adjust_yup_ylow_nonLTE = .false.
  double precision :: n_critical_CGS
  integer :: solve_method = 1
  !
  logical :: turn_off_dust = .false.
  !
  double precision :: maxx=0D0, maxy=0D0
  integer nf, nlam, nth, nx, ny
  double precision dist
  double precision, dimension(16) :: view_thetas
  !
end type type_mole_exc_conf


type :: type_molecule_exc
  type(type_mole_exc_conf) :: conf
  type(type_molecule_energy_set), pointer :: p => null()
  integer nlevel_keep, ntran_keep
  integer, dimension(:), allocatable :: ilv_keep, ilv_reverse
  integer, dimension(:), allocatable :: itr_keep, itr_reverse
end type type_molecule_exc


type :: type_fits_par
  character(len=256) :: filename
  integer stat, fU, blocksize, bitpix, naxis
  integer, dimension(3) :: naxes
  integer i, j, group, fpixel, nelements, decimals
  integer pcount, gcount
  logical simple, extend
  character(len=32) :: extname
  character(len=32) :: author, user
end type type_fits_par


type :: type_cell_rz_phy_basic
  integer(kind=LongInt) ab_count_dust, ab_count_water, &
                        sc_count_dust, sc_count_HI
  !
  integer ndustcompo
  !
  double precision :: &
    volume, surf_area, area_T, area_B, area_I, area_O, &
    Tgas, &
    Tdust, &
    grand_gas_abundance, &
    !
    n_gas, &
    !
    mgas_cell, &
    !
    Tdusts(MaxNumOfDustComponents), &
    rho_dusts(MaxNumOfDustComponents), &
    n_dusts(MaxNumOfDustComponents), &
    mp_dusts(MaxNumOfDustComponents), &
    mdusts_cell(MaxNumOfDustComponents), &
    dust2gas_scale_height_ratio(MaxNumOfDustComponents) = 1D0, &
    !
    en_exchange_per_vol(MaxNumOfDustComponents), &
    en_exchange(MaxNumOfDustComponents), &
    en_exchange_tot, &
    dEmit_dTd(MaxNumOfDustComponents), &
    !
    abso_wei(MaxNumOfDustComponents), &
    !
    sig_dusts(MaxNumOfDustComponents), &
    sig_dusts0(MaxNumOfDustComponents), &
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
    UV_G0_factor_background, &
    !
    ab_en_water, &
    !
    phflux_Lya, &
    !
    flux_UV_star_unatten, &
    flux_Lya_star_unatten, &
    !flux_Vis_star_unatten, &
    !
    G0_UV_toISM, &
    G0_UV_toStar, &
    G0_Lya_atten, &
    G0_UV_H2phd, &
    G0_UV_toStar_photoDesorb, &
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
    zeta_Xray_H2, &
    !
    R_H2_form_rate_coeff, &
    R_H2_form_rate, &
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
    n_mol_on_grain, &
    dust_depletion, &
    PAH_abundance, &
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
    ambipolar_f, &
    ion_charge, &
    Neufeld_G, &
    Neufeld_dv_dz, &
    !
    t_final, &
    tmax_this, &
    !
    X_H2, X_HI, X_CI, X_CII, X_OI, X_CO, &
    X_H2O, X_OH, X_E, X_Hplus, X_gH, X_Heplus, &
    X_NII, X_FeII, X_SiII, &
    flux_tot, flux_Xray, flux_UV, flux_Lya, &
    flux_Vis, flux_NIR, flux_MIR, flux_FIR, &
    !
    dir_tot_r=0D0, dir_tot_z =0D0, dir_Xray_r=0D0, dir_Xray_z=0D0, &
    dir_UV_r =0D0, dir_UV_z  =0D0, dir_Lya_r =0D0, dir_Lya_z =0D0, &
    dir_Vis_r=0D0, dir_Vis_z =0D0, dir_NIR_r =0D0, dir_NIR_z =0D0, &
    dir_MIR_r=0D0, dir_MIR_z =0D0, dir_FIR_r =0D0, dir_FIR_z =0D0, &
    aniso_tot=0D0, aniso_Xray=0D0, aniso_UV  =0D0, aniso_Lya =0D0, &
    aniso_Vis=0D0, aniso_NIR =0D0, aniso_MIR =0D0, aniso_FIR =0D0, &
    !
    pressure_thermal, gravity_z, gravity_acc_z
  !
end type type_cell_rz_phy_basic


type :: type_dust_MRN
  double precision :: rmin, rmax, n
  double precision :: rav, r2av, r3av
end type type_dust_MRN


type :: type_Andrews_disk
  logical :: useNumDens = .true.
  double precision :: particlemass = 1.4D0 * 1.67262158D-24
  double precision :: Md=0.00D0 ! Disk mass in Msun
  double precision :: rin=0.5D0, rout=200D0
  double precision :: rc=200D0  ! Disk outer boundary
  double precision :: hc=50D0   ! Scale height at outer boundary
  double precision :: rc_for_hc=0D0  ! For calculating hc; if zero then set equal to rc
  double precision :: gam=1D0   ! Power index for surface density
  double precision :: psi=1D0   ! Power index for scale height
  !
  double precision :: r0_in_exp = 0D0 ! Exponential taper inward this radius
  double precision :: rs_in_exp = 1D99 ! Scale length of the exponential taper
  double precision :: p_in_exp  = 1D0 ! Power index
  double precision :: f_in_exp  = 1D0 ! Prefactor
  double precision :: r0_out_exp = 1D99 ! Exponential taper outward this radius
  double precision :: rs_out_exp = 1D99 ! Scale length of the exponential taper
  double precision :: p_out_exp  = 1D0 ! Power index
  double precision :: f_out_exp  = 1D0 ! Prefactor
  !
  double precision :: r0_in_change = 0D0 ! Change vertical scale inward this radius
  double precision :: f_in_change = 1D0 ! Factor of vertical scale change
  double precision :: r0_out_change = 1D99 ! Change vertical scale outward this radius
  double precision :: f_out_change = 1D0 ! Factor of vertical scale change
  !
  double precision :: r_in_flatten = 0D0 ! Flatten the surface density inside this radius
  !
end type type_Andrews_disk


type :: type_a_dust_component
  integer itype
  double precision pmass_CGS
  double precision :: d2g_ratio = 0.01D0
  type(type_dust_MRN) :: mrn
  type(type_Andrews_disk) :: andrews
end type type_a_dust_component


type :: type_heating_cooling_rates_list
  double precision :: &
    hc_net_rate = 0D0, &
    heating_photoelectric_small_grain_rate = 0D0, &
    heating_formation_H2_rate = 0D0, &
    heating_cosmic_ray_rate = 0D0, &
    heating_vibrational_H2_rate = 0D0, &
    heating_ionization_CI_rate = 0D0, &
    heating_photodissociation_H2_rate = 0D0, &
    heating_photodissociation_H2O_rate = 0D0, &
    heating_photodissociation_OH_rate = 0D0, &
    heating_Xray_Bethell_rate = 0D0, &
    heating_viscosity_rate = 0D0, &
    heating_chem = 0D0, &
    cooling_photoelectric_small_grain_rate = 0D0, &
    cooling_vibrational_H2_rate = 0D0, &
    cooling_gas_grain_collision_rate = 0D0, &
    cooling_OI_rate = 0D0, &
    cooling_CII_rate = 0D0, &
    cooling_Neufeld_H2O_rate_rot = 0D0, &
    cooling_Neufeld_H2O_rate_vib = 0D0, &
    cooling_Neufeld_CO_rate_rot = 0D0, &
    cooling_Neufeld_CO_rate_vib = 0D0, &
    cooling_Neufeld_H2_rot_rate = 0D0, &
    cooling_LymanAlpha_rate = 0D0, &
    cooling_free_bound_rate = 0D0, &
    cooling_free_free_rate = 0D0, &
    cooling_NII_rate = 0D0, &
    cooling_SiII_rate = 0D0, &
    cooling_FeII_rate = 0D0, &
    cooling_OH_rot_rate = 0D0
end type type_heating_cooling_rates_list


type :: type_cell_ptr
  type(type_cell), pointer :: p
end type type_cell_ptr

type :: type_leaves
  integer :: nlen = 0
  type(type_cell_ptr), dimension(:), allocatable :: list
end type type_leaves

type :: type_neighbor
  integer :: n = 0
  integer, dimension(:), allocatable :: idx
  !double precision, dimension(:), allocatable :: fra
  !double precision :: fra_tot = 0D0
end type type_neighbor


type :: type_cell
  double precision :: xmin=0D0, xmax=0D0, ymin=0D0, ymax=0D0
  double precision, dimension(NumOfBasicParam) :: val
  logical :: using = .false., converged = .false.
  integer :: id = -1
  integer :: order=0, nChildren=0, nOffspring=0, nleaves=0
  type(type_cell), pointer :: parent => null()
  type(type_cell_ptr), dimension(:), allocatable :: children
  type(type_neighbor), pointer :: &
        inner => null(), outer => null(), &
        below => null(), above => null(), &
        around => null()
  type(type_cell_rz_phy_basic), pointer :: par => null()
  type(type_heating_cooling_rates_list), allocatable :: h_c_rates
  double precision, dimension(:), allocatable :: abundances
  double precision, dimension(:), allocatable :: col_den_toStar, col_den_toISM
  type(type_local_encounter_collection), allocatable :: optical
  type(type_continuum_lut), allocatable :: cont_lut
  type(type_mole_f_occ), allocatable :: focc
  integer :: iIter = 0
  integer :: quality = 0
end type type_cell


end module data_struct

