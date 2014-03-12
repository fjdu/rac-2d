## Install

Go to the ```src``` directory.  Inside it there is a ```makefile```.  You may
edit it for your own needs (but you don't have to).  Then run
```bash
    make
```
and an executable file with default name ```a.out``` will be generated in the
same directory.

### Dependency

    1. Gfortran higher than 4.6.2 or Intel Fortran higher than 12.1.5.
    2. The cfitsio library.


## Run the code

There are a few input data files that are needed for the code to run.

The following files are compulsary:

    1. Configuration file.
    2. Chemical network.
    3. Initial chemical composition.
    4. Dust optical properties.

The following files are optional:

    1. Density structure.
    2. Enthalpy of formation of species.
    3. Molecular ransition data.
    4. Stellar spectrum.
    5. Points to output the intermediate steps of chemcial evolution.
    6. Species to output the intermediate steps of chemcial evolution.
    7. Species to check for grid refinement.

By default all these files are in the ```inp``` directory, though some of them
do not have to.  Go to this directory.  Edit the file ```configure.dat```.  It
has nearly 200 entries.  Some of them are for setting up the physics and
chemistry of the model, some are for setting up the running environment, while
others are switches telling the code whether or not it should execute some
specific tasks.  Details for editing the configure file are included below.

After you have get the configre file ready, and have all the needed files in
place, then in a terminal you can go to the directory on top of ```inp```, and
type in
```
    ./src/a.out ./inp/configure.dat
```
to start running the code.

With the template files that are already there the configuration file should
let the code run without any modification needed.

## Contents of configure.dat

The configuration file is in the Fortran _namelist_ format, so you may want to
set the language type for syntax highlighting of your editor to Fortran.

At the end of the configuration file you can write down any notes as you want.
They will not be read in by the code.

```fortran
! All comments should be preceded by "!".
! Inline comments should be separated from the content by at least one blank
! space.
&grid_configure
  grid_config%rmin = 1D0   ! Grid inner boundary
  grid_config%rmax = 10D0  ! Grid outer boundary
  grid_config%zmin = 0D0   ! Grid lower boundary
  grid_config%zmax = 10D0  ! Grid upper boundary
  grid_config%dr0  = 1D-3  ! Width of the first r step
  grid_config%columnwise = .true.
  grid_config%ncol = 50
  grid_config%use_data_file_input = .false.
  grid_config%data_dir = './inp/'
  grid_config%data_filename = 'RADMC_density_temperature.dat'
  grid_config%analytical_to_use = 'Andrews'
  grid_config%interpolation_method = 'spline'
  grid_config%max_ratio_to_be_uniform = 1.2D0
  grid_config%density_scale     = 14D0
  grid_config%density_log_range = 6D0
  grid_config%max_val_considered = 1d19
  grid_config%min_val_considered = 1d4
  grid_config%very_small_len = 1D-6
  grid_config%smallest_cell_size = 1D-2
  grid_config%largest_cell_size  = 1D0
  grid_config%small_len_frac = 1D-3
/
&chemistry_configure
  chemsol_params%allow_stop_before_t_max     = .false.
  chemsol_params%dt_first_step               = 1D-8
  chemsol_params%t_max                       = 1D6
  chemsol_params%ratio_tstep                 = 1.1D0
  chemsol_params%max_runtime_allowed         = 30.0
  chemsol_params%RTOL                        = 1D-5
  chemsol_params%ATOL                        = 1D-30
  chemsol_params%mxstep_per_interval         = 6000
  chemsol_params%chem_files_dir              = './inp/'
  chemsol_params%filename_chemical_network   = 'rate06_dipole_reformated_again_withgrain.dat'
  chemsol_params%filename_initial_abundances = 'initial_condition_Garrod08_mod_waterice.dat'
  chemsol_params%filename_species_enthalpy   = 'Species_enthalpy.dat'
  chemsol_params%H2_form_use_moeq            = .false.
  chemsol_params%neutralize                  = .true.
  chemsol_params%flag_chem_evol_save         = .false.
/
&heating_cooling_configure
  ! When use_analytical_CII_OI is true, the two files will not be used.
  heating_cooling_config%use_analytical_CII_OI   = .true.
  heating_cooling_config%dir_transition_rates    = './inp/'
  heating_cooling_config%filename_Cplus          = 'C+.dat'
  heating_cooling_config%filename_OI             = 'Oatom.dat'
  heating_cooling_config%use_mygasgraincooling   = .true.
/
&montecarlo_configure
  mc_conf%nph                   = 1000000     ! Divide total star luminosity into this number.
  mc_conf%nmax_cross            = 1999999999  ! Max num of cell crossing before any ab or sc
  mc_conf%nmax_encounter        = 1999999999  ! Max num of absor and scat events
  mc_conf%min_ang               = 0D0  ! Angle range of rays
  mc_conf%max_ang               = 15D0
  mc_conf%refine_UV             = 0.1D0
  mc_conf%refine_LyA            = 0.01D0
  mc_conf%savephoton            = .false.
  mc_conf%mc_dir_in             = './inp/'
  mc_conf%mc_dir_out            = 'mc/'
  mc_conf%fname_photons         = 'escaped_photons.dat'
  mc_conf%fname_water           = 'H2O.photoxs'
  mc_conf%fname_star            = 'tw_hya_spec_combined.dat'
  mc_conf%use_blackbody_star    = .false.
/
&dustmix_configure
  dustmix_info%nmixture = 1  ! Number of mixtures you want to make
  dustmix_info%lam_min = 0.09D0 ! Minimum wavelength (micron) to be considered
  dustmix_info%lam_max = 1D3    ! Maximum ...
  !
  dustmix_info%mix(1)%id = 1  ! Mixture 1
  dustmix_info%mix(1)%nrawdust = 3  ! Number of raw material for mixing
  dustmix_info%mix(1)%rho = 3D0 ! Dust material density in g cm-3
  dustmix_info%mix(1)%dir          = './inp/'
  dustmix_info%mix(1)%filenames(1) = 'suvSil_81'  ! Filename of raw material 1
  dustmix_info%mix(1)%filenames(2) = 'Gra_81'
  dustmix_info%mix(1)%filenames(3) = 'SiC_81'
  dustmix_info%mix(1)%weights(1) = 0.7D0  ! Weight of raw material 1
  dustmix_info%mix(1)%weights(2) = 0.2D0
  dustmix_info%mix(1)%weights(3) = 0.1D0
  !
  !dustmix_info%mix(2)%id = 2  ! Mixture 2
  !dustmix_info%mix(2)%nrawdust = 3
  !dustmix_info%mix(2)%dir = '/Users/fdu/Dropbox/work/protoplanetary_disk/inp/'
  !dustmix_info%mix(2)%filenames(1) = 'suvSil_81'
  !dustmix_info%mix(2)%filenames(2) = 'Gra_81'
  !dustmix_info%mix(2)%filenames(3) = 'SiC_81'
  !dustmix_info%mix(2)%weights(1) = 0.7D0
  !dustmix_info%mix(2)%weights(2) = 0.2D0
  !dustmix_info%mix(2)%weights(3) = 0.1D0
/
&disk_configure
  a_disk%star_luminosity_in_Lsun   = 1D0 ! Not used for montecarlo; only for xray lumi.
  a_disk%star_mass_in_Msun         = 0.6D0
  a_disk%star_radius_in_Rsun       = 1D0
  a_disk%star_temperature          = 4000D0
  a_disk%starpos_r                 = 0D0
  a_disk%starpos_z                 = 0D0
  !
  a_disk%use_fixed_alpha_visc      = .false.
  a_disk%allow_gas_dust_en_exch    = .false.
  a_disk%base_alpha                = 0.01D0
  !
  a_disk%waterShieldWithRadTran    = .false.
  !
  !a_disk%ratio_uv2total            = 1D-2  ! Not used because this ratio is calcualted from the stellar spectrum.
  !a_disk%ratio_lyman2uv            = 0.1D0 ! Not used.
  a_disk%ratio_xray2total          = 1D-3  ! In use, for heating.
  !a_disk%geometric_factor_UV       = 1D-2  ! Not used.
  a_disk%geometric_factor_Xray     = 1D0  ! In use.
  !
  a_disk%backup_src                = .true. !
  a_disk%backup_src_cmd            = 'find src/*.f90 src/*.f src/makefile inp/*dat | cpio -pdm '
  !
  a_disk%andrews_gas%useNumDens    = .true.  ! Gas distribution
  a_disk%andrews_gas%Md            = 5D-2  ! Total mass in Msun
  a_disk%andrews_gas%rin           = 1D0    ! Inner edge
  a_disk%andrews_gas%rout          = 140D0  ! Outer edge
  a_disk%andrews_gas%rc            = 100D0  ! Characteristic radius
  a_disk%andrews_gas%hc            = 10D0  ! Scale height at boundary
  a_disk%andrews_gas%gam           = 1.0D0  ! Power index
  a_disk%andrews_gas%psi           = 1.0D0  ! 
  !
  a_disk%ndustcompo                = 2  ! Number of dust components
  !
  a_disk%dustcompo(1)%itype        = 1  ! Dust mixture type to use
  a_disk%dustcompo(1)%mrn%rmin     = 1D0  ! Min radius
  a_disk%dustcompo(1)%mrn%rmax     = 100D0  ! Max radius
  a_disk%dustcompo(1)%mrn%n        = 3.5D0  ! Power index
  a_disk%dustcompo(1)%andrews%useNumDens = .false.  ! Use mass density insteady of number density
  a_disk%dustcompo(1)%andrews%Md         = 5D-6
  a_disk%dustcompo(1)%andrews%rin        = 1D0
  a_disk%dustcompo(1)%andrews%rout       = 140D0
  a_disk%dustcompo(1)%andrews%rc         = 100D0
  a_disk%dustcompo(1)%andrews%hc         = 10D0
  a_disk%dustcompo(1)%andrews%gam        = 1.0D0
  a_disk%dustcompo(1)%andrews%psi        = 1.0D0
  !
  a_disk%dustcompo(2)%itype        = 1  ! Dust mixture type to use
  a_disk%dustcompo(2)%mrn%rmin     = 0.01D0  ! Min radius
  a_disk%dustcompo(2)%mrn%rmax     = 1D0  ! Max radius
  a_disk%dustcompo(2)%mrn%n        = 3.5D0  ! Power index
  a_disk%dustcompo(2)%andrews%useNumDens = .false.  ! Use mass density insteady of number density
  a_disk%dustcompo(2)%andrews%Md         = 1D-8
  a_disk%dustcompo(2)%andrews%rin        = 1D0
  a_disk%dustcompo(2)%andrews%rout       = 140D0
  a_disk%dustcompo(2)%andrews%rc         = 100D0
  a_disk%dustcompo(2)%andrews%hc         = 10D0
  a_disk%dustcompo(2)%andrews%gam        = 1.0D0
  a_disk%dustcompo(2)%andrews%psi        = 1.0D0
/
&mole_line_configure
  mole_line_conf%dirname_mol_data  = './transitions/'
  mole_line_conf%fname_mol_data    = 'oh2o@rovib.dat'
  mole_line_conf%nfreq_window      = 1
  mole_line_conf%freq_mins         = 375.0D9
  mole_line_conf%freq_maxs         = 385.0D9
  mole_line_conf%E_min             = 1D2
  mole_line_conf%E_max             = 5D3
  mole_line_conf%abundance_factor  = 2D-3 ! H2-18O
  mole_line_conf%useLTE            = .true.
  mole_line_conf%nf                = 100
  mole_line_conf%nth               = 4
  mole_line_conf%nx                = 101
  mole_line_conf%ny                = 101
  mole_line_conf%dist              = 55.0
/
&cell_configure
  cell_params_ini%omega_albedo              = 0.5D0  ! Dust albedo, only for chemistry
  cell_params_ini%UV_G0_factor_background   = 1.0D0  ! ISM UV
  cell_params_ini%zeta_cosmicray_H2         = 1.36D-17  ! Cosmic ray intensity
  cell_params_ini%GrainMaterialDensity_CGS  = 2D0  ! Density of dust material
  cell_params_ini%MeanMolWeight             = 1.4D0
  cell_params_ini%alpha_viscosity           = 0.01D0
/
&analyse_configure
  ! Do chemical analysis for some species at some locations
  a_disk_ana_params%do_analyse                      = .true.
  a_disk_ana_params%analyse_points_inp_dir          = './inp/'
  a_disk_ana_params%file_list_analyse_points        = 'points_to_analyse.dat'
  a_disk_ana_params%file_list_analyse_species       = 'Species_to_analyse.dat'
  a_disk_ana_params%file_analyse_res_ele            = 'elemental_reservoir.dat'
  a_disk_ana_params%file_analyse_res_contri         = 'contributions.dat'
/
&iteration_configure
  a_disk_iter_params%n_iter                         = 10
  a_disk_iter_params%nlocal_iter                    = 4
  a_disk_iter_params%rtol_abun                      = 0.2D0  ! For checking convergency of a cell
  a_disk_iter_params%atol_abun                      = 1D-12
  a_disk_iter_params%converged_cell_percentage_stop = 0.95  ! Stop if so many cells have converged.
  a_disk_iter_params%filename_list_check_refine     = 'species_check_refine.dat'  ! Do refinement based on species in this file.
  a_disk_iter_params%threshold_ratio_refine         = 10D0  ! Do refinement if the gradient (ratio) is so large
  a_disk_iter_params%nMax_refine                    = -100  ! Max times of refinments
  a_disk_iter_params%redo_montecarlo                = .true.  ! Redo Monte Carlo after each full chemcial run.
  a_disk_iter_params%flag_save_rates                = .false.  ! Whether save the calculated rates.
  a_disk_iter_params%flag_shortcut_ini              = .false.  ! Whether use short-cut initial condition
  a_disk_iter_params%do_line_transfer               = .false.
  a_disk_iter_params%iter_files_dir                 = './res/20140311_nowaterRad_lessDust_1/'
/
!-Notes------------------------------------------------------------------
! 1. Based on ./res/20140218_a/ and ./res/20140131_a/.
! 2. With many iterations.
! 3. Water shielding not included in the radiative transfer.
! 4. Redo.
```
