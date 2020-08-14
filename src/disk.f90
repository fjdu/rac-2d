module disk

use data_struct
use grid
use chemistry
use heating_cooling
use montecarlo
use load_Draine_dusts
use data_dump
use vertical_structure
use my_timer
use statistic_equilibrium


implicit none

type :: type_disk_basic_params
  double precision star_mass_in_Msun, star_radius_in_Rsun, star_temperature
  double precision :: distance = 0D0
  double precision :: T_Xray=1D7, E0_Xray=0.1D0, E1_Xray=10D0, lumi_Xray=1D30
  double precision :: starpos_r = 0D0, starpos_z = 0D0
  !
  double precision disk_mass_in_Msun
  !
  integer :: ngascompo=1, ndustcompo=1
  !type(type_Andrews_disk), dimension(MaxNumOfGasComponents) :: andrews_gas
  type(type_Andrews_disk) :: andrews_gas
  type(type_a_dust_component), dimension(MaxNumOfDustComponents) :: dustcompo
  !
  logical :: use_fixed_alpha_visc=.true.
  logical :: allow_gas_dust_en_exch=.false.
  double precision :: base_alpha = 0.01D0
  !
  logical :: waterShieldWithRadTran = .true.
  logical :: Tdust_iter_tandem = .false.
end type type_disk_basic_params


type :: type_disk_iter_params
  integer :: n_iter=128, n_iter_used = 0
  integer :: nlocal_iter = 2
  !
  double precision :: rtol_T = 0.1D0,    atol_T = 2D0
  double precision :: rtol_abun = 0.2D0, atol_abun = 1D-12
  !
  double precision :: Tgas_crazy = 3D4
  double precision :: n_gas_thrsh_noTEvol = 1D15
  double precision :: minimum_Tdust = 5D0
  !
  logical flag_converged
  integer n_cell_converged
  real converged_cell_percentage_stop
  !
  integer :: max_num_of_cells = 10000
  !
  logical :: rerun_whole = .false.
  logical :: rerun_whole_noMonteCarlo = .true.
  logical :: redo_something = .false.
  double precision :: single_p_x=0D0, single_p_y=0D0
  !
  double precision :: deplete_oxygen_carbon_tstart = 0D0
  logical :: deplete_oxygen_carbon = .false.
  logical :: deplete_nitrogen = .false.
  logical :: deplete_nitrogen_as_carbon = .false.
  character(len=16) :: deplete_oxygen_carbon_method = ''
  character(len=16) :: deplete_oxygen_method = ''
  character(len=16) :: deplete_carbon_method = ''
  double precision a_O, b_O, gam_O, vfac_O
  double precision a_C, b_C, gam_C, vfac_C
  double precision :: f_depl_O = 1D0, f_depl_C = 1D0
  integer nvfacs_O, nvfacs_C
  double precision,dimension(16) :: rmins_O, rmaxs_O, vfacs_O, &
                                    rmins_C, rmaxs_C, vfacs_C
  !
  double precision :: gval_O=1D-4, gval_C=1D-4
  double precision :: tads_O=1D2, tads_C=1D2, tsed_O=1D5, tsed_C=1D5, &
                      r0_O=50D0, r0_C=50D0, &
                      k_O=1D0, k_C=1D0, &
                      p_O=1D0, p_C=1D0, &
                      rin_O=0D0, rin_C=0D0, &
                      fin_O=1D0, fin_C=1D0, &
                      rout_O=1D9, rout_C=1D9, &
                      fout_O=1D0, fout_C=1D0
  double precision :: tanh_r_O = 0D0, tanh_scale_O = 1D2, tanh_minval_O = 0.999D0, tanh_maxval_O = 1D0, &
                      tanh_r_C = 0D0, tanh_scale_C = 1D2, tanh_minval_C = 0.999D0, tanh_maxval_C = 1D0
  double precision :: tanh_OC_enhance_max_O = 10D0, tanh_OC_enhance_max_C = 10D0, tanh_OC_enhance_max = 1D99
  double precision :: O_to_C_ISM = 2.285714D0  ! = 3.2/1.4
  double precision :: C_to_O_ratio = 0D0
  double precision :: dep_zscale = 0D0
  !
  logical :: use_fixed_tmax = .false.
  double precision :: nOrbit_tmax = 1D5
  double precision :: minDust2GasNumRatioAllowed = 1D-15
  logical :: vertical_structure_fix_grid = .true.
  logical :: vertical_structure_fix_dust = .false.
  logical :: calc_Av_toStar_from_Ncol = .false.
  logical :: calc_zetaXray_from_Ncol  = .false.
  double precision :: dust2gas_mass_ratio_deflt = 0.01D0  ! Will be obsolete; use d2g_ratio instead
  logical :: rescale_ngas_2_rhodust = .false.
  logical :: do_gas_vertical_simple = .false.
  logical :: do_simple_chem_before_mc = .false.
  integer :: do_vertical_every = 3
  integer :: nVertIterTdust = 6
  logical :: do_vertical_struct = .false.
  logical :: do_vertical_with_Tdust = .false.
  logical :: redo_montecarlo = .true.
  logical :: flag_save_rates = .false.
  !
  integer :: nSpecies_check_refine = 0
  integer :: ncell_refine = 0, count_refine = 0
  integer :: nMax_refine = 2
  double precision :: threshold_ratio_refine = 10D0
  character(len=128) filename_list_check_refine
  !
  character(len=128) iter_files_dir
  character(len=256) filename_exe
  logical            :: backup_src = .true.
  character(len=512) :: backup_src_cmd = &
    'find *.f90 *.f *.py makefile | cpio -pdm --insecure '
  !
  logical :: do_line_transfer=.false., do_continuum_transfer=.false.
  !
  character(len=128) :: dump_common_dir=''
  character(len=128) :: dump_sub_dir_in='', dump_sub_dir_out=''
  character(len=128) :: dump_filename_optical='', dump_filename_chemical='', &
                        dump_filename_physical='', dump_filename_physical_aux='', &
                        dump_filename_grid=''
  logical :: use_backup_optical_data  = .false.
  logical :: use_backup_chemical_data = .false.
  logical :: use_backup_physical_data = .false.
  logical :: use_backup_grid_data  = .false.
  !
end type type_disk_iter_params


type :: type_ana_params
  logical :: do_analyse = .false.
  integer ana_i_incr
  character(len=128) analyse_points_inp_dir, analyse_out_dir
  character(len=128) file_list_analyse_points, file_list_analyse_species, &
                     file_analyse_res_ele, file_analyse_res_contri
  type(type_cell_rz_phy_basic) chempar
end type type_ana_params


type :: type_disk_iter_storage
  double precision, dimension(:), allocatable :: T_s
  double precision, dimension(:,:), allocatable :: abundances
end type type_disk_iter_storage


type :: type_book_keeping
  integer fU
  character(len=128) dir, filename_log
end type type_book_keeping

! For logging
type(type_book_keeping) a_book_keeping


type(date_time), private :: a_date_time

! Basic config params of the disk model
type(type_disk_basic_params) a_disk

! Iteration params for the disk model
type(type_disk_iter_params) a_disk_iter_params

! Storage for the disk model
type(type_disk_iter_storage) a_iter_stor

! Initialization params for each cell, which are the same for all the cells.
type(type_cell_rz_phy_basic) cell_params_ini

! Params for doing analysis
type(type_ana_params) a_disk_ana_params

! Point (location) list and species list for analysis
type(type_simple_integer_list) :: ana_ptlist, ana_splist

! Index list of the cells that are being calculated
integer, dimension(:), allocatable :: calculating_cells
integer n_calculating_cells, n_calculating_cells_max

! Filename for saving the results for each iteration; the name will be changed
! for different iterations.
character(len=128) :: filename_save_results
integer fU_save_results

! Perez-Becker 2011
double precision, parameter, private :: beta_ion_neutral_colli = 2D-9

! Index of species that are to be used for check whether cell refinement is
! needed.
integer, dimension(:), allocatable, private :: idx_Species_check_refine
double precision, dimension(:), allocatable, private :: &
    thr_Species_check_refine

! For displaying some text to the screen
character(len=256) str_disp

character(len=256), private :: dump_dir_in, dump_dir_out

! Field length for certain output
integer, parameter :: len_item=14

! Namelists for reading config params
namelist /disk_configure/ &
  a_disk

namelist /cell_configure/ &
  cell_params_ini

namelist /iteration_configure/ &
  a_disk_iter_params

namelist /analyse_configure/ &
  a_disk_ana_params


contains


subroutine disk_iteration
  integer i, i0, i_count, l_count, ii, iVertIter
  !
  dump_dir_in = combine_dir_filename(a_disk_iter_params%dump_common_dir, &
                                  a_disk_iter_params%dump_sub_dir_in)
  dump_dir_out = combine_dir_filename(a_disk_iter_params%dump_common_dir, &
                                  a_disk_iter_params%dump_sub_dir_out)
  call my_mkdir(dump_dir_out)
  write(*, '(2A)') 'Binary files are saved in: ', dump_dir_out
  !
  call disk_iteration_prepare
  !
  call montecarlo_prep
  !
  call save_post_config_params
  !
  call do_vertical_struct_with_Tdust
  !
  ! Now start the major big loop.
  !
  a_disk_iter_params%count_refine = 0
  !
  write(str_disp, '("! ", A)') "Iteration begins."
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '(A)') '! Current time: ' // &
    trim(a_date_time%date_time_str())
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  do ii = 1, a_disk_iter_params%n_iter
    !
    a_disk_iter_params%n_iter_used = ii
    !
    ! Load or save grid information.
    call do_grid_stuff(ii, overwrite=.true.)
    !
    ! Do the optical stuff AFTER loading the grid.
    write(*, '(A)') 'Doing optical stuff.'
    call do_optical_stuff(ii, overwrite=.true.)
    !
    ! Write header to the file
    call disk_save_results_pre
    !
    call prepare_cont_lut_for_line_cooling
    !
    call do_chemical_stuff(ii)
    !
    write(str_disp, '("! ", A, I4, A)') "Global iteration ", ii, " finished."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // &
        trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! At this point all the layers have been walked through,
    ! so check whether the global iteration has converged.
    call check_convergency_whole_disk
    !
    do i=1, leaves%nlen
      a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
      a_iter_stor%abundances(:,i) = &
        leaves%list(i)%p%abundances(chem_idx_some_spe%idx)
    end do
    !
    write(fU_save_results, '(A, L6)') &
      '! flag_converged = ', a_disk_iter_params%flag_converged
    write(fU_save_results, '(A)') &
      '! Finish saving ' // trim(filename_save_results)
    write(fU_save_results, '(A)') '! at ' // trim(a_date_time%date_time_str())
    flush(fU_save_results)
    close(fU_save_results)
    !
    if (a_disk_iter_params%flag_converged) then
      ! Converged.  Finish iteration.
      exit
    end if
    !
    if (a_disk_iter_params%redo_montecarlo) then
      ! If you don't redo Monte Carlo, there is no point to adjust the vertical
      ! structure.
      if (a_disk_iter_params%do_vertical_struct .and. &
          (ii .lt. a_disk_iter_params%n_iter) .and. &
          (mod(ii, a_disk_iter_params%do_vertical_every) .eq. &
            (a_disk_iter_params%do_vertical_every-1))) then
        !
        write(str_disp, '(A)') '! Adjusting the vertical structure.'
        call display_string_both(str_disp, a_book_keeping%fU)
        !
        if (a_disk_iter_params%do_gas_vertical_simple) then
          call vertical_pressure_gravity_balance_simple
        else
          if (a_disk_iter_params%vertical_structure_fix_grid) then
            call vertical_pressure_gravity_balance_alt(a_disk%star_mass_in_Msun, &
                fix_dust_struct=a_disk_iter_params%vertical_structure_fix_dust, &
                disk_gas_mass_preset = a_disk%disk_mass_in_Msun)
          else
            call vertical_pressure_gravity_balance
          end if
        end if
        !
        call remake_index
        !
        call post_vertical_structure_adj
        !
        write(str_disp, '(A)') '! Current time: ' // &
                               trim(a_date_time%date_time_str())
        call display_string_both(str_disp, a_book_keeping%fU)
        !
        !cycle
      end if
    end if
    !
    if (a_disk_iter_params%count_refine .gt. &
             a_disk_iter_params%nMax_refine) then
      write(str_disp, '(A, I4, " > ", I4)') &
        '! Will not refine anymore. count_refine: ', &
        a_disk_iter_params%count_refine, a_disk_iter_params%nMax_refine
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    end if
    !
    if (ii .eq. a_disk_iter_params%n_iter) then
      ! Do no refinement after the final iteration
      exit
    end if
    write(*, '(/A)') 'Doing refinements where necessary.'
    !
    if (leaves%nlen .gt. a_disk_iter_params%max_num_of_cells) then
      cycle
    end if
    !
    call do_refine
    !
    if (a_disk_iter_params%ncell_refine .ge. 1) then
      !
      a_disk_iter_params%count_refine = a_disk_iter_params%count_refine + 1
      a_disk_iter_params%flag_converged = .false.
      !
      write(*, '(I5, " out of ", I5, " cells are refined.", /)') &
        a_disk_iter_params%ncell_refine, leaves%nlen
      !
      call post_vertical_structure_adj
      !
      if (allocated(calculating_cells)) then
        deallocate(calculating_cells)
      end if
      n_calculating_cells_max = leaves%nlen
      allocate(calculating_cells(n_calculating_cells_max))
      !
      write(str_disp, '("!", A, 2X, I8)') 'New number of cells (leaf):', leaves%nlen
      call display_string_both(str_disp, a_book_keeping%fU)
      write(str_disp, '("!", A, 2X, I8)') 'New number of cells (total):', &
          root%nOffspring
      call display_string_both(str_disp, a_book_keeping%fU)
    else
      write(str_disp, '("! ", A)') "No further refinement needed."
      call display_string_both(str_disp, a_book_keeping%fU)
    end if
  end do
  !
  if (a_disk_iter_params%flag_converged) then
    write(*, '(A/)') "Iteration has finished and converged!"
  else
    write(*, '(A/)') "Iteration has finished but not converged yet."
  end if
  !
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=iteration_configure)
  end if
  write(str_disp, '("!Final number of cells =", I10)') leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  !
end subroutine disk_iteration



subroutine montecarlo_prep
  integer i, i1, i2, j
  double precision lam_max
  type(type_stellar_params) star_tmp
  double precision, parameter :: T_Lya=1000D0, lam_max_0=1D8
  integer, parameter :: n_interval_star=2000, n_interval_Xray=200
  character(len=64) strTMP
  character(len=256) strTMP_L
  !
  write(*, '(A/)') 'Preparing for the Monte Carlo.'
  !
  mc_conf%starpos_r = a_disk%starpos_r
  mc_conf%starpos_z = a_disk%starpos_z
  !
  if (mc_conf%ph_init_symmetric) then
    mc_conf%maxw = get_surf_max_angle(a_disk%starpos_r, a_disk%starpos_z)
    mc_conf%minw = -mc_conf%maxw
  else
    mc_conf%maxw = get_surf_max_angle(a_disk%starpos_r, a_disk%starpos_z)
    mc_conf%minw = get_bott_min_angle(a_disk%starpos_r, a_disk%starpos_z)
  end if
  !
  write(*,'(A, 2ES12.4)') 'Star location r,z = ', &
        mc_conf%starpos_r, mc_conf%starpos_z
  write(*,'(A, 2ES12.4/)') 'minw,maxw = ', mc_conf%minw, mc_conf%maxw
  !
  call load_H2O_ab_crosssection( &
    combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_water), &
    water_0)
  !
  if (.not. a_disk%waterShieldWithRadTran) then
    water_0%ab = 0D0
    water_0%sc = 0D0
    water_0%g  = 0D0
  end if
  !
  call make_H_Lya(T_Lya, HI_0)
  !
  ! Let the dust and H2O, H data share the same wavelength axis
  call align_optical_data
  !
  ! Create the lookup table for finding temperature and new wavelength during
  ! the Monte Carlo
  call make_luts
  !
  call openFileSequentialWrite(i1, &
       combine_dir_filename(a_book_keeping%dir, 'optical_parameters.dat'), &
       999, getu=1)
  do i=1, HI_0%n
    strTMP_L = ''
    do j=1, dusts%n
      write(strTMP, '(3ES15.5E3)') &
        dusts%list(j)%ab(i), dusts%list(j)%sc(i), dusts%list(j)%g(i)
      strTMP_L = trim(strTMP_L) // strTMP
    end do
    write(i1, '(I8, 7ES15.5E3, "  ", A)') &
        i, HI_0%lam(i), &
        HI_0%ab(i),          HI_0%sc(i),          HI_0%g(i), &
        water_0%ab(i),       water_0%sc(i),       water_0%g(i), &
        trim(strTMP_L)
  end do
  close(i1)
  !
  ! Prepare for the stellar spectrum
  write(*, '(A)') 'Preparing for the stellar spectrum.'
  !
  lam_max = min(lam_max_0, dust_0%lam(dust_0%n)) ! in angstrom
  !
  ! Black body
  call make_stellar_spectrum(dust_0%lam(1), &
    lam_max, n_interval_star, a_star)
  call openFileSequentialWrite(i1, &
       combine_dir_filename(a_book_keeping%dir, 'stellar_spectrum_blackbody.dat'), &
       99, getu=1)
  do i=1, a_star%n
    write(i1, '(I8, 2ES20.10E3)') i, a_star%lam(i), a_star%vals(i)
  end do
  close(i1)
  !
  ! Thermal X-ray
  star_tmp%T_Xray = a_star%T_Xray
  star_tmp%lumi_Xray = a_star%lumi_Xray
  star_tmp%E0_Xray = a_star%E0_Xray
  star_tmp%E1_Xray = max(a_star%E1_Xray, &  ! The 1.01D0 here is for a stupid reason.
    phy_hPlanck_CGS*phy_SpeedOfLight_CGS/(a_star%lam(1)*1D-8)/phy_eV2erg/1D3*1.01D0)
  !
  call make_stellar_spectrum_Xray(n_interval_Xray, star_tmp)
  call openFileSequentialWrite(i1, &
       combine_dir_filename(a_book_keeping%dir, 'stellar_spectrum_Xray.dat'), &
       99, getu=1)
  do i=1, star_tmp%n
    write(i1, '(I8, 2ES20.10E3)') i, star_tmp%lam(i), star_tmp%vals(i)
  end do
  close(i1)
  !
  call merge_stellar_spectrum(star_tmp, a_star)
  !
  ! From observation
  if (.not. mc_conf%use_blackbody_star) then
    call load_stellar_spectrum( &
      trim(combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_star)), &
      star_tmp)
    !
    call openFileSequentialWrite(i1, &
         combine_dir_filename(a_book_keeping%dir, 'stellar_spectrum_fromfile.dat'), &
         99, getu=1)
    do i=1, star_tmp%n
      write(i1, '(I8, 2ES20.10E3)') i, star_tmp%lam(i), star_tmp%vals(i)
    end do
    close(i1)
    !
    call merge_stellar_spectrum(star_tmp, a_star)
  end if
  !
  ! Rescale the UV part of the stellar spectrum
  ! The edge index is extended by one to make the transition smoother
  do i=1, a_star%n
    if (a_star%lam(i+1) .gt. lam_range_UV(1)) then
      i1 = i
      exit
    end if
  end do
  do i=1, a_star%n
    if (a_star%lam(i) .gt. lam_range_UV(2)) then
      i2 = i
      exit
    end if
  end do
  do i=i1, i2
    a_star%vals(i) = a_star%vals(i) * mc_conf%stellar_spectr_UV_rescale_factor
  end do
  !
  call openFileSequentialWrite(i1, &
       combine_dir_filename(a_book_keeping%dir, 'stellar_spectrum_merged.dat'), &
       99, getu=1)
  do i=1, a_star%n
    write(i1, '(I8, 2ES20.10E3)') i, a_star%lam(i), a_star%vals(i)
  end do
  close(i1)
  !
  a_star%lumi = get_stellar_luminosity(a_star)
  a_star%lumi_UV = get_stellar_luminosity(a_star, lam_range_UV(1), &
    lam_range_UV(2))
  a_star%lumi_Lya = get_stellar_luminosity(a_star, lam_range_LyA(1), &
    lam_range_LyA(2))
  a_star%lumi_Vis = get_stellar_luminosity(a_star, lam_range_Vis(1), &
    lam_range_Vis(2))
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar total luminosity: ', a_star%lumi, ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar X-ray luminosity: ', &
    get_stellar_luminosity(a_star, lam_range_Xray(1), lam_range_Xray(2)), &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar X-ray luminosity (given): ', a_star%lumi_Xray, &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar UV luminosity: ', a_star%lumi_UV, ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') & ! Schindhelm 2012
    'Stellar UV luminosity (1160 to 1695): ', &
    get_stellar_luminosity(a_star, 1160D0, 1695D0), &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar Lyman alpha luminosity: ', a_star%lumi_Lya, &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar visual luminosity: ', a_star%lumi_Vis, &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar NIR luminosity: ', &
    get_stellar_luminosity(a_star, lam_range_NIR(1), lam_range_NIR(2)), &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar MIR luminosity: ', &
    get_stellar_luminosity(a_star, lam_range_MIR(1), lam_range_MIR(2)), &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  write(str_disp, '(A, ES16.6, A)') &
    'Stellar FIR luminosity: ', &
    get_stellar_luminosity(a_star, lam_range_FIR(1), lam_range_FIR(2)), &
    ' erg s-1.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  allocate(a_star%vals0(a_star%n))
  a_star%vals0 = a_star%vals
  a_star%lumi0 = a_star%lumi
  a_star%lumi_UV0 =  a_star%lumi_UV
  !
  ! The lumi and vals of a_star will be changed here.
  call get_mc_stellar_par(a_star, mc_conf)
  !
  ! Global optical property collection
  call make_global_coll
  !
  p4lam%n = luts%list(1)%m
  allocate(p4lam%pvals(0:p4lam%n))
  !
  if (mc_conf%collect_photon) then
    call set_up_collector(collector, &
      minlam=mc_conf%collect_lam_min, maxlam=mc_conf%collect_lam_max, &
      dmu=mc_conf%collect_dmu, nmu=mc_conf%collect_nmu, nr=mc_conf%collect_nr, &
      nphi=mc_conf%collect_nphi, &
      ang_mins=mc_conf%collect_ang_mins, ang_maxs=mc_conf%collect_ang_maxs)
  end if
  !
  write(*, '(A, I10)') 'Number of spectral channels:', dust_0%n
  write(*, '(A, I10)') 'Vectors per cell of this size (at least):', &
    opmaterials%ntype*2+7+3
  !
end subroutine montecarlo_prep



subroutine merge_stellar_spectrum(s1, s2)
  ! Merge spectra s1 and s2.  The result is stored in s2.
  ! In overlapping regions, s1 has a higher priority over s2, meaning that if
  ! both s1 and s2 has value at the same x value, the value of s2 will be
  ! overwritten by s1.
  type(type_stellar_params), intent(in) :: s1
  type(type_stellar_params), intent(inout) :: s2
  double precision, dimension(:), allocatable :: v1, v2
  integer n
  n = s1%n + s2%n
  allocate(v1(n), v2(n))
  call merge_vec(s1%n, s1%lam, s2%n, s2%lam, n, v1)
  call transfer_value(s2%n, s2%lam, s2%vals, n, v1, v2)
  call transfer_value(s1%n, s1%lam, s1%vals, n, v1, v2, keep=.true.)
  deallocate(s2%lam, s2%vals)
  allocate(s2%lam(n), s2%vals(n))
  s2%n = n
  s2%lam = v1
  s2%vals = v2
  deallocate(v1, v2)
end subroutine merge_stellar_spectrum



subroutine make_dusts_data
  integer i, j, itype, nradius, nradius_prev, nlam
  double precision rmin, rmax, ind, swei, m
  double precision, dimension(:), allocatable :: t1, t2
  !
  dusts%n = a_disk%ndustcompo
  allocate(dusts%list(dusts%n))
  nradius = 0
  !
  do i=1, dusts%n
    !
    call calc_dust_MRN_par(a_disk%dustcompo(i)%mrn)
    !
    rmin  = a_disk%dustcompo(i)%mrn%rmin
    rmax  = a_disk%dustcompo(i)%mrn%rmax
    ind   = a_disk%dustcompo(i)%mrn%n
    itype = a_disk%dustcompo(i)%itype
    !
    rmax = max(rmax, rmin*1.0001D0)
    !
    nradius_prev = nradius
    nradius = dustmix_data%list(itype)%nradius
    nlam = dustmix_data%list(itype)%nlam
    !
    if (i .gt. 1) then
      if (dusts%list(i-1)%n .ne. nlam) then
        write(*,'(A)') 'In make_dusts_data:'
        write(*,'(A)') 'Arrays for different dust types not '// &
                       'having the same dimension!'
        call error_stop()
      end if
      if (nradius .ne. nradius_prev) then
        write(*, '(A)') 'In make_dusts_data:'
        write(*, '(A)') 'Inconsistent radius array size!'
        call error_stop()
      end if
    end if
    !
    dusts%list(i)%n = nlam
    !
    allocate(dusts%list(i)%lam(nlam), &
             dusts%list(i)%ab(nlam), &
             dusts%list(i)%sc(nlam), &
             dusts%list(i)%g(nlam))
    !
    if (.not. allocated(t1)) then
      allocate(t1(nradius), t2(nradius))
    end if
    !
    do j=1, nlam
      !
      dusts%list(i)%lam(j) = &
        dustmix_data%list(itype)%w(j) / phy_Angstrom2micron
      !
      t1 = exp(-ind*log(dustmix_data%list(itype)%r)) ! = r**(-ind)
      !
      swei = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t1, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%ab(j, :)
      dusts%list(i)%ab(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%sc(j, :)
      dusts%list(i)%sc(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%g(j, :)
      dusts%list(i)%g(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      m = 4D0*phy_Pi/3D0 * a_disk%dustcompo(i)%mrn%r3av * &
          phy_micron2cm**3 * dustmix_info%mix(itype)%rho
      a_disk%dustcompo(i)%pmass_CGS = m ! dust particle mass in gram
      !
      ! Now the unit of ab and sc become cm2 g-1.
      dusts%list(i)%ab(j) = dusts%list(i)%ab(j) / swei * phy_micron2cm**2 / m
      dusts%list(i)%sc(j) = dusts%list(i)%sc(j) / swei * phy_micron2cm**2 / m
      dusts%list(i)%g(j) = dusts%list(i)%g(j) / swei
      !write(*, '(2I4, 10ES12.4)') i, j, dusts%list(i)%lam(j), &
      !  dusts%list(i)%ab(j), dusts%list(i)%sc(j), dusts%list(i)%g(j), &
      !  m, a_disk%dustcompo(i)%mrn%r3av, swei, rmin, rmax, maxval(dustmix_data%list(itype)%g(j, :))
    end do
  end do
  deallocate(t1, t2)
end subroutine make_dusts_data



subroutine do_grid_stuff(iiter, overwrite)
  integer, intent(in) :: iiter
  logical, intent(in), optional :: overwrite
  logical ovw
  if (present(overwrite)) then
    ovw = overwrite
  else
    ovw = .false.
  end if
  if ((iiter .gt. 1) .or. (.not. a_disk_iter_params%use_backup_grid_data)) then
    write(str_disp, '("! ", A)') "Dumping grid data..."
    call back_grid_info(dump_dir_out, iiter=iiter, dump=.true., overwrite=ovw)
    !
  end if
end subroutine do_grid_stuff


subroutine load_grid_grom_file
  allocate(root)
  call cell_init(root)
  !
  call back_grid_info(dump_dir_in, &
         fname=a_disk_iter_params%dump_filename_grid, dump=.false.)
  !
  call get_number_of_leaves(root)
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
end subroutine load_grid_grom_file



subroutine do_optical_stuff(iiter, overwrite)
  integer, intent(in) :: iiter
  logical, intent(in), optional :: overwrite
  logical ovw
  !
  if (present(overwrite)) then
    ovw = overwrite
  else
    ovw = .false.
  end if
  !
  if ((.not. a_disk_iter_params%use_backup_optical_data) .or. &
      (iiter .gt. 1)) then
    !
    if (a_disk_iter_params%rerun_whole .and. &
        a_disk_iter_params%rerun_whole_noMonteCarlo) then
      return
    end if
    !
    write(*, '(A)') 'Preparing optical data for all the cells.'
    call montecarlo_reset_cells
    !
    write(str_disp, '("! ", A)') "Monte Carlo begins."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // &
        trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call montecarlo_do(mc_conf, root)
    !
    write(str_disp, '("! ", A)') "Monte Carlo finished. Doing post-Monte Carlo."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! Retrieve physical parameters from the monte carlo results
    call post_montecarlo
    !
    write(str_disp, '("! ", A, I16)') 'Number of photon packets used:', mc_conf%icount
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // &
        trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    write(str_disp, '(A, I6)') '! iiter = ', iiter
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    if (mc_conf%collect_photon) then
      ! Save the spectrum generated from the collected photons that have escaped.
      call save_collected_photons_iter(iiter)
    end if
    !
    write(str_disp, '("! ", A)') "Dumping optical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call back_cells_optical_data(dump_dir_out, iiter=iiter, &
        dump=.true., overwrite=ovw)
    !
    write(str_disp, '("! ", A)') "Dumping preliminary physical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call back_cells_physical_data(dump_dir_out, iiter=iiter, &
        dump=.true., preliminary=.true., overwrite=ovw)
    !
    write(str_disp, '("! ", 2A)') 'Data dumped in ', trim(dump_dir_out)
    call display_string_both(str_disp, a_book_keeping%fU)
    !
  else
    write(str_disp, '("! ", A)') "Loading backuped optical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    write(*, '(A)') 'Preparing optical data for all the cells.'
    call montecarlo_reset_cells
    !
    call back_cells_optical_data(dump_dir_in, &
           fname=a_disk_iter_params%dump_filename_optical, dump=.false.)
    !
    if ((len_trim(a_disk_iter_params%dump_filename_physical) .gt. 0) .and. &
        a_disk_iter_params%use_backup_physical_data .and. &
        (.not. a_disk_iter_params%use_backup_chemical_data)) then
      ! When intending to use previous optical data for further chemical
      ! calculation, some associated physical data must also be loaded for the
      ! chemical part to run properly.
      call back_cells_physical_data(dump_dir_in, &
           fname=a_disk_iter_params%dump_filename_physical, dump=.false.)
    end if
    !
  end if
end subroutine do_optical_stuff



subroutine do_chemical_stuff(iiter)
  integer, intent(in) :: iiter
  integer i, i0, i_count, l_count
  !
  if ((.not. a_disk_iter_params%use_backup_chemical_data) .or. &
      (iiter .gt. 1)) then
    !
    if ((iiter .eq. 1) .and. (a_disk_iter_params%use_backup_physical_data)) then
      write(str_disp, '("! ", A)') &
          "Loading backuped physical data..."
      call display_string_both(str_disp, a_book_keeping%fU)
      call back_cells_physical_data(dump_dir_in, &
             fname=a_disk_iter_params%dump_filename_physical, dump=.false.)
      !
      if (len_trim(a_disk_iter_params%dump_filename_physical_aux) .gt. 0) then
        call back_cells_physical_data_aux(dump_dir_in, &
             fname=a_disk_iter_params%dump_filename_physical_aux, dump=.false.)
      end if
      !
      call check_consistency_of_loaded_data_phy
    end if
    ! Start from the inner edge layer.
    n_calculating_cells = columns_idx(1)%nlen
    calculating_cells(1:n_calculating_cells) = columns_idx(1)%vals
    !
    i_count = 0 ! Counter for cells
    l_count = 0 ! Counter for layers
    !
    do
      l_count = l_count + 1
      !
      do i=1, n_calculating_cells
        i_count = i_count + 1
        i0 = calculating_cells(i)
        !
        write(*, '(4(A, I5, A, I5, ",", 2X), A, 4ES12.3)') &
          "Iter:", a_disk_iter_params%n_iter_used, "/", a_disk_iter_params%n_iter, &
          "Cell:", i_count, '/', leaves%nlen, &
          "Row:", i, '/', n_calculating_cells, &
          "Col:", l_count, '/', n_columns, &
          'rz:', &
          leaves%list(i0)%p%xmin, &
          leaves%list(i0)%p%xmax, &
          leaves%list(i0)%p%ymin, &
          leaves%list(i0)%p%ymax
        write(*, '(2(A, ES10.3, 2X), 2X, 2A, 2X, 2A)') &
          'n_gas: ', leaves%list(i0)%p%par%n_gas, &
          'Tdust: ', leaves%list(i0)%p%par%Tdust, &
          'exe: ', trim(a_disk_iter_params%filename_exe), &
          'dir: ', trim(a_disk_iter_params%iter_files_dir)
        !
        call calc_this_cell(i0)
        !
        call check_convergency_cell(i0)
        !
        write(*, '(12X, 10A10)') chem_idx_some_spe%names(1:10)
        write(*, '(A, 2X, 10ES10.3, L3/)') 'Abundances:',  &
          leaves%list(i0)%p%abundances(chem_idx_some_spe%idx(1:10)), &
          leaves%list(i0)%p%converged
        !
        call disk_save_results_write(fU_save_results, leaves%list(i0)%p)
        flush(fU_save_results)
      end do
      !
      !call update_calculating_cells
      if (l_count .ge. n_columns) then
        exit
      end if
      n_calculating_cells = columns_idx(l_count+1)%nlen
      calculating_cells(1:n_calculating_cells) = columns_idx(l_count+1)%vals
      !
      if (n_calculating_cells .le. 0) then
        exit
      end if
    end do
    !
    write(str_disp, '(A)') '! Current time: ' // &
        trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    write(str_disp, '(A, I6)') '! iiter = ', iiter
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    write(str_disp, '("! ", A)') "Dumping chemical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call back_cells_chemical_data(dump_dir_out, iiter=iiter, dump=.true.)
    !
    write(str_disp, '("! ", 2A)') 'Data dumped in ', trim(dump_dir_out)
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    write(str_disp, '("! ", A)') "Dumping physical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call back_cells_physical_data(dump_dir_out, iiter=iiter, dump=.true.)
    !
    call back_cells_physical_data_aux(dump_dir_out, iiter=iiter, dump=.true.)
    !
    write(str_disp, '("! ", 2A)') 'Data dumped in ', trim(dump_dir_out)
    call display_string_both(str_disp, a_book_keeping%fU)
  else
    write(str_disp, '("! ", A)') &
        "Loading backuped chemical and physical data..."
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call back_cells_chemical_data(dump_dir_in, &
           fname=a_disk_iter_params%dump_filename_chemical, dump=.false.)
    !
    call back_cells_physical_data(dump_dir_in, &
           fname=a_disk_iter_params%dump_filename_physical, dump=.false.)
    !
    call back_cells_physical_data_aux(dump_dir_in, &
           fname=a_disk_iter_params%dump_filename_physical_aux, dump=.false.)
    !
    call check_consistency_of_loaded_data_phy
  end if
end subroutine do_chemical_stuff



subroutine do_vertical_struct_with_Tdust
  logical vertIterCvg
  integer iVertIter
  double precision fr_min, fr_max, nd_min, ng_min
  !
  if (.not. a_disk_iter_params%do_vertical_with_Tdust) then
    return
  end if
  !
  fr_max = 1D99
  fr_min = 1D-99
  !
  write(str_disp, '(A)') &
    'Doing vertical structure calculation based on dust temperature.'
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '(A)') '! Current time: ' // &
    trim(a_date_time%date_time_str())
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  vertIterCvg = .false.
  !
  ! Calculate the column gravty force.
  call post_vertical_structure_adj
  !
  call do_save_only_structure
  !
  do iVertIter=1, a_disk_iter_params%nVertIterTdust
    write(str_disp, '(A, I4)') 'Vertical structure with Tdust.  Iter ', iVertIter
    call display_string_both(str_disp, a_book_keeping%fU, beginwithnewline=.true.)
    !
    if (a_disk_iter_params%do_simple_chem_before_mc) then
      call do_simple_chemistry(iVertIter)
    end if
    !
    call do_optical_stuff(1, overwrite=.true.)
    !
    if (a_disk_iter_params%vertical_structure_fix_grid) then
      if ((iVertIter .le. max(a_disk_iter_params%nVertIterTdust-4, a_disk_iter_params%nVertIterTdust*2/3)) .and. &
          ((fr_max .gt. 10D0) .or. (fr_min .lt. 0.1D0))) then
        nd_min = grid_config%min_val_considered*1D-40
        ng_min = grid_config%min_val_considered*1D-20
      else if ((fr_max .le. 10D0) .and. (fr_min .ge. 0.1D0)) then
        nd_min = grid_config%min_val_considered_use*1D-17
        ng_min = grid_config%min_val_considered_use
      else
        nd_min = grid_config%min_val_considered_use*1D-30
        ng_min = grid_config%min_val_considered_use*1D-10
      end if
      call vertical_pressure_gravity_balance_alt(a_disk%star_mass_in_Msun, &
        useTdust=.true., Tdust_lowerlimit=a_disk_iter_params%minimum_Tdust, &
        ngas_lowerlimit=ng_min, &
        ndust_lowerlimit=nd_min, &
        fix_dust_struct=a_disk_iter_params%vertical_structure_fix_dust, &
        disk_gas_mass_preset = a_disk%disk_mass_in_Msun, &
        maxfac=fr_max, minfac=fr_min)
      !
      ! The number of using cells may have changed.
      call remake_index
      !
      if ((iVertIter .ge. 1) .and. &
          ((fr_max .gt. 2D0) .or. (fr_min .lt. 5D-1))) then
        write(*, '(A)') 'Merging cells.'
        call merge_cells
        write(*, '(A)') 'Refining the vertical structure.'
        call refine_after_vertical
      end if
      !
    else
      call vertical_pressure_gravity_balance(frescale_max=fr_max, &
        frescale_min=fr_min, useTdust=.true.)
      !
      call remake_index
      !
      write(*, '(A, 4ES16.6)') 'New bounding box:', &
          root%xmin, root%xmax, root%ymin, root%ymax
      write(*, '(A, 4ES16.6)') 'Inner top cell:', &
        leaves%list(surf_cells%idx(1))%p%xmin, &
        leaves%list(surf_cells%idx(1))%p%xmax, &
        leaves%list(surf_cells%idx(1))%p%ymin, &
        leaves%list(surf_cells%idx(1))%p%ymax
      write(*, '(A, 4ES16.6)') 'Outer top cell:', &
        leaves%list(surf_cells%idx(surf_cells%nlen))%p%xmin, &
        leaves%list(surf_cells%idx(surf_cells%nlen))%p%xmax, &
        leaves%list(surf_cells%idx(surf_cells%nlen))%p%ymin, &
        leaves%list(surf_cells%idx(surf_cells%nlen))%p%ymax
    end if
    !
    write(*, '(A, I8)') 'Number of bottom cells:', bott_cells%nlen
    write(*, '(A, I8)') 'Number of surface cells:', surf_cells%nlen
    write(str_disp, '(A, 2ES16.6)') 'rescale_max, rescale_min: ', fr_max, fr_min
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call post_vertical_structure_adj
    !
    write(*, '(A)') 'Saving structure information in ascii.'
    call do_save_only_structure
    !
    if ((fr_max .le. 2D0) .and. (fr_min .ge. 5D-1)) then
      vertIterCvg = .true.
      exit
    end if
  end do
  !
  if (vertIterCvg) then
    write(str_disp, '(A)') 'Vertical structure converged (with Tdust).'
  else
    write(str_disp, '(A)') 'Vertical structure has not converged (with Tdust).'
  end if
  call display_string_both(str_disp, a_book_keeping%fU)
end subroutine do_vertical_struct_with_Tdust



subroutine allocate_iter_stor
  integer i
  !
  if (allocated(a_iter_stor%T_s)) then
    deallocate(a_iter_stor%T_s, a_iter_stor%abundances)
  end if
  allocate(a_iter_stor%T_s(leaves%nlen), &
           a_iter_stor%abundances(chem_idx_some_spe%nItem, &
                                          leaves%nlen))
  !
  do i=1, leaves%nlen
    !leaves%list(i)%p%abundances = chemsol_stor%y(1:chem_species%nSpecies)
    a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
    a_iter_stor%abundances(:, i) = chemsol_stor%y(chem_idx_some_spe%idx)
  end do
end subroutine allocate_iter_stor



subroutine post_vertical_structure_adj
  integer i
  mc_conf%maxw = min(get_surf_max_angle(), 0.99D0)
  if (mc_conf%ph_init_symmetric) then
    mc_conf%minw = -mc_conf%maxw
  end if
  !
  write(str_disp, '(A, 2F9.4)') 'Using minw,maxw: ', mc_conf%minw, mc_conf%maxw
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  call get_mc_stellar_par(a_star, mc_conf)
  !
  do i=1, leaves%nlen
    call calc_gravity_single_cell(leaves%list(i)%p)
  end do
  do i=1, leaves%nlen
    call calc_gravity_column(leaves%list(i)%p)
  end do
end subroutine post_vertical_structure_adj



subroutine do_simple_chemistry(iiter)
  ! Calculate the atomic hydrogen abundance analytically
  integer, intent(in) :: iiter
  integer i, j, id, iH2PHD
  double precision kf, kd
  !
  iH2PHD = 0
  do i=1, chem_net%nReactions
    if ((chem_net%ctype(i) .eq. 'PH') .and. &
        (chem_net%reac_names(1, i) .eq. 'H2')) then
      iH2PHD = i
      exit
    end if
  end do
  if (iH2PHD .eq. 0) then
    write(*, '(A)') 'Cannot find H2 photodissociation reaction!'
    return
  end if
  !
  do i=1, n_columns
    do j=1, columns_idx(i)%nlen
      id = columns_idx(i)%vals(j)
      !
      call update_params_above_alt(id)
      !
      associate(c => leaves%list(id)%p)
        if (iiter .eq. 1) then
          c%par%Tgas = 100D0
          c%par%flux_UV_star_unatten = &
            (a_star%lumi_UV0 - a_star%lumi_Lya) / &
            (4D0*phy_Pi * &
              ((((c%xmax+c%xmin)*0.5D0)**2 + &
                ((c%ymax+c%ymin)*0.5D0)**2) * phy_AU2cm**2))
          c%par%G0_UV_toStar = c%par%flux_UV_star_unatten / &
            phy_Habing_energy_flux_CGS
          c%par%G0_UV_toISM = c%par%UV_G0_factor_background
          call calc_Ncol_to_ISM(c)
          call calc_Ncol_to_Star(c)
          c%par%Av_toISM  = phy_colDen2Av_coeff * c%par%Ncol_toISM
          c%par%Av_toStar = phy_colDen2Av_coeff * c%par%Ncol_toStar
        else
          c%par%Tgas = c%par%Tdust
        end if
        !
        call calc_local_dynamics(leaves%list(id)%p)
        !
        kf = 0.5D0 * c%par%sigdust_ave * c%par%ndust_tot * &
          sqrt(phy_kBoltzmann_CGS * c%par%Tgas / phy_mProton_CGS)
        kd = c%par%zeta_cosmicray_H2 * &
             exp(-chem_params%Ncol_toISM / const_cosmicray_attenuate_N) &
             + &
             chem_net%ABC(1, iH2PHD) * ( &
               c%par%G0_UV_toISM &
                 * exp(-chem_net%ABC(3, iH2PHD) * c%par%Av_toISM) &
                 * c%par%f_selfshielding_toISM_H2 + &
               c%par%G0_UV_toStar &
                 * exp(-chem_net%ABC(3, iH2PHD) * c%par%Av_toStar) &
                 * c%par%f_selfshielding_toStar_H2)
        c%abundances(chem_idx_some_spe%i_H2) = kf / (2D0*kf + kd)
        c%abundances(chem_idx_some_spe%i_HI) = kd / (2D0*kf + kd)
      end associate
    end do
  end do
end subroutine do_simple_chemistry


subroutine post_montecarlo
  integer i, j
  integer i1, i2
  double precision vx, vy, vz
  double precision RR, tmp, tmp0, tmp1, tmp2, tmp3, tmp4
  !
  tmp2 = 1D99
  tmp4 = 0D0
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      tmp = 0D0
      tmp0 = 0D0
      do j=1, dusts%n
        if (c%par%mdusts_cell(j) .le. 1D-50) then
          c%par%Tdusts(j) = 0D0
          cycle
        end if
        !
        tmp3 = (c%par%en_gains(j) + c%par%en_exchange(j)) &
               / (4*phy_Pi*c%par%mdusts_cell(j))
        if (isnan(tmp3)) then
          write(*, '(A)') 'NaN in post_montecarlo:'
          write(*, '(A, I4, 7ES12.4)') &
            'j, egains, eexs, mds, nds, rhods, abs, vol', j, &
            c%par%en_gains(j), c%par%en_exchange(j), c%par%mdusts_cell(j), &
            c%par%n_dusts(j), c%par%rho_dusts(j), c%par%abso_wei(j), c%par%volume
          write(*, '(A, 2ES12.4)') 'xmin,ymin: ', c%xmin, c%ymin
          call error_stop()
        end if
        c%par%Tdusts(j) = get_Tdust_from_LUT(tmp3, luts%list(j), i1)
        !
        if (c%par%Tdusts(j) .le. 0D0) then
          write(*, '(A, 4ES16.6, I4, 2ES16.6)') 'Tdusts(j)=0: ', &
            c%xmin, c%xmax, c%ymin, c%ymax, j, c%par%en_gains(j), c%par%en_exchange(j)
        end if
        !
        tmp1 = c%par%n_dusts(j) * a_disk%dustcompo(j)%mrn%r2av
        tmp0 = tmp0 + c%par%Tdusts(j) * tmp1
        tmp = tmp + tmp1
        !
      end do
      !else
      !end if
      if (tmp .gt. 0D0) then
        c%par%Tdust = max(tmp0 / tmp, a_disk_iter_params%minimum_Tdust)
      else
        c%par%Tdust = a_disk_iter_params%minimum_Tdust
      end if
      !
      c%par%en_gain_tot = sum(c%par%en_gains)
      c%par%en_gain_abso_tot = sum(c%par%en_gains_abso)
      !
      ! Flux of each cell as a function of wavelength
      c%optical%flux = c%optical%flux * (phy_AU2cm / c%par%volume)
      !
      if (mc_conf%do_fill_blank) then
        call fill_blank(dust_0%lam, c%optical%flux, c%optical%phc, &
                        c%optical%nlam, mc_conf%fill_blank_threshold, 3+c%optical%nlam/100)
      end if
      !
      ! Get some properties of the radiation field
      ! Only use the wavelength vector of dust_0
      !
      ! Total
      i1 = 1
      i2 = dust_0%n
      c%par%flux_tot = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_tot)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_tot)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_tot)
      c%par%dir_tot_r = vx
      c%par%dir_tot_z = vz
      c%par%aniso_tot = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! X-ray
      i1 = max(1, get_idx_for_kappa(lam_range_Xray(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_Xray(2), dust_0)) - 1
      c%par%flux_Xray = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Xray)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Xray)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Xray)
      c%par%dir_Xray_r = vx
      c%par%dir_Xray_z = vz
      c%par%aniso_Xray = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! UV
      i1 = max(1, get_idx_for_kappa(lam_range_UV(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_UV(2), dust_0)) - 1
      c%par%flux_UV = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      c%par%dir_UV_r = vx
      c%par%dir_UV_z = vz
      c%par%aniso_UV = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! Lya
      i1 = max(1, get_idx_for_kappa(lam_range_LyA(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_LyA(2), dust_0)) - 1
      c%par%flux_Lya = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      c%par%dir_Lya_r = vx
      c%par%dir_Lya_z = vz
      c%par%aniso_Lya = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! Visual
      i1 = max(1, get_idx_for_kappa(lam_range_Vis(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_Vis(2), dust_0)) - 1
      c%par%flux_Vis = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Vis)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Vis)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Vis)
      c%par%dir_Vis_r = vx
      c%par%dir_Vis_z = vz
      c%par%aniso_Vis = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! NIR
      i1 = max(1, get_idx_for_kappa(lam_range_NIR(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_NIR(2), dust_0)) - 1
      c%par%flux_NIR = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      c%par%dir_NIR_r = vx
      c%par%dir_NIR_z = vz
      c%par%aniso_NIR = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! MIR
      i1 = max(1, get_idx_for_kappa(lam_range_MIR(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_MIR(2), dust_0)) - 1
      c%par%flux_MIR = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      c%par%dir_MIR_r = vx
      c%par%dir_MIR_z = vz
      c%par%aniso_MIR = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! FIR
      i1 = max(1, get_idx_for_kappa(lam_range_FIR(1), dust_0)) + 1
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_FIR(2), dust_0)) - 1
      c%par%flux_FIR = sum(c%optical%flux(i1:i2))
#ifdef SAVE_PHOTON_FIELD_DIR
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      c%par%dir_FIR_r = vx
      c%par%dir_FIR_z = vz
      c%par%aniso_FIR = sqrt(vx**2 + vy**2 + vz**2)
#endif
      !
      ! Local number flux of Lyman alpha
      c%par%phflux_Lya = c%par%flux_Lya / phy_LyAlpha_energy_CGS
      c%par%G0_Lya_atten = c%par%flux_Lya / phy_Habing_energy_flux_CGS
      !
      ! Calculate the total column density to the star and to the ISM
      call calc_Ncol_to_ISM(leaves%list(i)%p)
      call calc_Ncol_to_Star(leaves%list(i)%p)
      !
      RR = (((c%xmax+c%xmin)*0.5D0)**2 + &
            ((c%ymax+c%ymin)*0.5D0)**2) * phy_AU2cm**2
      c%par%flux_UV_star_unatten = a_star%lumi_UV0 / (4D0*phy_Pi*RR)
      c%par%flux_Lya_star_unatten = a_star%lumi_Lya / (4D0*phy_Pi*RR)
      !
      ! 2014-06-18 Wed 23:47:06
      ! Now only UV continuum is included.  The Lya line is subtracted.
      c%par%flux_UV_star_unatten = c%par%flux_UV_star_unatten - &
                                   c%par%flux_Lya_star_unatten
      c%par%flux_UV = c%par%flux_UV - c%par%flux_Lya
      if (c%par%flux_UV_star_unatten .le. 0D0) then
        write(*,*) c%par%flux_Lya_star_unatten, c%par%flux_UV_star_unatten
        write(*,*) 'No UV!'
        stop
      end if
      !
      ! Calculate the G0 factors
      ! The G0 is the unattenuated one, so a further
      ! exp(-k*Av) should be applied.
      c%par%G0_UV_toStar = c%par%flux_UV_star_unatten / phy_Habing_energy_flux_CGS
      c%par%G0_UV_toISM  = c%par%UV_G0_factor_background
      !
      if (a_disk_iter_params%calc_Av_toStar_from_Ncol) then
        c%par%Av_toStar = 1.086D0 * &
          calc_Ncol_from_cell_to_point(c, 0D0, 0D0, -6, fromCellCenter=.false.)
        c%par%G0_UV_toStar_photoDesorb = &
            c%par%G0_UV_toStar * exp(-c%par%Av_toStar/1.086D0*phy_UVext2Av)
        !
        ! UV to dissociate H2
        c%par%G0_UV_H2phd = &
          get_stellar_luminosity(a_star, lam_range_UV_H2phd(1), & 
            lam_range_UV_H2phd(2)) &
          / (4D0*phy_Pi*RR) / phy_Habing_energy_flux_CGS &
          * exp(-c%par%Av_toStar/1.086D0*phy_UVext2Av)
      else
        if ((c%par%flux_UV .le. 0D0) .or. (c%par%flux_UV_star_unatten .le. 0D0)) then
          c%par%Av_toStar = 1D99
        else
          c%par%Av_toStar = min(1D99, max(0D0, &
            -1.086D0 * log(c%par%flux_UV / c%par%flux_UV_star_unatten) / phy_UVext2Av))
        end if
        c%par%G0_UV_toStar_photoDesorb = c%par%flux_UV / phy_Habing_energy_flux_CGS
        !
        ! UV to dissociate H2
        i1 = max(1, get_idx_for_kappa(lam_range_UV_H2phd(1), dust_0)) + 1
        i2 = min(dust_0%n, get_idx_for_kappa(lam_range_UV_H2phd(2), dust_0)) - 1
        c%par%G0_UV_H2phd = sum(c%optical%flux(i1:i2)) / phy_Habing_energy_flux_CGS
      end if
      !
      ! The Av to ISM is a simple scaling of the dust column density
      ! The factor 2 is to account for the scattering.
      c%par%Av_toISM = 1.086D0 * &
        calc_Ncol_from_cell_to_point(c, (c%xmin+c%xmax)*0.5D0, &
          root%ymax*2D0, -6, fromCellCenter=.false.)
      tmp2 = min(tmp2, c%par%flux_UV)
      tmp4 = max(tmp4, c%par%flux_UV)
    end associate
  end do
  write(str_disp, '(A, ES12.2)') 'Min G0_UV: ', tmp2 / phy_Habing_energy_flux_CGS
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '(A, ES12.2)') 'Max G0_UV: ', tmp4 / phy_Habing_energy_flux_CGS
  call display_string_both(str_disp, a_book_keeping%fU)
end subroutine post_montecarlo


pure subroutine fill_blank(x, v, mask, n, nth, nrange)
  integer, intent(in) :: n, nth, nrange
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(inout) :: v
  integer, dimension(n), intent(in) :: mask
  integer i, j, jmin, jmax
  double precision s, smean
  do i=1, n
    if (mask(i) .lt. nth) then
      jmin = n
      jmax = 1
      do j=i-1, 1, -1
        if (mask(j) .ge. nth) then
          jmin = j
          exit
        end if
      end do
      do j=i+1, n
        if (mask(j) .ge. nth) then
          jmax = j
          exit
        end if
      end do
      jmin = min(jmin, max(1, i-nrange))
      jmax = max(jmax, min(n, i+nrange))
      s = 0D0
      do j=jmin, jmax-1
        s = s + v(j)
      end do
      smean = s / abs(x(jmax) - x(jmin))
      do j=jmin, jmax-1
        v(j) = smean * abs(x(j+1) - x(j))
      end do
    end if
  end do
end subroutine fill_blank



subroutine montecarlo_reset_cells
  integer i
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      !
      if ((.not. associated(c%par)) .or. &
          (.not. allocated(c%abundances))) then
        write(*, '(A)') 'Par or abundances not allocated!'
        write(*, '(A)') 'In montecarlo_reset_cells.'
        write(*, '(A, I8)') 'i = ', i
        call error_stop()
      end if
      !
      c%par%Tdusts = 0D0
      !
      c%par%X_HI  = c%abundances(chem_idx_some_spe%i_HI)
      c%par%X_H2O = c%abundances(chem_idx_some_spe%i_H2O)
      !
      !write(*, '(A, I8)') 'Allocating optical arrays.', i
      call allocate_local_optics(leaves%list(i)%p, &
                                 opmaterials%ntype, dust_0%n)
      !write(*, '(A, I8)') 'Resetting optical data.', i
      call reset_local_optics(leaves%list(i)%p)
    end associate
  end do
  !
  if (allocated(collector%energy)) then
    collector%energy = 0D0
    collector%counts = 0
  end if
  !
end subroutine montecarlo_reset_cells



subroutine disk_iteration_prepare
  integer i
  !
  ! The density structure is needed for making the grid.
  a_andrews_4ini = a_disk%andrews_gas
  a_andrews_4ini%particlemass = 1.4D0 * phy_mProton_CGS
  !
  if (abs(a_andrews_4ini%gam - 2D0) .le. 1D-4) then
    write(*, '(A)') 'a_andrews_4ini%gam cannot be 2.0!'
    call error_stop()
  end if
  do i=1, a_disk%ndustcompo
    if (abs(a_disk%dustcompo(i)%andrews%gam-2D0) .le. 1D-4) then
      write(*, '(A, I4)') 'a_disk%dustcompo(i)%andrews%gam cannot be 2.0!', i
      call error_stop()
    end if
  end do
  !
  a_star%mass      = a_disk%star_mass_in_Msun
  a_star%radius    = a_disk%star_radius_in_Rsun
  a_star%T         = a_disk%star_temperature
  !
  a_star%T_Xray    = a_disk%T_Xray
  a_star%E0_Xray   = a_disk%E0_Xray
  a_star%E1_Xray   = a_disk%E1_Xray
  a_star%lumi_Xray = a_disk%lumi_Xray
  !
  if (a_disk_iter_params%use_backup_grid_data) then
    write(*, '(A/)') 'Loading grid from file.'
    call load_grid_grom_file
  else
    write(*, '(A/)') 'Making grid.'
    call make_grid
  end if
  !
  n_calculating_cells_max = leaves%nlen
  allocate(calculating_cells(n_calculating_cells_max))
  !
  write(*, '(A/)') 'Preparing dust data.'
  call prep_dust_data ! Load dust data, and create the mixtures
  !
  call make_dusts_data ! Prepare the dust optical data for use
  !
  ! Prepare the chemical stuff
  chemsol_params%fU_log = a_book_keeping%fU
  !
  write(*, '(A/)') 'Loading chemcial reactions.'
  call chem_read_reactions()
  call chem_load_reactions()
  call chem_parse_reactions()
  call chem_get_dupli_reactions()
  call chem_get_idx_for_special_species()
  call chem_load_species_enthalpies
  call chem_get_reaction_heat
  call save_chem_rates(0) ! Save back the imported network.
  !
  call load_refine_check_species
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare_run_once
  !
  call chem_load_initial_abundances
  !
  write(str_disp, '("!", A, 2X, I8)') 'Number of cells (leaf):', leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I8)') 'Number of cells (total):', root%nOffspring
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I8)') 'Number of reactions:', chem_net%nReactions
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I8)') 'Number of species ', chem_species%nSpecies
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  !! Set the disk and cell parameters
  !write(*, '(A/)') 'Setting disk parameters.'
  !call disk_set_disk_params
  !
  write(*, '(A/)') 'Setting cell parameters.'
  call disk_set_gridcell_params_runonce
  !
  call make_columns
  !
  call allocate_iter_stor
  !
  call disk_calc_disk_mass
  !
  write(*, '(A/)') 'Preparing for the heating cooling parameters.'
  !
  call heating_cooling_prepare
  !
  !if (a_disk%allow_gas_dust_en_exch) then
  !  heating_cooling_config%use_mygasgraincooling = .true.
  !end if
  !
  if (a_disk_ana_params%do_analyse) then
    call load_ana_snippet
  end if
  !
  mc_conf%mc_dir_out = trim( &
    combine_dir_filename( &
      a_disk_iter_params%iter_files_dir, &
      mc_conf%mc_dir_out))
  if (.not. dir_exist(mc_conf%mc_dir_out)) then
    call my_mkdir(mc_conf%mc_dir_out)
  end if
  !
end subroutine disk_iteration_prepare



subroutine calc_this_cell(id)
  integer, intent(in) :: id
  integer i, isav, j
  double precision tmp, Tgas0
  double precision, dimension(const_nElement) :: ele_bef, ele_aft
  !
  leaves%list(id)%p%iIter = a_disk_iter_params%n_iter_used
  !
  if (chemsol_params%flag_chem_evol_save) then
    chemsol_params%chem_evol_save_filename = &
      trim(combine_dir_filename( &
        a_disk_iter_params%iter_files_dir, 'chem_evol_tmp.dat'))
  end if
  !
  hc_params => chem_params
  !
  ! Set the initial condition for chemical evolution
  call set_initial_condition_4solver(id, 1, leaves%list(id)%p%iIter)
  !
  Tgas0 = 0D0
  leaves%list(id)%p%par%t_final = 0D0
  !
  do j=1, a_disk_iter_params%nlocal_iter
    !
    chem_params = leaves%list(id)%p%par
    !
    write(*, '("Local iter: ", I4, " of ", I4)') j, &
        a_disk_iter_params%nlocal_iter
    !
    if (j .gt. 1) then
      !
      call set_initial_condition_4solver_continue(id, j)
      !
      write(*, '(A, ES16.6)') 'Continue running from tout=', chemsol_params%t0
      if ((abs(Tgas0 - chem_params%Tgas) .le. 1D-2 * Tgas0) .and. &
          (leaves%list(id)%p%par%t_final .ge. 0.2D0*chemsol_params%t_max) .and. &
          (heating_cooling_is_very_slow())) then
        chemsol_params%evolT = .false.
      end if
    end if
    Tgas0 = chem_params%Tgas
    !
    call chem_set_solver_flags_alt(j)
    !
    write(*, '(4X, A, F12.3/)') 'Tgas_old: ', leaves%list(id)%p%par%Tgas
    !
    call update_params_above_alt(id)
    call calc_local_dynamics(leaves%list(id)%p)
    !
    call set_hc_chem_params_from_cell(id)
    !
    ! For checking the elemental conservation.
    call get_elemental_abundance(chemsol_stor%y, chemsol_params%NEQ, &
        ele_bef, const_nElement)
    !
    current_cell_ptr => leaves%list(id)%p
    !
    call chem_evol_solve
    !
    call get_elemental_abundance(chemsol_stor%y, chemsol_params%NEQ, &
        ele_aft, const_nElement)
    !
    write(*, '(A)') 'Elemental abundances (before, after, diff):'
    do i=1, const_nElement
      if ((abs(ele_aft(i)) .le. 1D-90) .and. (abs(ele_bef(i)) .le. 1D-90)) then
        tmp = 0D0
      else
        tmp = (ele_aft(i) - ele_bef(i)) / (ele_bef(i) + ele_aft(i))
      end if
      if (abs(tmp) .ge. 1D-6) then
        write(*, '(4X, A8, 3ES16.6)') const_nameElements(i), &
          ele_bef(i), ele_aft(i), tmp
      end if
    end do
    !
    if ((j .gt. 1) .and. &
        (chemsol_stor%touts(chemsol_params%n_record_real) .le. &
         leaves%list(id)%p%par%t_final)) then
      write(str_disp, '(A)') 'Local iteration does not proceed!!!'
      call display_string_both(str_disp, a_book_keeping%fU)
      !
      write(str_disp, '(A, ES16.6)') 'Stop local iteration.  Tgas=', &
          chemsol_stor%y(chem_species%nSpecies+1)
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    end if
    !
    do isav=chemsol_params%n_record_real, 1, -1
      if ((.not. isnan(chemsol_stor%record(chem_species%nSpecies+1, isav))) .and. &
          (.not. isnan(chemsol_stor%record(chem_idx_some_spe%i_H2, isav)))) then
        exit
      end if
    end do
    !
    leaves%list(id)%p%quality = chemsol_params%quality
    !
    if (isav .le. 1) then
      write(str_disp, '(A)') 'No useful data produced!'
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    end if
    !
    leaves%list(id)%p%abundances = chemsol_stor%record(1:chem_species%nSpecies, isav)
    leaves%list(id)%p%par%Tgas = chemsol_stor%record(chem_species%nSpecies+1, isav)
    leaves%list(id)%p%par%t_final = chemsol_stor%touts(isav)
    !
    ! Update chem_params%n_mol_on_grain, as a side effect of the function call; this is bad.
    tmp = get_ice_coverage(chem_species%nSpecies, leaves%list(id)%p%abundances)
    !
    leaves%list(id)%p%par%Tdusts = chem_params%Tdusts
    leaves%list(id)%p%par%Tdust  = chem_params%Tdust
    leaves%list(id)%p%par%R_H2_form_rate  = chem_params%R_H2_form_rate
    leaves%list(id)%p%par%n_mol_on_grain  = chem_params%n_mol_on_grain
    !
    if (isnan(leaves%list(id)%p%par%Tgas) .or. &
        (leaves%list(id)%p%par%Tgas .le. 0D0) .or. &
        (leaves%list(id)%p%par%Tgas .gt. a_disk_iter_params%Tgas_crazy)) then
      write(str_disp, '(A, ES16.6)') 'Tgas is crazy! ', leaves%list(id)%p%par%Tgas
      call display_string_both(str_disp, a_book_keeping%fU)
    end if
    !
    write(*, '(2X, 2(2X, A, F12.3))') 'Tgas_new: ', leaves%list(id)%p%par%Tgas, &
      'Tdust: ', leaves%list(id)%p%par%Tdust
    !
    call update_params_above_alt(id)
    call calc_local_dynamics(leaves%list(id)%p)
    !
    call set_hc_chem_params_from_cell(id)
    !
    tmp = heating_minus_cooling()
    call update_en_exchange_with_dust
    !
    leaves%list(id)%p%h_c_rates = heating_cooling_rates
    leaves%list(id)%p%par%en_exchange = hc_params%en_exchange
    leaves%list(id)%p%par%en_exchange_tot = hc_params%en_exchange_tot
    leaves%list(id)%p%par%en_exchange_per_vol = hc_params%en_exchange_per_vol
    !
    if (a_disk_iter_params%flag_save_rates) then
      call save_chem_rates(id)
    end if
    !
    if (isnan(leaves%list(id)%p%par%pressure_thermal)) then
      write(str_disp, '(A, 4ES16.9)') 'Tgas,ngas,X_HI,X_H2:', &
        leaves%list(id)%p%par%Tgas, leaves%list(id)%p%par%n_gas, &
        leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI), &
        leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)
      call display_string_both(str_disp, a_book_keeping%fU)
    end if
    !
    if ((j .eq. 1) .and. a_disk_iter_params%deplete_oxygen_carbon .and. &
        (a_disk_iter_params%deplete_oxygen_carbon_tstart .gt. 0D0) .and. &
        (a_disk_iter_params%deplete_oxygen_carbon_tstart .lt. &
            leaves%list(id)%p%par%tmax_this)) then
      cycle
    end if
    !
    if ((leaves%list(id)%p%quality .eq. 0) .or. &
        (leaves%list(id)%p%par%t_final .ge. &
         0.5D0*leaves%list(id)%p%par%tmax_this)) then
      exit
    end if
    !
  end do
  !
  if (a_disk_ana_params%do_analyse) then
    if ((a_disk_iter_params%n_iter_used .gt. 0) .and. &
        (is_in_list_int(id, ana_ptlist%nlen, ana_ptlist%vals))) then
      call chem_analyse(id)
    end if
  end if
  !
  nullify(hc_params)
end subroutine calc_this_cell


subroutine update_en_exchange_with_dust
  integer i
  !
  hc_params%en_exchange_tot = cooling_gas_grain_collision() * hc_params%volume
  !
  ! This is the energy that the gas transfer to each type of dust per cell.
  ! Can be negative.
  if (a_disk%allow_gas_dust_en_exch) then
    do i=1, a_disk%ndustcompo
      hc_params%en_exchange(i) = &
        max(hc_params%en_exchange_per_vol(i) * hc_params%volume, &
            -frac_dust_lose_en * hc_params%en_gains(i))
        !
    end do
  end if
end subroutine update_en_exchange_with_dust



subroutine update_params_above_alt(i0)
  use load_Visser_CO_selfshielding
  integer, intent(in) :: i0
  integer i
  !
  ! Calculate the column density of a few species to the star and to the ISM
  do i=1, chem_idx_some_spe%nItem
    call calc_Ncol_to_ISM(leaves%list(i0)%p, i)
    call calc_Ncol_to_Star(leaves%list(i0)%p, i)
  end do
  associate(c => leaves%list(i0)%p)
    ! Kwok eq 10.20
    ! c%par%Av_toISM = 1.086D0 * c%par%ratioDust2HnucNum * &
    !   (phy_Pi * c%par%GrainRadius_CGS**2) * 2D0 * c%par%Ncol_toISM
    ! c%par%Av_toStar = 1.086D0 * c%par%ratioDust2HnucNum * &
    !   (phy_Pi * c%par%GrainRadius_CGS**2) * 2D0 * c%par%Ncol_toStar
    !
    c%par%f_selfshielding_toISM_H2  = min(1D0, get_H2_self_shielding( &
      c%col_den_toISM(chem_idx_some_spe%iiH2), c%par%velo_width_turb))
    c%par%f_selfshielding_toStar_H2  = min(1D0, get_H2_self_shielding( &
      c%col_den_toStar(chem_idx_some_spe%iiH2), c%par%velo_width_turb))
    !
    ! H2O and OH self shielding are already taken into account in the radiative transfer.
    ! Only for output; not used.
    c%par%f_selfshielding_toISM_H2O = &
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiH2O) * &
               const_LyAlpha_cross_H2O)))
    c%par%f_selfshielding_toStar_H2O = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiH2O) * &
               const_LyAlpha_cross_H2O)))
    !
    c%par%f_selfshielding_toISM_OH = &
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiOH) * &
               const_LyAlpha_cross_OH)))
    c%par%f_selfshielding_toStar_OH = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiOH) * &
               const_LyAlpha_cross_OH)))
    !
    c%par%f_selfshielding_toISM_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toISM(chem_idx_some_spe%iiH2), &
      c%col_den_toISM(chem_idx_some_spe%iiCO))))
    !
    c%par%f_selfshielding_toStar_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toStar(chem_idx_some_spe%iiH2), &
      c%col_den_toStar(chem_idx_some_spe%iiCO))))
    !
    c%par%zeta_Xray_H2 = calc_Xray_ionization_rate(c)
    !
    ! Calculate the gravitational force from above
    call calc_gravity_column(c)
    !
    if (.not. a_disk%waterShieldWithRadTran) then
      c%par%phflux_Lya = &
        c%par%flux_Lya_star_unatten / phy_LyAlpha_energy_CGS &
        * exp(-phy_UVext2Av * c%par%Av_toISM) &
        * c%par%f_selfshielding_toStar_H2O &
        * c%par%f_selfshielding_toStar_OH
      c%par%G0_Lya_atten = c%par%phflux_Lya / phy_Habing_photon_flux_CGS
    end if
  end associate
end subroutine update_params_above_alt



function get_H2_self_shielding(N_H2, dv_turb)
  ! Draine 1996, equation 37
  double precision get_H2_self_shielding
  double precision, intent(in) :: N_H2, dv_turb
  double precision x, b5, tmp
  x = N_H2 / 5D14
  b5 = dv_turb / 1D5
  tmp = sqrt(1D0 + x)
  get_H2_self_shielding = 0.965D0 / (1D0 + x/b5)**2 + &
    0.035 / tmp * exp(-8.5D-4 * tmp)
end function get_H2_self_shielding



subroutine check_convergency_cell(i0)
  integer, intent(in) :: i0
  ! Temperature is not considered.
  if (maxval(abs(leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                 - a_iter_stor%abundances(:, i0)) &
             - (a_disk_iter_params%atol_abun + &
                a_disk_iter_params%rtol_abun * &
                abs(leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                  + a_iter_stor%abundances(:, i0))) &
            ) .le. 0D0) then
    leaves%list(i0)%p%converged = .true.
  else
    leaves%list(i0)%p%converged = .false.
  end if
end subroutine check_convergency_cell



subroutine check_convergency_whole_disk
  integer i
  a_disk_iter_params%n_cell_converged = 0
  do i=1, leaves%nlen
    if (leaves%list(i)%p%converged) then
      a_disk_iter_params%n_cell_converged = a_disk_iter_params%n_cell_converged + 1
    end if
  end do
  a_disk_iter_params%flag_converged = &
    a_disk_iter_params%n_cell_converged .ge. &
    int(a_disk_iter_params%converged_cell_percentage_stop * real(leaves%nlen))
  write(str_disp, '("! Iter", I8, 4X, "Number of cells converged: ", I8, "/", I8)') &
    a_disk_iter_params%n_iter_used, a_disk_iter_params%n_cell_converged, leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
end subroutine check_convergency_whole_disk



subroutine update_calculating_cells
  integer, dimension(:), allocatable :: list_tmp
  integer i, i0, itmp, j, k, n
  logical flag_notyet
  allocate(list_tmp(n_calculating_cells_max))
  n = 0
  do i=1, n_calculating_cells
    i0 = calculating_cells(i)
    do j=1, leaves%list(i0)%p%below%n
      itmp = leaves%list(i0)%p%below%idx(j)
      flag_notyet = .true.
      do k=1, n
        if (list_tmp(k) .eq. itmp) then
          flag_notyet = .false.
          exit
        end if
      end do
      if (flag_notyet) then
        n = n + 1
        list_tmp(n) = itmp
      end if
    end do
  end do
  n_calculating_cells = n
  if (n .ge. 1) then
    calculating_cells(1:n) = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine update_calculating_cells



function calc_Xray_ionization_rate(c) result(z_Xray)
  ! Ionization rate per H atom
  use load_Bethell_Xray_cross
  double precision z_Xray
  type(type_cell), intent(in) :: c
  integer i, i1, i2
  double precision lam, en, sig, local_flux, rr, dlam
  double precision, parameter :: en_per_ion = 37D0 ! eV
  if (.not. allocated(c%optical)) then
    z_Xray = 0D0
    return
  end if
  if (.not. allocated(c%optical%flux)) then
    z_Xray = 0D0
    return
  end if
  i1 = max(1, get_idx_for_kappa(lam_range_Xray(1), dust_0)) + 1
  i2 = min(dust_0%n, get_idx_for_kappa(lam_range_Xray(2), dust_0)) - 1
  z_Xray = 0D0
  do i=i1, i2
    lam = dust_0%lam(i)
    en = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
         / (lam * 1D-8) / phy_eV2erg / 1D3 ! in keV
    sig = sigma_Xray_Bethell(en, c%par%dust_depletion, &
          c%par%ratioDust2HnucNum, c%par%GrainRadius_CGS)
    if (a_disk_iter_params%calc_zetaXray_from_Ncol) then
      rr = ((c%xmin + c%xmax)**2 + (c%ymin + c%ymax)**2) * 0.25D0 * phy_AU2cm**2
      if (i .lt. i2) then
        dlam = a_star%lam(i+1) - a_star%lam(i)
      else
        dlam = a_star%lam(i) - a_star%lam(i-1)
      end if
      local_flux = a_star%vals(i) * dlam * exp(-sig * c%par%Ncol_toStar) &
                   / (4D0 * phy_Pi * rr)
    else
      local_flux = c%optical%flux(i)
    end if
    z_Xray = z_Xray + &
             local_flux / (en * 1D3 * phy_eV2erg) * &
             sig * (en * 1D3 / en_per_ion)
  end do
end function calc_Xray_ionization_rate



subroutine set_initial_condition_4solver(id, j, iiter)
  integer, intent(in) :: id, j, iiter
  integer flag
  double precision t_min
  t_min = 1D2
  !
  flag = 0
  if (iiter .gt. 1) then
    flag = 10
  end if
  if (j .gt. 1) then
    flag = flag + 1
  end if
  select case(flag)
    case(0) ! iiter=1, j=1
      if (leaves%list(id)%p%above%n .gt. 0) then
        leaves%list(id)%p%par%Tgas = &
          leaves%list(leaves%list(id)%p%above%idx(1))%p%par%Tgas
      else
        if (leaves%list(id)%p%par%Tgas .le. 0D0) then
          leaves%list(id)%p%par%Tgas = &
            maxval(leaves%list(id)%p%par%Tdusts)*1.1D0 + 10D0
        end if
      end if
    case(1, 11) ! iiter=1, j>1
      if (leaves%list(id)%p%above%n .gt. 0) then
        leaves%list(id)%p%par%Tgas = &
          leaves%list(leaves%list(id)%p%above%idx(1))%p%par%Tgas
      else
        leaves%list(id)%p%par%Tgas = dble(j+1) * leaves%list(id)%p%par%Tdust + 10D0
      end if
    case(10) ! iiter>1, j=1
      ! Do nothing
  end select
  !
  if (isnan(leaves%list(id)%p%par%Tgas)) then
    leaves%list(id)%p%par%Tgas = dble(iiter+1) * leaves%list(id)%p%par%Tdust + 10D0
  end if
  !
  if (leaves%list(id)%p%par%Tgas .gt. a_disk_iter_params%Tgas_crazy) then
    leaves%list(id)%p%par%Tgas = a_disk_iter_params%Tgas_crazy
  end if
  !
  chemsol_stor%y(1:chem_species%nSpecies) = &
    chemsol_stor%y0(1:chem_species%nSpecies)
  !
  if (chem_idx_some_spe%i_Grain0 .ne. 0) then
    ! Initial abundance of *neutral* dust
    ! There should be no dust in the input initial abundance file
    chemsol_stor%y(chem_idx_some_spe%i_Grain0) = &
        leaves%list(id)%p%par%ratioDust2HnucNum
  end if
  !
  chemsol_stor%y(chem_species%nSpecies+1) = leaves%list(id)%p%par%Tgas
  !
  chemsol_params%evolT = .true.
  chemsol_params%maySwitchT = .true.
  if (leaves%list(id)%p%par%en_gain_tot .le. 0D0) then
    chemsol_params%evolT = .false.
  end if
  !
  chemsol_params%t0 = 0D0
  chemsol_params%dt_first_step = chemsol_params%dt_first_step0
  !
  if (a_disk_iter_params%use_fixed_tmax) then
    chemsol_params%t_max = chemsol_params%t_max0
  else
    chemsol_params%t_max = min(chemsol_params%t_max0, &
      max(t_min, a_disk_iter_params%nOrbit_tmax * &
      phy_2Pi/leaves%list(id)%p%par%omega_Kepler/phy_SecondsPerYear))
  end if
  leaves%list(id)%p%par%tmax_this = chemsol_params%t_max
  !
  if (a_disk_iter_params%deplete_oxygen_carbon .and. &
      (a_disk_iter_params%deplete_oxygen_carbon_tstart .le. 0D0)) then
    call deplete_oxygen_carbon_adhoc(id, chemsol_stor%y)
  end if
  if (a_disk_iter_params%deplete_oxygen_carbon .and. &
      (a_disk_iter_params%deplete_oxygen_carbon_tstart .gt. 0D0)) then
    chemsol_params%t_max = min(chemsol_params%t_max, &
                               a_disk_iter_params%deplete_oxygen_carbon_tstart)
  end if
  !
  call chem_evol_solve_prepare_ongoing
  !
end subroutine set_initial_condition_4solver



subroutine set_initial_condition_4solver_continue(id, j)
  integer, intent(in) :: id, j
  associate(c => leaves%list(id)%p)
    !
    !if ((c%par%t_final .ge. 2D0 * chemsol_params%t0) .and. &
    !    (abs(c%par%Tgas - chemsol_stor%y(chem_species%nSpecies+1)) .le. &
    !     1D-2*(c%par%Tgas + chemsol_stor%y(chem_species%nSpecies+1)))) then
    !  chemsol_params%evolT = .false.
    !  chemsol_params%maySwitchT = .true.
    !else
    chemsol_params%evolT = .true.
    chemsol_params%maySwitchT = .true.
    if (leaves%list(id)%p%par%en_gain_tot .le. 0D0) then
      chemsol_params%evolT = .false.
    end if
    !end if
    !
    chemsol_stor%y(1:chem_species%nSpecies) = c%abundances
    chemsol_stor%y(chem_species%nSpecies+1) = c%par%Tgas
    !if (j .ge. 3) then
    !  chemsol_stor%y(chem_species%nSpecies+1) = maxval(c%par%Tdusts)*1.1D0
    !end if
    !
    call rectify_abundances(chemsol_params%NEQ, chemsol_stor%y)
    !
    chemsol_params%t0 = c%par%t_final
    chemsol_params%dt_first_step = &
      max(chemsol_params%dt_first_step0, chemsol_params%t0*1D-3)
    !
    if ((j .gt. 1) .and. &
        a_disk_iter_params%deplete_oxygen_carbon .and. &
        (a_disk_iter_params%deplete_oxygen_carbon_tstart .gt. 0D0) .and. &
        (c%par%t_final .le. a_disk_iter_params%deplete_oxygen_carbon_tstart) .and. &
        (c%par%t_final .ge. 0.9D0*a_disk_iter_params%deplete_oxygen_carbon_tstart)) then
      !
      call deplete_oxygen_carbon_adhoc(id, chemsol_stor%y, 2)
      chemsol_params%t_max = c%par%tmax_this
      !
  end if
  end associate
  !
  call chem_evol_solve_prepare_ongoing
  !
end subroutine set_initial_condition_4solver_continue



subroutine deplete_oxygen_carbon_adhoc(id, y, flag)
  integer, intent(in) :: id
  integer, intent(in), optional :: flag
  double precision, intent(inout), dimension(:) :: y
  double precision r0, x_O, x_C, dep_O, dep_C, dep_N
  double precision X_O_0, X_C_0, X_N_0
  double precision Tgas, ngas
  integer, parameter :: iele_C = 7, iele_N = 8, iele_O = 9
  integer i, i0
  !
  r0 = (leaves%list(id)%p%xmin + leaves%list(id)%p%xmax) * 0.5D0
  x_O = r0 / a_disk_iter_params%r0_O
  x_C = r0 / a_disk_iter_params%r0_C
  !
  if (a_disk_iter_params%deplete_oxygen_carbon_method .eq. 'radial') then
    ! The degree of depletion is only a function of radius.
    dep_O = depl_f(x_O, a_disk_iter_params%a_O, a_disk_iter_params%b_O, &
                 a_disk_iter_params%gam_O)
    dep_C = depl_f(x_C, a_disk_iter_params%a_C, a_disk_iter_params%b_C, &
                 a_disk_iter_params%gam_C)
  else if (a_disk_iter_params%deplete_oxygen_carbon_method .eq. 'vscale') then
    ! The degree of depletion is a function of both r and z.
    if (a_disk_iter_params%deplete_oxygen_method .eq. 'table') then
      dep_O = depl_h(id, &
        depl_vfac_tab(r0, a_disk_iter_params%rmins_O, a_disk_iter_params%rmaxs_O, &
        a_disk_iter_params%vfacs_O, a_disk_iter_params%nvfacs_O), &
        a_disk_iter_params%gval_O)
    !
    else if (a_disk_iter_params%deplete_oxygen_method .eq. 'tanh') then
      dep_O = depl_h(id, &
        depl_vfac_tanh(r0, a_disk_iter_params%tanh_r_O, &
                       a_disk_iter_params%tanh_scale_O, &
                       a_disk_iter_params%tanh_minval_O, &
                       a_disk_iter_params%tanh_maxval_O &
                      ), &
        a_disk_iter_params%gval_O)
    !
    else
      dep_O = depl_h(id, &
          a_disk_iter_params%vfac_O * depl_vfac(x_O, a_disk_iter_params%p_O) &
          + a_disk_iter_params%k_O, &
          a_disk_iter_params%gval_O)
    end if
    !
    if (a_disk_iter_params%deplete_carbon_method .eq. 'table') then
      dep_C = depl_h(id, &
        depl_vfac_tab(r0, a_disk_iter_params%rmins_C, a_disk_iter_params%rmaxs_C, &
        a_disk_iter_params%vfacs_C, a_disk_iter_params%nvfacs_C), &
        a_disk_iter_params%gval_C)
    !
    else if (a_disk_iter_params%deplete_carbon_method .eq. 'tanh') then
      dep_C = depl_h(id, &
        depl_vfac_tanh(r0, a_disk_iter_params%tanh_r_C, &
                       a_disk_iter_params%tanh_scale_C, &
                       a_disk_iter_params%tanh_minval_C, &
                       a_disk_iter_params%tanh_maxval_C &
                      ), &
        a_disk_iter_params%gval_C)
    !
    else
      dep_C = depl_h(id, &
          a_disk_iter_params%vfac_C * depl_vfac(x_C, a_disk_iter_params%p_C) &
          + a_disk_iter_params%k_C, &
          a_disk_iter_params%gval_C)
    end if
  else if (a_disk_iter_params%deplete_oxygen_carbon_method .eq. 'vertical') then
    ! The deree of depletion is only a function of z.
    Tgas = leaves%list(id)%p%par%Tgas
    ngas = leaves%list(id)%p%par%n_gas
    dep_O = depl_g(chemsol_params%t_max, a_disk_iter_params%gval_O, &
        a_disk_iter_params%tads_O, &
        a_disk_iter_params%tsed_O, &
        a_disk_iter_params%r0_O, &
        a_disk_iter_params%k_O, &
        a_disk_iter_params%p_O, &
        Tgas, ngas, &
        r0, a_disk%star_mass_in_Msun)
    dep_C = depl_g(chemsol_params%t_max, a_disk_iter_params%gval_C, &
        a_disk_iter_params%tads_C, &
        a_disk_iter_params%tsed_C, &
        a_disk_iter_params%r0_C, &
        a_disk_iter_params%k_C, &
        a_disk_iter_params%p_C, &
        Tgas, ngas, &
        r0, a_disk%star_mass_in_Msun)
  else if (a_disk_iter_params%deplete_oxygen_carbon_method .eq. 'C/O-ratio') then
    ! First set the oxygen depletion factor, then set the carbon depletion
    ! factor with a prescribed distribution of the C/O ratio.
    dep_O = depl_h(id, &
        a_disk_iter_params%vfac_O * depl_vfac(x_O, a_disk_iter_params%p_O) &
        + a_disk_iter_params%k_O, &
        a_disk_iter_params%gval_O)
    if (abs(a_disk_iter_params%dep_zscale) .ge. 1D-10) then
      dep_C = min(1D0, dep_O * (1D0 + &
          a_disk_iter_params%O_to_C_ISM * &
          leaves%list(id)%p%ymin / a_disk_iter_params%dep_zscale))
    else
      dep_C = min(1D0, a_disk_iter_params%C_to_O_ratio * dep_O * a_disk_iter_params%O_to_C_ISM)
    end if
  else if (a_disk_iter_params%deplete_oxygen_carbon_method .eq. 'uniform') then
    dep_O = a_disk_iter_params%f_depl_O
    dep_C = a_disk_iter_params%f_depl_C
  end if
  !
  if (r0 .le. a_disk_iter_params%rin_O) then
    dep_O = dep_O * a_disk_iter_params%fin_O
  end if
  if (r0 .le. a_disk_iter_params%rin_C) then
    dep_C = dep_C * a_disk_iter_params%fin_C
  end if
  if (r0 .ge. a_disk_iter_params%rout_O) then
    dep_O = dep_O * a_disk_iter_params%fout_O
  end if
  if (r0 .ge. a_disk_iter_params%rout_C) then
    dep_C = dep_C * a_disk_iter_params%fout_C
  end if
  !
  dep_O = min(dep_O, a_disk_iter_params%tanh_OC_enhance_max_O)
  dep_C = min(dep_C, a_disk_iter_params%tanh_OC_enhance_max_C)
  !
  if (  (abs(dep_O - 1D0) .le. 1D-3) .and. &
        (abs(dep_C - 1D0) .le. 1D-3)) then
    ! Do nothing if no depletion to be applied
    return
  end if
  !
  if (a_disk_iter_params%deplete_nitrogen .and. a_disk_iter_params%deplete_nitrogen_as_carbon) then
    dep_N = dep_C
  else
    dep_N = 1D0
  end if
  !
  if (.not. present(flag)) then
    !y(chem_idx_some_spe%i_gH2O) = 3.2D-4 * dep_O
    !y(chem_idx_some_spe%i_CO)   = 0D0
    !y(chem_idx_some_spe%i_CI)   = 1.4D-4 * dep_C
    X_O_0 = y(chem_idx_some_spe%i_gH2O) + y(chem_idx_some_spe%i_H2O) + y(chem_idx_some_spe%i_OI) + y(chem_idx_some_spe%i_CO)
    X_C_0 = y(chem_idx_some_spe%i_CO) + y(chem_idx_some_spe%i_CI) + y(chem_idx_some_spe%i_CII)
    X_N_0 = y(chem_idx_some_spe%i_NI)
    !
    y(chem_idx_some_spe%i_gH2O) = X_O_0 * dep_O / 3D0 ! X_O_0 * dep_O  !1.8D-4 * dep_O
    y(chem_idx_some_spe%i_H2O)  = X_O_0 * dep_O / 3D0 ! X_O_0 * dep_O  !1.8D-4 * dep_O
    y(chem_idx_some_spe%i_CO)   = min(X_O_0 * dep_O / 3D0, X_C_0 * dep_C) ! 0D0            !1.4D-4 * min(dep_O, dep_C)
    y(chem_idx_some_spe%i_CI) = max(0D0, X_C_0 * dep_C - y(chem_idx_some_spe%i_CO))  ! X_C_0 * dep_C  !max(0D0, 1.4D-4 * dep_C - y(chem_idx_some_spe%i_CO))
    y(chem_idx_some_spe%i_NI) = X_N_0 * dep_N  !max(0D0, 7.5D-5 * dep_N)
  else if (flag .eq. 1) then
    y(chem_idx_some_spe%i_gH2O) = y(chem_idx_some_spe%i_gH2O) * dep_O
    y(chem_idx_some_spe%i_H2O) = y(chem_idx_some_spe%i_H2O) * dep_O
    y(chem_idx_some_spe%i_OI) = y(chem_idx_some_spe%i_OI) * dep_O
    y(chem_idx_some_spe%i_gCO)  = y(chem_idx_some_spe%i_gCO) * dep_C
    y(chem_idx_some_spe%i_CO)  = y(chem_idx_some_spe%i_CO) * dep_C
    y(chem_idx_some_spe%i_gCO2) = y(chem_idx_some_spe%i_gCO2) * dep_C
    y(chem_idx_some_spe%i_CI)  = y(chem_idx_some_spe%i_CI) * dep_C
    y(chem_idx_some_spe%i_CII)  = y(chem_idx_some_spe%i_CII) * dep_C
  else
    chemsol_stor%y = y
    call chem_elemental_residence
    do i=1, chem_ele_resi(iele_C)%n_nonzero
      i0 = chem_ele_resi(iele_C)%iSpecies(i)
      if (chem_species%elements(iele_O, i0) .gt. 0) then
        y(i0) = y(i0) * min(dep_C, dep_O)
      else
        y(i0) = y(i0) * dep_C
      end if
    end do
    do i=1, chem_ele_resi(iele_O)%n_nonzero
      i0 = chem_ele_resi(iele_O)%iSpecies(i)
      if (chem_species%elements(iele_C, i0) .eq. 0) then
        y(i0) = y(i0) * dep_O
      end if
    end do
    do i=1, chem_ele_resi(iele_N)%n_nonzero
      i0 = chem_ele_resi(iele_N)%iSpecies(i)
      if ((chem_species%elements(iele_C, i0) .eq. 0) .and. &
          (chem_species%elements(iele_O, i0) .eq. 0)) then
        y(i0) = y(i0) * dep_N
      end if
    end do
  end if
end subroutine deplete_oxygen_carbon_adhoc



function depl_f(x, a, b, gam)
  double precision depl_f
  double precision, intent(in) :: x, a, b, gam
  depl_f = (x**gam * a + b) / (x**gam + 1D0)
end function depl_f


function depl_g(t_evol, ground_val, t0_ads, t0_sed, r0, k, p, &
                Tgas, n_gas, RtoStar_AU, MstarInMsun)
  ! t0_ads = 1D2 yr
  ! t0_sed = 1D5 yr
  double precision depl_g
  double precision, intent(in) :: &
    ground_val, t0_ads, t0_sed, r0, k, p, t_evol, Tgas, n_gas, RtoStar_AU, MstarInMsun
  double precision t_ads, t_sed, tmp
  tmp = sqrt(Tgas/1D2) * (n_gas/1D7)
  t_ads = t0_ads / tmp
  t_sed = t0_sed * (RtoStar_AU/1D2)**3 / (MstarInMsun) * tmp
  depl_g = ground_val + 1D0/(k + (RtoStar_AU/r0)**p) * exp(-t_evol / (t_ads + t_sed))
end function depl_g


function depl_h(id, vfac, gval)
  double precision depl_h
  integer, intent(in) :: id
  double precision, intent(in) :: vfac, gval
  double precision vscal_factor
  integer icol
  logical found
  type(type_cell), pointer :: cthis
  found = .false.
  do icol=1, n_columns
    cthis => leaves%list(bott_cells%idx(icol))%p
    if ((abs(cthis%xmin - leaves%list(id)%p%xmin) .le. &
         1D-3*(cthis%xmin + leaves%list(id)%p%xmin)) .and. &
        (abs(cthis%xmax - leaves%list(id)%p%xmax) .le. &
         1D-3*(cthis%xmax + leaves%list(id)%p%xmax))) then
      found = .true.
      exit
    end if
  end do
  if (.not. found) then
    write(*, '(A)') 'In depl_h:'
    write(*, '(A)') 'Cannot find column index!'
    write(*, '(A, 4ES14.4)') 'xmin,xmax,ymin,ymax = ', &
        leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
        leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
    call error_stop()
  end if
  vscal_factor = leaves%list(id)%p%par%n_gas / cthis%par%n_gas
  depl_h = vscal_factor**vfac + gval
end function depl_h


function depl_vfac(x, p)
  double precision depl_vfac
  double precision, intent(in) :: x, p
  double precision tmp
  tmp = x**p
  depl_vfac = tmp / (1D0 + tmp)
end function depl_vfac


function depl_vfac_tanh(x, xshift, xscale, minv, maxv)
    double precision depl_vfac_tanh
    double precision, intent(in) :: x, xshift, xscale, minv, maxv
    double precision y
    double precision, parameter :: ymin = -1D0, ymax = 1D0
    y = -tanh((x-xshift)/xscale)
    y = (y - ymin) * ((maxv - minv) / (ymax - ymin)) + minv
    depl_vfac_tanh = 1D0/(y*y) - 1D0
end function depl_vfac_tanh


function depl_vfac_tab(r, rmins, rmaxs, vs, n)
  double precision depl_vfac_tab
  integer, intent(in) :: n
  double precision, intent(in) :: r
  double precision, dimension(:), intent(in) :: rmins, rmaxs, vs
  integer i
  do i=1, n
    if ((rmins(i) .le. r) .and. (r .le. rmaxs(i))) then
        depl_vfac_tab = vs(i)
        return
    end if
  end do
  depl_vfac_tab = 0D0
end function depl_vfac_tab


subroutine deallocate_columns
  integer i, j
#ifdef DIAGNOSIS_TRACK_FUNC_CALL
  write(*,*) 'Running deallocate_columns.'
#endif
  if (allocated(columns)) then
    do i=1, n_columns
      if (allocated(columns(i)%list)) then
        do j=1, columns(i)%nlen
          nullify(columns(i)%list(j)%p)
        end do
        deallocate(columns(i)%list)
        if (allocated(columns_idx(i)%vals)) then
          deallocate(columns_idx(i)%vals)
        end if
      end if
    end do
    deallocate(columns, columns_idx)
  end if
end subroutine deallocate_columns



subroutine make_columns
  integer i, j, i1
  type(type_cell), pointer :: cthis, cnext
  double precision length, r, z, eps
  integer dirtype, n_using
  logical found
  type(type_ray) ray
  !
#ifdef DIAGNOSIS_TRACK_FUNC_CALL
  write(*,*) 'Running make_columns.'
#endif
  n_columns = bott_cells%nlen
  allocate(columns(n_columns), columns_idx(n_columns))
  do i=1, n_columns
    cthis => leaves%list(bott_cells%idx(i))%p
    r = (cthis%xmin + cthis%xmax) * 0.5D0
    z = root%ymax * 2D0
    columns(i)%nlen = int(calc_Ncol_from_cell_to_point(cthis, r, z, -2))
    if (.not. allocated(columns(i)%list)) then
      allocate(columns(i)%list(columns(i)%nlen))
    end if
    !
    ray%x = (cthis%xmin + cthis%xmax) * 0.5D0
    ray%y = 0D0
    ray%z = (cthis%ymin + cthis%ymax) * 0.5D0
    ray%vx = 0D0
    ray%vy = 0D0
    ray%vz = 1D0
    !
    n_using = 0
    j = 0
    ! First make a list of all the cells (including null ones) above the bottom
    ! (mid-plane) one.  They are ordered from bottom to top.
    do
      call calc_intersection_ray_cell(ray, cthis, &
        length, r, z, eps, found, dirtype)
      if (.not. found) then
        write(str_disp,'(A)') 'In make_columns:'
        call display_string_both(str_disp, a_book_keeping%fU)
        write(str_disp,'(A, 9ES16.6)') 'ray does not intersect cthis: ', &
          sqrt(ray%x**2+ray%y**2), ray%z, ray%vx, ray%vy, ray%vz, &
            cthis%xmin, cthis%xmax, cthis%ymin, cthis%ymax
        call display_string_both(str_disp, a_book_keeping%fU)
        write(*, *)
        return
      end if
      !
      j = j + 1
      columns(i)%list(j)%p => cthis
      !
      if (columns(i)%list(j)%p%using) then
        n_using = n_using + 1
      end if
      !
      ray%x = ray%x + ray%vx * (length + eps)
      ray%y = ray%y + ray%vy * (length + eps)
      ray%z = ray%z + ray%vz * (length + eps)
      !
      call locate_photon_cell_alt(r, z, cthis, dirtype, cnext, found)
      if (found) then
        cthis => cnext
      else
        exit
      end if
    end do
    !
    columns_idx(i)%nlen = n_using
    allocate(columns_idx(i)%vals(n_using))
    i1 = n_using + 1
    do j=1, columns(i)%nlen
      if (columns(i)%list(j)%p%using) then
         i1 = i1 - 1
        columns_idx(i)%vals(i1) = columns(i)%list(j)%p%id
      end if
    end do
  end do
end subroutine make_columns



subroutine calc_Ncol_to_ISM(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout) :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toISM(iSp) = calc_Ncol_from_cell_to_point( &
      c, (c%xmin+c%xmax)*0.5D0, root%ymax * 2D0, &
      chem_idx_some_spe%idx(iSp), fromCellCenter=.false.)
  else
    c%par%Ncol_toISM = calc_Ncol_from_cell_to_point( &
      c, (c%xmin+c%xmax)*0.5D0, root%ymax * 2D0, fromCellCenter=.false.)
  end if
end subroutine calc_Ncol_to_ISM



subroutine calc_Ncol_to_Star(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout) :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toStar(iSp) = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0, chem_idx_some_spe%idx(iSp), fromCellCenter=.false.)
  else
    c%par%Ncol_toStar = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0, fromCellCenter=.false.)
  end if
end subroutine calc_Ncol_to_Star



function calc_Ncol_from_cell_to_point(c, r, z, iSpe, fromCellCenter) result(N)
  double precision N
  type(type_cell), intent(in), target :: c ! Todo: pointer or not?
  double precision, intent(in) :: r, z
  integer, intent(in), optional :: iSpe
  logical, intent(in), optional :: fromCellCenter
  logical fromCC
  type(type_ray) ray
  type(type_cell), pointer :: cthis, cnext
  double precision t, length, r1, z1, eps, dx, dy
  logical found
  integer dirtype, iloc(1)
  double precision, parameter :: small_dist = 1D-50, small_frac = 1D-6
  !
  if (present(fromCellCenter)) then
    fromCC = fromCellCenter
  else
    fromCC = .true.
  end if
  !
  N = 0D0
  !
  if (fromCC) then
    ray%x = (c%xmin+c%xmax)*0.5D0
    ray%y = 0D0
    ray%z = (c%ymin+c%ymax)*0.5D0
  else
    iloc = minloc( &
      (/ (r-c%xmin)**2 + (z-c%ymin)**2, &
         (r-c%xmin)**2 + (z-c%ymax)**2, &
         (r-c%xmax)**2 + (z-c%ymin)**2, &
         (r-c%xmax)**2 + (z-c%ymax)**2, &
         (r-0.5D0*(c%xmin+c%xmax))**2 + (z-0.5D0*(c%ymin+c%ymax))**2 &
      /))
    dx = c%xmax - c%xmin
    dy = c%ymax - c%ymin
    select case(iloc(1))
      case(1)
        ray%x = c%xmin + dx * small_frac
        ray%z = c%ymin + dy * small_frac
      case(2)
        ray%x = c%xmin + dx * small_frac
        ray%z = c%ymax - dy * small_frac
      case(3)
        ray%x = c%xmax - dx * small_frac
        ray%z = c%ymin + dy * small_frac
      case(4)
        ray%x = c%xmax - dx * small_frac
        ray%z = c%ymax - dy * small_frac
      case(5)
        ray%x = (c%xmin+c%xmax)*0.5D0
        ray%z = (c%ymin+c%ymax)*0.5D0
    end select
    ray%y = 0D0
  end if
  !
  ray%vx = r - ray%x
  ray%vy = 0D0
  ray%vz = z - ray%z
  t = sqrt(ray%vx**2 + ray%vy**2 + ray%vz**2)
  if (t .le. small_dist) then
    return
  end if
  ray%vx = ray%vx / t
  ray%vy = ray%vy / t
  ray%vz = ray%vz / t
  !
  cthis => c
  !
  do
    call calc_intersection_ray_cell(ray, cthis, &
        length, r1, z1, eps, found, dirtype)
    if (.not. found) then
      write(str_disp,'(A, 9ES13.4)') &
        'In calc_Ncol_from_cell_to_point: ray does not intersect cthis: ', &
        sqrt(ray%x**2+ray%y**2), ray%z, ray%vx, ray%vy, ray%vz, &
            cthis%xmin, cthis%xmax, cthis%ymin, cthis%ymax
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    end if
    !
    if (present(iSpe)) then
      if ((iSpe .ge. 1) .and. (iSpe .le. chem_species%nSpecies)) then
        ! Calculate column density of a species
        if (cthis%using) then
          N = N + cthis%par%n_gas * cthis%abundances(iSpe) * &
              (length * phy_AU2cm)
        end if
      !
      else if (iSpe .eq. 0) then
        ! Calculate column density
        if (cthis%using) then
          N = N + cthis%par%n_gas * length * phy_AU2cm
        end if
      !
      else if (iSpe .eq. -1) then
        ! Calculate gravity force (accumulated), not including the gravity of
        ! the starting cell.
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          if (cthis%using) then
            N = N + cthis%par%gravity_z
          end if
        end if
      !
      else if (iSpe .eq. -2) then
        ! Calculate the number of cells crossed by the ray.
        N = N + 1D0
      !
      else if (iSpe .eq. -3) then
        ! Calculate column density, including null cells.
        ! Density is taken from the initial value
        N = N + cthis%val(1) * length * phy_AU2cm
      !
      else if (iSpe .eq. -4) then
        ! Calculate column density, including null cells, but not itself.
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          N = N + cthis%val(1) * length * phy_AU2cm
        end if
      !
      else if (iSpe .eq. -5) then
        ! Calculate the dust column density
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          if (cthis%using) then
            N = N + cthis%par%ndust_tot * length * phy_AU2cm
          end if
        end if
      else if (iSpe .eq. -6) then
        ! Calculate the dust column density
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          if (cthis%using) then
            N = N + cthis%par%ndust_tot * &
                    (phy_Pi * c%par%GrainRadius_CGS**2 * 2D0) * &
                    length * phy_AU2cm
          end if
        end if
      else
        write(str_disp,'(A)') 'I do not know what to do!'
        call display_string_both(str_disp, a_book_keeping%fU)
        call error_stop()
      end if
    else
      ! Calculate column density by default
      if (cthis%using) then
        N = N + cthis%par%n_gas * length * phy_AU2cm
      end if
    end if
    !
    ray%x = ray%x + ray%vx * (length + eps)
    ray%y = ray%y + ray%vy * (length + eps)
    ray%z = ray%z + ray%vz * (length + eps)
    !
    call locate_photon_cell_alt(r1, z1, cthis, dirtype, cnext, found)
    if (found) then
      cthis => cnext
    else
      exit
    end if
  end do
  !
  if (present(iSpe)) then
    if (iSpe .eq. -1) then
      if (cthis%using) then
        N = N + cthis%par%gravity_acc_z
      end if
    end if
  end if
end function calc_Ncol_from_cell_to_point



subroutine disk_save_results_pre
  write(filename_save_results, '("iter_", I4.4, ".dat")') &
    a_disk_iter_params%n_iter_used
  filename_save_results = trim(combine_dir_filename( &
    a_disk_iter_params%iter_files_dir, filename_save_results))
  call openFileSequentialWrite(fU_save_results, filename_save_results, &
       99999, getu=1)
  !
  call write_header(fU_save_results)
end subroutine disk_save_results_pre


subroutine write_header(fU)
  integer, intent(in) :: fU
  character(len=64) fmt_str
  character(len=10000) tmp_str
  write(fmt_str, '("(", I4, "A14)")') chem_species%nSpecies
  write(tmp_str, fmt_str) chem_species%names
  write(fU, '(A)') &
    '!' // &
    str_pad_to_len('cvg',  4) // &
    str_pad_to_len('qual', 5) // &
    str_pad_to_len('cr_count',len_item) // &
    str_pad_to_len('abc_dus', len_item) // &
    str_pad_to_len('scc_HI',  len_item) // &
    str_pad_to_len('abc_wat', len_item) // &
    str_pad_to_len('t_final', len_item) // &
    str_pad_to_len('rmin',    len_item) // &
    str_pad_to_len('rmax',    len_item) // &
    str_pad_to_len('zmin',    len_item) // &
    str_pad_to_len('zmax',    len_item) // &
    str_pad_to_len('n_gas',   len_item) // &
    str_pad_to_len('Tgas',    len_item) // &
    str_pad_to_len('Tdust',   len_item) // &
    str_pad_to_len('Tdust1',  len_item) // &
    str_pad_to_len('Tdust2',  len_item) // &
    str_pad_to_len('Tdust3',  len_item) // &
    str_pad_to_len('Tdust4',  len_item) // &
    str_pad_to_len('ndust_1', len_item) // &
    str_pad_to_len('ndust_2', len_item) // &
    str_pad_to_len('ndust_3', len_item) // &
    str_pad_to_len('ndust_4', len_item) // &
    str_pad_to_len('ndust_t', len_item) // &
    str_pad_to_len('rhodus_1', len_item) // &
    str_pad_to_len('rhodus_2', len_item) // &
    str_pad_to_len('rhodus_3', len_item) // &
    str_pad_to_len('rhodus_4', len_item) // &
    str_pad_to_len('sigdus_1', len_item) // &
    str_pad_to_len('sigdus_2', len_item) // &
    str_pad_to_len('sigdus_3', len_item) // &
    str_pad_to_len('sigdus_4', len_item) // &
    str_pad_to_len('sigd_av', len_item) // &
    str_pad_to_len('d2gmas',  len_item) // &
    str_pad_to_len('d2gnum',  len_item) // &
    str_pad_to_len('deplet',  len_item) // &
    str_pad_to_len('mg_cell', len_item) // &
    str_pad_to_len('md_cell', len_item) // &
    str_pad_to_len('presr_t', len_item) // &
    str_pad_to_len('presr_g', len_item) // &
    str_pad_to_len('egain_d', len_item) // &
    str_pad_to_len('egain_ab',len_item) // &
    str_pad_to_len('egain_e', len_item) // &
    str_pad_to_len('egain_d1', len_item) // &
    str_pad_to_len('egain_e1', len_item) // &
    str_pad_to_len('egain_d2', len_item) // &
    str_pad_to_len('egain_e2', len_item) // &
    str_pad_to_len('egain_d3', len_item) // &
    str_pad_to_len('egain_e3', len_item) // &
    str_pad_to_len('egain_d4', len_item) // &
    str_pad_to_len('egain_e4', len_item) // &
    str_pad_to_len('flx_tot', len_item) // &
    str_pad_to_len('flx_Xray',len_item) // &
    str_pad_to_len('G0_UV',   len_item) // &
    str_pad_to_len('flx_Lya', len_item) // &
    str_pad_to_len('flx_Vis', len_item) // &
    str_pad_to_len('flx_NIR', len_item) // &
    str_pad_to_len('flx_MIR', len_item) // &
    str_pad_to_len('flx_FIR', len_item) // &
    str_pad_to_len('vr_tot',  len_item) // &
    str_pad_to_len('vz_tot',  len_item) // &
    str_pad_to_len('ani_tot', len_item) // &
    str_pad_to_len('vr_Xray', len_item) // &
    str_pad_to_len('vz_Xray', len_item) // &
    str_pad_to_len('ani_Xray',len_item) // &
    str_pad_to_len('vr_UV',   len_item) // &
    str_pad_to_len('vz_UV',   len_item) // &
    str_pad_to_len('ani_UV',  len_item) // &
    str_pad_to_len('vr_Lya',  len_item) // &
    str_pad_to_len('vz_Lya',  len_item) // &
    str_pad_to_len('ani_Lya', len_item) // &
    str_pad_to_len('vr_Vis',  len_item) // &
    str_pad_to_len('vz_Vis',  len_item) // &
    str_pad_to_len('ani_Vis', len_item) // &
    str_pad_to_len('vr_NIR',  len_item) // &
    str_pad_to_len('vz_NIR',  len_item) // &
    str_pad_to_len('ani_NIR', len_item) // &
    str_pad_to_len('vr_MIR',  len_item) // &
    str_pad_to_len('vz_MIR',  len_item) // &
    str_pad_to_len('ani_MIR', len_item) // &
    str_pad_to_len('vr_FIR',  len_item) // &
    str_pad_to_len('vz_FIR',  len_item) // &
    str_pad_to_len('ani_FIR', len_item) // &
    str_pad_to_len('Av_ISM',  len_item) // &
    str_pad_to_len('Av_Star', len_item) // &
    str_pad_to_len('UV_G0_I', len_item) // &
    str_pad_to_len('UV_G0_S', len_item) // &
    str_pad_to_len('LyAG0_a', len_item) // &
    str_pad_to_len('LyANF0',  len_item) // &
    str_pad_to_len('zeta_X',  len_item) // &
    str_pad_to_len('Ncol_I',  len_item) // &
    str_pad_to_len('Ncol_S',  len_item) // &
    str_pad_to_len('N_H2_I',  len_item) // &
    str_pad_to_len('N_H2O_I', len_item) // &
    str_pad_to_len('N_OH_I',  len_item) // &
    str_pad_to_len('N_CO_I',  len_item) // &
    str_pad_to_len('N_H2_S',  len_item) // &
    str_pad_to_len('N_H2O_S', len_item) // &
    str_pad_to_len('N_OH_S',  len_item) // &
    str_pad_to_len('N_CO_S',  len_item) // &
    str_pad_to_len('f_H2_I',  len_item) // &
    str_pad_to_len('f_H2O_I', len_item) // &
    str_pad_to_len('f_OH_I',  len_item) // &
    str_pad_to_len('f_CO_I',  len_item) // &
    str_pad_to_len('f_H2_S',  len_item) // &
    str_pad_to_len('f_H2O_S', len_item) // &
    str_pad_to_len('f_OH_S',  len_item) // &
    str_pad_to_len('f_CO_S',  len_item) // &
    str_pad_to_len('R_H2_fo', len_item) // &
    str_pad_to_len('hc_net',  len_item) // &
    str_pad_to_len('h_ph_gr', len_item) // &
    str_pad_to_len('h_fo_H2', len_item) // &
    str_pad_to_len('h_cosmi', len_item) // &
    str_pad_to_len('h_vi_H2', len_item) // &
    str_pad_to_len('h_io_CI', len_item) // &
    str_pad_to_len('h_ph_H2', len_item) // &
    str_pad_to_len('h_ph_wa', len_item) // &
    str_pad_to_len('h_ph_OH', len_item) // &
    str_pad_to_len('h_Xray ', len_item) // &
    str_pad_to_len('h_visco', len_item) // &
    str_pad_to_len('h_chem',  len_item) // &
    str_pad_to_len('c_el_gr', len_item) // &
    str_pad_to_len('c_vi_H2', len_item) // &
    str_pad_to_len('c_gg_co', len_item) // &
    str_pad_to_len('c_OI   ', len_item) // &
    str_pad_to_len('c_CII  ', len_item) // &
    str_pad_to_len('c_NII  ', len_item) // &
    str_pad_to_len('c_SiII ', len_item) // &
    str_pad_to_len('c_FeII ', len_item) // &
    str_pad_to_len('c_OH_ro', len_item) // &
    str_pad_to_len('c_wa_ro', len_item) // &
    str_pad_to_len('c_wa_vi', len_item) // &
    str_pad_to_len('c_CO_ro', len_item) // &
    str_pad_to_len('c_CO_vi', len_item) // &
    str_pad_to_len('c_H2_ro', len_item) // &
    str_pad_to_len('c_LyAlp', len_item) // &
    str_pad_to_len('c_fb   ', len_item) // &
    str_pad_to_len('c_ff   ', len_item) // &
    str_pad_to_len('alpha  ', len_item) // &
    str_pad_to_len('am     ', len_item) // &
    str_pad_to_len('ion_cha', len_item) // &
    str_pad_to_len('v_Kep',   len_item) // &
    str_pad_to_len('w_Kep',   len_item) // &
    str_pad_to_len('dv_dr',   len_item) // &
    str_pad_to_len('c_sound', len_item) // &
    str_pad_to_len('dv_turb', len_item) // &
    str_pad_to_len('l_coher', len_item) // &
    str_pad_to_len('nsit_gr', len_item) // &
    str_pad_to_len('nmol_gr', len_item) // &
    trim(tmp_str)
end subroutine write_header


subroutine disk_save_results_write(fU, c)
  character(len=64) fmt_str
  integer, intent(in) :: fU
  type(type_cell), pointer, intent(in) :: c
  integer converged, crct
  !
  write(fmt_str, '(", ", I4, "ES14.5E3)")') chem_species%nSpecies
  if (c%converged) then
    converged = 1
  else
    converged = 0
  end if
  !
  if (allocated(c%optical)) then
    crct = c%optical%cr_count
  else
    crct = 0
  end if
  write(fU, '(2I5, 4I14, 142ES14.5E3' // trim(fmt_str)) &
  converged                                              , &
  c%quality                                              , &
  crct                                                   , &
  c%par%ab_count_dust                                    , &
  c%par%sc_count_HI                                      , &
  c%par%ab_count_water                                   , &
  c%par%t_final                                          , &
  c%xmin                                                 , &
  c%xmax                                                 , &
  c%ymin                                                 , &
  c%ymax                                                 , &
  c%par%n_gas                                            , &
  c%par%Tgas                                             , &
  c%par%Tdust                                            , &
  c%par%Tdusts(1)                                        , &
  c%par%Tdusts(2)                                        , &
  c%par%Tdusts(3)                                        , &
  c%par%Tdusts(4)                                        , &
  c%par%n_dusts(1)                                       , &
  c%par%n_dusts(2)                                       , &
  c%par%n_dusts(3)                                       , &
  c%par%n_dusts(4)                                       , &
  c%par%ndust_tot                                        , &
  c%par%rho_dusts(1)                                     , &
  c%par%rho_dusts(2)                                     , &
  c%par%rho_dusts(3)                                     , &
  c%par%rho_dusts(4)                                     , &
  c%par%sig_dusts(1)                                     , &
  c%par%sig_dusts(2)                                     , &
  c%par%sig_dusts(3)                                     , &
  c%par%sig_dusts(4)                                     , &
  c%par%sigdust_ave                                      , &
  c%par%ratioDust2GasMass                                , &
  c%par%ratioDust2HnucNum                                , &
  c%par%dust_depletion                                   , &
  c%par%mgas_cell                                        , &
  c%par%mdust_tot                                        , &
  c%par%pressure_thermal                                 , &
  c%par%gravity_acc_z / c%par%area_T                     , &
  c%par%en_gain_tot                                      , &
  c%par%en_gain_abso_tot                                 , &
  c%par%en_exchange_tot                                  , &
  c%par%en_gains(1)                                      , &
  c%par%en_exchange(1)                                   , &
  c%par%en_gains(2)                                      , &
  c%par%en_exchange(2)                                   , &
  c%par%en_gains(3)                                      , &
  c%par%en_exchange(3)                                   , &
  c%par%en_gains(4)                                      , &
  c%par%en_exchange(4)                                   , &
  c%par%flux_tot                                         , &
  c%par%flux_Xray                                        , &
  c%par%flux_UV/phy_Habing_energy_flux_CGS               , &
  c%par%flux_Lya                                         , &
  c%par%flux_Vis                                         , &
  c%par%flux_NIR                                         , &
  c%par%flux_MIR                                         , &
  c%par%flux_FIR                                         , &
  c%par%dir_tot_r                                        , &
  c%par%dir_tot_z                                        , &
  c%par%aniso_tot                                        , &
  c%par%dir_Xray_r                                       , &
  c%par%dir_Xray_z                                       , &
  c%par%aniso_Xray                                       , &
  c%par%dir_UV_r                                         , &
  c%par%dir_UV_z                                         , &
  c%par%aniso_UV                                         , &
  c%par%dir_Lya_r                                        , &
  c%par%dir_Lya_z                                        , &
  c%par%aniso_Lya                                        , &
  c%par%dir_Vis_r                                        , &
  c%par%dir_Vis_z                                        , &
  c%par%aniso_Vis                                        , &
  c%par%dir_NIR_r                                        , &
  c%par%dir_NIR_z                                        , &
  c%par%aniso_NIR                                        , &
  c%par%dir_MIR_r                                        , &
  c%par%dir_MIR_z                                        , &
  c%par%aniso_MIR                                        , &
  c%par%dir_FIR_r                                        , &
  c%par%dir_FIR_z                                        , &
  c%par%aniso_FIR                                        , &
  c%par%Av_toISM                                         , &
  c%par%Av_toStar                                        , &
  c%par%G0_UV_toISM                                      , &
  c%par%G0_UV_toStar                                     , &
  c%par%G0_Lya_atten                                     , &
  c%par%phflux_Lya                                       , &
  c%par%zeta_Xray_H2                                     , &
  c%par%Ncol_toISM                                       , &
  c%par%Ncol_toStar                                      , &
  c%col_den_toISM(chem_idx_some_spe%iiH2)                , &
  c%col_den_toISM(chem_idx_some_spe%iiH2O)               , &
  c%col_den_toISM(chem_idx_some_spe%iiOH)                , &
  c%col_den_toISM(chem_idx_some_spe%iiCO)                , &
  c%col_den_toStar(chem_idx_some_spe%iiH2)               , &
  c%col_den_toStar(chem_idx_some_spe%iiH2O)              , &
  c%col_den_toStar(chem_idx_some_spe%iiOH)               , &
  c%col_den_toStar(chem_idx_some_spe%iiCO)               , &
  c%par%f_selfshielding_toISM_H2                         , &
  c%par%f_selfshielding_toISM_H2O                        , &
  c%par%f_selfshielding_toISM_OH                         , &
  c%par%f_selfshielding_toISM_CO                         , &
  c%par%f_selfshielding_toStar_H2                        , &
  c%par%f_selfshielding_toStar_H2O                       , &
  c%par%f_selfshielding_toStar_OH                        , &
  c%par%f_selfshielding_toStar_CO                        , &
  c%par%R_H2_form_rate                                   , &
  c%h_c_rates%hc_net_rate                                , &
  c%h_c_rates%heating_photoelectric_small_grain_rate     , &
  c%h_c_rates%heating_formation_H2_rate                  , &
  c%h_c_rates%heating_cosmic_ray_rate                    , &
  c%h_c_rates%heating_vibrational_H2_rate                , &
  c%h_c_rates%heating_ionization_CI_rate                 , &
  c%h_c_rates%heating_photodissociation_H2_rate          , &
  c%h_c_rates%heating_photodissociation_H2O_rate         , &
  c%h_c_rates%heating_photodissociation_OH_rate          , &
  c%h_c_rates%heating_Xray_Bethell_rate                  , &
  c%h_c_rates%heating_viscosity_rate                     , &
  c%h_c_rates%heating_chem                               , &
  c%h_c_rates%cooling_photoelectric_small_grain_rate     , &
  c%h_c_rates%cooling_vibrational_H2_rate                , &
  c%h_c_rates%cooling_gas_grain_collision_rate           , &
  c%h_c_rates%cooling_OI_rate                            , &
  c%h_c_rates%cooling_CII_rate                           , &
  c%h_c_rates%cooling_NII_rate                           , &
  c%h_c_rates%cooling_SiII_rate                          , &
  c%h_c_rates%cooling_FeII_rate                          , &
  c%h_c_rates%cooling_OH_rot_rate                        , &
  c%h_c_rates%cooling_Neufeld_H2O_rate_rot               , &
  c%h_c_rates%cooling_Neufeld_H2O_rate_vib               , &
  c%h_c_rates%cooling_Neufeld_CO_rate_rot                , &
  c%h_c_rates%cooling_Neufeld_CO_rate_vib                , &
  c%h_c_rates%cooling_Neufeld_H2_rot_rate                , &
  c%h_c_rates%cooling_LymanAlpha_rate                    , &
  c%h_c_rates%cooling_free_bound_rate                    , &
  c%h_c_rates%cooling_free_free_rate                     , &
  c%par%alpha_viscosity                                  , &
  c%par%ambipolar_f                                      , &
  c%par%ion_charge                                       , &
  c%par%velo_Kepler                                      , &
  c%par%omega_Kepler                                     , &
  c%par%velo_gradient                                    , &
  c%par%sound_speed                                      , &
  c%par%velo_width_turb                                  , &
  c%par%coherent_length                                  , &
  c%par%SitesPerGrain                                    , &
  c%par%n_mol_on_grain                                   , &
  c%abundances
end subroutine disk_save_results_write


subroutine disk_calc_disk_mass
  integer i
  double precision m
  m = 0D0
  do i=1, leaves%nlen
    associate(p => leaves%list(i)%p%par)
      m = m + p%mgas_cell
        !p%n_gas * p%MeanMolWeight * (phy_2Pi * p%rcen * p%dr * p%dz)
    end associate
  end do
  a_disk%disk_mass_in_Msun = m * 2D0 / phy_Msun_CGS ! Accounts for the z<0 side
#ifdef DIAGNOSIS_TRACK_FUNC_CALL
  write(*,*) 'Running disk_calc_disk_mass'
  write(*,*) 'a_disk%disk_mass_in_Msun = ', a_disk%disk_mass_in_Msun
#endif
end subroutine disk_calc_disk_mass


subroutine disk_calc_dust_components_mass(mdisk_dusts_in_Msun)
  double precision, dimension(:), intent(out) :: mdisk_dusts_in_Msun
  integer i, j
  double precision vol
  type(type_cell), pointer :: c
  !
  mdisk_dusts_in_Msun = 0D0
  !
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    vol = phy_Pi * (c%xmax + c%xmin) * &
          (c%xmax - c%xmin) * (c%ymax-c%ymin) * phy_AU2cm**3
    do j=1, a_disk%ndustcompo
      mdisk_dusts_in_Msun(j) = mdisk_dusts_in_Msun(j) + vol * c%par%rho_dusts(j)
    end do
  end do
  mdisk_dusts_in_Msun = mdisk_dusts_in_Msun * 2D0 &
                        / phy_Msun_CGS
end subroutine disk_calc_dust_components_mass


subroutine set_hc_chem_params_from_cell(id)
  integer id
  !
  hc_params = leaves%list(id)%p%par
  !
  hc_params%X_H2    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CI)
  hc_params%X_CII   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CII)
  hc_params%X_OI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OI)
  hc_params%X_NII   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_NII)
  hc_params%X_FeII  = leaves%list(id)%p%abundances(chem_idx_some_spe%i_FeII)
  hc_params%X_SiII  = leaves%list(id)%p%abundances(chem_idx_some_spe%i_SiII)
  hc_params%X_CO    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OH)
  hc_params%X_E     = leaves%list(id)%p%abundances(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = leaves%list(id)%p%abundances(chem_idx_some_spe%i_Hplus)
  hc_params%X_Heplus= leaves%list(id)%p%abundances(chem_idx_some_spe%i_Heplus)
  if (chem_idx_some_spe%i_gH .gt. 0) then
    hc_params%X_gH   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_gH)
  end if
  !
  call chem_cal_rates
  !
  hc_params%R_H2_form_rate = &
    get_H2_form_rate( &
      hc_params%R_H2_form_rate_coeff, &
      hc_params%X_gH, &
      hc_params%X_HI, &
      hc_params%n_gas)
  !
  hc_Tgas  = leaves%list(id)%p%par%Tgas
  hc_Tdust = leaves%list(id)%p%par%Tdust
  !
end subroutine set_hc_chem_params_from_cell


subroutine disk_set_a_cell_params(c, cell_params_copy, asCopied, only_rescal)
  integer i
  type(type_cell), intent(inout) :: c
  type(type_cell_rz_phy_basic), intent(in), optional :: cell_params_copy
  logical, intent(in), optional :: asCopied, only_rescal
  logical asCop, onlyrescal
  integer stat
  double precision, parameter :: n_gas_min_val = 1D-30
  if (present(asCopied)) then
    asCop = asCopied
  else
    asCop = .false.
  end if
  if (present(only_rescal)) then
    onlyrescal = only_rescal
  else
    onlyrescal = .false.
  end if
  !
  if (.not. associated(c%par)) then
    allocate(c%par)
  end if
  if (.not. allocated(c%h_c_rates)) then
    allocate(c%h_c_rates)
  end if
  !
  if (.not. allocated(c%abundances)) then
    allocate(c%abundances(chem_species%nSpecies), &
             c%col_den_toISM(chem_idx_some_spe%nItem), &
             c%col_den_toStar(chem_idx_some_spe%nItem), stat=stat)
    if (stat .ne. 0) then
      write(*, '(A)') 'Error in allocating array!'
      write(*, '(A)') 'In disk_set_a_cell_params.'
      write(*, '(A, I16/)') 'STAT = ', stat
      call error_stop()
    end if
  end if
  !
  c%iIter = 0
  c%quality = 0
  !
  if (present(cell_params_copy)) then
    c%par = cell_params_copy
  end if
  !
  !c%abundances(1:chem_species%nSpecies) = chemsol_stor%y0(1:chem_species%nSpecies)
  c%abundances(1:chem_species%nSpecies) = 0D0
  c%col_den_toISM = 0D0
  c%col_den_toStar = 0D0
  !
  c%par%volume = phy_Pi * (c%xmax + c%xmin) * &
                 (c%xmax - c%xmin) * (c%ymax-c%ymin) * phy_AU2cm**3
  c%par%area_T = phy_Pi * (c%xmax + c%xmin) * (c%xmax-c%xmin) * phy_AU2cm**2
  c%par%area_B = phy_Pi * (c%xmax + c%xmin) * (c%xmax-c%xmin) * phy_AU2cm**2
  c%par%area_I = phy_2Pi * c%xmin * (c%ymax-c%ymin) * phy_AU2cm**2
  c%par%area_O = phy_2Pi * c%xmax * (c%ymax-c%ymin) * phy_AU2cm**2
  c%par%surf_area = c%par%area_T + c%par%area_B + c%par%area_I + c%par%area_O
  !
  a_disk%andrews_gas%particlemass = c%par%MeanMolWeight * phy_mProton_CGS
  !
  c%par%ndustcompo = a_disk%ndustcompo
  !
  if (.not. asCop) then
    if (.not. onlyrescal) then
      c%par%rho_dusts = 0D0
      c%par%mp_dusts = 0D0
      c%par%n_dusts = 0D0
      c%par%ndust_tot = 0D0
      c%par%sig_dusts = 0D0
      c%par%sigdust_ave = 0D0
      c%par%mdusts_cell = 0D0
      c%par%mdust_tot = 0D0
      ! When doing refinement, no need to recalculate the densities.
      do i=1, a_disk%ndustcompo
        ! Dust mass density
        c%par%rho_dusts(i) = get_ave_val_analytic( &
                c%xmin, c%xmax, c%ymin, c%ymax, &
                a_disk%dustcompo(i)%andrews)
        c%par%mp_dusts(i) = a_disk%dustcompo(i)%pmass_CGS ! Dust particle mass in gram
        !
        c%par%sig_dusts(i) = phy_Pi * a_disk%dustcompo(i)%mrn%r2av * &
                             phy_micron2cm**2
        !write(str_disp, '(A, 2ES16.6, A, I2, A, ES12.3, A, ES12.3)') &
        !    'xymin:', c%xmin, c%ymin, ' Dust:', i, ' n_d:', c%par%n_dusts(i), &
        !    ' sig:', c%par%sig_dusts(i)
        !call display_string_both(str_disp, a_book_keeping%fU, onlyfile=.true.)
      end do
    end if
    !
    if (a_disk_iter_params%rescale_ngas_2_rhodust) then
      if (a_disk_iter_params%dust2gas_mass_ratio_deflt .gt. 0D0) then
        c%par%n_gas = sum(c%par%rho_dusts(1:a_disk%ndustcompo)) / &
          a_disk_iter_params%dust2gas_mass_ratio_deflt / &
          (phy_mProton_CGS*c%par%MeanMolWeight)
      else
        c%par%n_gas = 0D0
        do i=1, a_disk%ndustcompo
          c%par%n_gas = c%par%n_gas + &
            c%par%rho_dusts(i) / &
            a_disk%dustcompo(i)%d2g_ratio / &
            (phy_mProton_CGS*c%par%MeanMolWeight)
        end do
      end if
    else
      if (.not. onlyrescal) then
        c%par%n_gas = get_ave_val_analytic(c%xmin, c%xmax, c%ymin, c%ymax, &
                                         a_disk%andrews_gas)
      end if
    end if
    !
    if (c%par%n_gas .le. 0D0) then
      c%par%n_gas = n_gas_min_val
    end if
    !
  end if
  !
  c%par%sig_dusts0 = c%par%sig_dusts
  c%par%en_exchange = 0D0
  c%par%en_exchange_per_vol = 0D0
  c%par%en_exchange_tot = 0D0
  !
  call calc_dustgas_struct_snippet1(c)
  !
  call calc_dustgas_struct_snippet2(c)
  !
  !write(str_disp, '(A, ES12.3, A, ES12.3, A, ES12.3)') &
  !  'nd_tot:', c%par%ndust_tot, ' sig_d_ave:', c%par%sigdust_ave, &
  !  ' r_d:', c%par%GrainRadius_CGS
  !call display_string_both(str_disp, a_book_keeping%fU, onlyfile=.true.)
  !
  if (grid_config%use_data_file_input) then
    c%par%Tgas    = c%val(2)
    c%par%Tdust   = c%val(2)
  else
    !c%par%Tgas    = 0D0
    !c%par%Tdust   = 0D0
    ! 2014-06-11 Wed 03:28:07
    ! Another deeply hidden trivial-looking bug!
    ! The initial Tgas cannot be zero, because it will be used to calculate the
    ! HI scattering cross section.
    c%par%Tgas    = 600D0 / (1D0 + (c%xmin+c%xmax)*0.5D0) * &
                    (1D0 + (c%ymin+c%ymax)*0.5D0)
    c%par%Tdust   = 0D0 ! instead of c%par%Tgas
  end if
  !
  c%par%pressure_thermal = 0D0
  c%par%gravity_z = 0D0
  c%par%gravity_acc_z = 0D0
  !
  ! Calculate the local velocity gradient, thermal velocity width, turbulent
  ! width, coherent length
  !
  call calc_local_dynamics(c, init=.true.)
  !
end subroutine disk_set_a_cell_params


subroutine calc_local_dynamics(c, init)
  type(type_cell), intent(inout) :: c
  logical, intent(in), optional :: init
  double precision r, total_gas_abundance
  logical is_init
  integer i
  !
  if (present(init)) then
    is_init = init
  else
    is_init = .false.
  end if
  !
  r = (c%xmin + c%xmax) * 0.5D0
  c%par%velo_Kepler = &
    sqrt( &
      (phy_GravitationConst_CGS * a_star%mass * &
      phy_Msun_CGS) / (r * phy_AU2cm))
  c%par%omega_Kepler = c%par%velo_Kepler / (r * phy_AU2cm)
  c%par%velo_gradient = 0.5D0 * c%par%velo_Kepler / &
                        (r * phy_AU2cm)
  c%par%Neufeld_dv_dz = c%par%velo_gradient * 1D-5 ! cm s-1 to km s-1
  c%par%Neufeld_G     = 1D0
  !
  if (c%par%Tgas .le. 0D0) then
    if (is_init) then
      ! Use a guess value for initialization
      c%par%sound_speed = sqrt(phy_kBoltzmann_CGS * 100D0 / &
            (phy_mProton_CGS * c%par%MeanMolWeight*2D0))
    else
      c%par%sound_speed = 0D0
    end if
  else
    c%par%sound_speed = sqrt(phy_kBoltzmann_CGS*c%par%Tgas / &
            (phy_mProton_CGS * c%par%MeanMolWeight*2D0))
  end if
  !
  c%par%velo_width_turb = c%par%sound_speed
  !
  c%par%coherent_length = c%par%velo_width_turb / c%par%velo_gradient
  !
  call calc_gravity_single_cell(c)
  !
  c%par%ion_charge = get_ion_charge(c)
  c%par%ambipolar_f = c%par%n_gas * c%par%ion_charge *  &
    beta_ion_neutral_colli / c%par%omega_Kepler
  !
  if (.not. a_disk%use_fixed_alpha_visc) then
    if (is_init) then
      c%par%alpha_viscosity = a_disk%base_alpha
    else
      c%par%alpha_viscosity = get_alpha_viscosity(c%par%ambipolar_f) * &
        a_disk%base_alpha
      if (isnan(c%par%alpha_viscosity)) then
        write(str_disp, '(4ES16.6)') c%par%alpha_viscosity, c%par%ambipolar_f, &
          c%par%ion_charge, c%par%omega_Kepler
        call display_string_both(str_disp, a_book_keeping%fU, onlyfile=.true.)
      end if
    end if
  end if
  !
  total_gas_abundance = 0D0
  do i=1, chem_species%nSpecies
    if ((chem_species%names(i)(1:1) .eq. 'g') .or. &
        (chem_species%names(i)(1:4) .eq. 'Grai')) then
      cycle
    end if
    total_gas_abundance = total_gas_abundance + c%abundances(i)
  end do
  c%par%pressure_thermal = &
    c%par%n_gas * c%par%Tgas * phy_kBoltzmann_CGS * &
    total_gas_abundance
end subroutine calc_local_dynamics


subroutine calc_gravity_single_cell(c)
  type(type_cell), intent(inout) :: c
  double precision R3
  R3 = (sqrt(((c%xmax+c%xmin)*0.5D0)**2 + &
             ((c%ymax+c%ymin)*0.5D0)**2))**3
  c%par%gravity_z = &
      phy_GravitationConst_CGS * a_star%mass * phy_Msun_CGS * &
      (c%par%mgas_cell + &
       c%par%mdust_tot) * &
      (-((c%ymax+c%ymin)*0.5D0) / R3 / (phy_AU2cm**2))
end subroutine calc_gravity_single_cell


subroutine calc_gravity_column(c)
  type(type_cell), intent(inout) :: c
  if (c%above%n .eq. 0) then
    c%par%gravity_acc_z = &
      phy_GravitationConst_CGS * (a_star%mass * phy_Msun_CGS) * &
        (calc_Ncol_from_cell_to_point(c, (c%xmax+c%xmin)*0.5D0, &
            root%ymax*2D0, -4) * &
         c%par%area_T * phy_mProton_CGS * c%par%MeanMolWeight) * &
        (-c%ymax / (sqrt(c%xmax**2 + c%ymax**2))**3 / (phy_AU2cm**2))
  else
    c%par%gravity_acc_z = &
      calc_Ncol_from_cell_to_point(c, (c%xmax+c%xmin)*0.5D0, &
            root%ymax*2D0, -1)
  end if
end subroutine calc_gravity_column


pure subroutine get_alpha_viscosity_alt(cpar, y, n, base_val)
  type(type_cell_rz_phy_basic), intent(inout) :: cpar
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: y
  double precision, intent(in) :: base_val
  !
  cpar%ion_charge = get_ion_charge_y(y, n)
  !
  cpar%ambipolar_f = cpar%n_gas * cpar%ion_charge * &
    beta_ion_neutral_colli / cpar%omega_Kepler
  !
  cpar%alpha_viscosity = base_val * get_alpha_viscosity(cpar%ambipolar_f)
end subroutine get_alpha_viscosity_alt


pure function get_ion_charge(c)
  double precision get_ion_charge
  type(type_cell), intent(in) :: c
  integer i
  integer, parameter :: iCharge = 1
  get_ion_charge = 0D0
  do i=1, chem_species%nSpecies
    if (chem_species%elements(iCharge, i) .gt. 0) then
      get_ion_charge = get_ion_charge + &
        dble(chem_species%elements(iCharge, i)) * c%abundances(i)
    end if
  end do
end function get_ion_charge


pure function get_ion_charge_y(y, n)
  double precision get_ion_charge_y
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: y
  integer i
  integer, parameter :: iCharge = 1
  get_ion_charge_y = 0D0
  do i=1, chem_species%nSpecies
    if ((y(i) .ge. 1D-30) .and. (chem_species%elements(iCharge, i) .gt. 0)) then
      get_ion_charge_y = get_ion_charge_y + &
        dble(chem_species%elements(iCharge, i)) * y(i)
    end if
  end do
end function get_ion_charge_y


pure function get_alpha_viscosity(am)
  double precision get_alpha_viscosity
  double precision, intent(in) :: am
  double precision, parameter :: smallnum = 1D-20
  double precision tmp, tmp1, tmp2
  if (am .le. smallnum) then
    get_alpha_viscosity = 0D0
  else
    tmp = log(am)
    tmp1 = exp(-2.4D0 * tmp)
    tmp2 = exp(-0.3D0 * tmp)
    get_alpha_viscosity = &
      0.5D0 / sqrt(2500D0*tmp1 + (8D0*tmp2 + 1D0)**2)
  end if
end function get_alpha_viscosity


subroutine disk_set_gridcell_params_runonce
  integer i
  double precision, dimension(MaxNumOfDustComponents) :: mdisk_dusts
  double precision, dimension(MaxNumOfDustComponents) :: f_rescal_dust
  double precision mdisk_gas, f_rescal_gas
  do i=1, leaves%nlen
    call disk_set_a_cell_params(leaves%list(i)%p, cell_params_ini, &
            asCopied=.false., only_rescal=.false.)
  end do
  !
  call disk_calc_dust_components_mass(mdisk_dusts)
  mdisk_gas = calc_disk_gas_mass()
  write(*,*) 'mdisk_dusts: ', mdisk_dusts
  write(*,*) 'mdisk_gas: ', mdisk_gas
  !
  do i=1, a_disk%ndustcompo
    f_rescal_dust(i) = a_disk%dustcompo(i)%andrews%Md / mdisk_dusts(i)
  end do
  f_rescal_gas = a_disk%andrews_gas%Md / mdisk_gas
  !
  do i=1, leaves%nlen
    call rescale_dust_density(leaves%list(i)%p, f_rescal_dust)
    if (.not. a_disk_iter_params%rescale_ngas_2_rhodust) then
      leaves%list(i)%p%par%n_gas = leaves%list(i)%p%par%n_gas * f_rescal_gas
    end if
    call disk_set_a_cell_params(leaves%list(i)%p, &
            asCopied=.false., only_rescal=.true.)
  end do
end subroutine disk_set_gridcell_params_runonce



subroutine rescale_dust_density(c, fresc)
  type(type_cell), intent(inout) :: c
  double precision, dimension(:), intent(in) :: fresc
  integer i
  do i=1, a_disk%ndustcompo
    c%par%rho_dusts(i) = c%par%rho_dusts(i) * fresc(i)
  end do
end subroutine rescale_dust_density




subroutine calc_dust_MRN_par(mrn)
  type(type_dust_MRN), intent(inout) :: mrn
  double precision tmp1, tmp2, norm
  double precision, parameter :: smallnum = 1D-6
  !
  tmp1 = exp((1D0 - mrn%n) * log(mrn%rmin)) ! = rmin**(1-n)
  tmp2 = exp((1D0 - mrn%n) * log(mrn%rmax)) ! = rmax**(1-n)
  if (abs(mrn%n - 1D0) .le. smallnum) then
    norm = log(mrn%rmax/mrn%rmin)
  else
    norm = (tmp2 - tmp1) / (1D0 - mrn%n)
  end if
  if (abs(mrn%n - 2D0) .le. smallnum) then
    mrn%rav = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%rav  = (tmp2 * mrn%rmax    - tmp1 * mrn%rmin) &
        / ((2D0 - mrn%n) * norm)
  end if
  if (abs(mrn%n - 3D0) .le. smallnum) then
    mrn%r2av = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%r2av = (tmp2 * mrn%rmax**2 - tmp1 * mrn%rmin**2) &
        / ((3D0 - mrn%n) * norm)
  end if
  if (abs(mrn%n - 4D0) .le. smallnum) then
    mrn%r3av = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%r3av = (tmp2 * mrn%rmax**3 - tmp1 * mrn%rmin**3) &
        / ((4D0 - mrn%n) * norm)
  end if
end subroutine calc_dust_MRN_par


subroutine save_chem_rates(i0)
  integer, intent(in) :: i0
  integer fU, k
  character(len=128) filename, dir
  character(len=9) strtmp
  ! Use namelist for output some logging infomation.
  ! Not very readable, but easy to implement.
  namelist /cell_par_log/ hc_params
  !
  write(filename, '("reac_rates_cell_", I4.4, ".dat")') i0
  dir = trim(combine_dir_filename(a_book_keeping%dir, 'rates_log/'))
  if (.NOT. dir_exist(dir)) then
    call my_mkdir(dir)
  end if
  call openFileSequentialWrite(fU, combine_dir_filename(dir, filename), &
       99999, getu=1)
  !
  if (associated(hc_params)) then
    write(fU, nml=cell_par_log)
  end if
  !
  do k=1, chem_net%nReactions
    !write(fU, '(A135, ES16.4E4)') chem_reac_str%list(k), chem_net%rates(k)
    call double2str(strtmp, chem_net%ABC(3, k), 9, 1)
    write(fU, &
      '(7(A12), ES9.2, F9.2, A9, 2I6, I3, X, A1, X, A2, ES16.6E3)') &
      chem_net%reac_names(:,k), &
      chem_net%prod_names(:,k), &
      chem_net%ABC(1:2,k), &
      strtmp, &
      int(chem_net%T_range(:,k)), &
      chem_net%itype(k), &
      chem_net%reliability(k), &
      chem_net%ctype(k), &
      chem_net%rates(k)
  end do
  close(fU)
end subroutine save_chem_rates


subroutine save_post_config_params
  type(type_disk_basic_params) disk_params_tmp
  namelist /disk_params_log/ disk_params_tmp
  disk_params_tmp = a_disk
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=disk_params_log)
    flush(a_book_keeping%fU)
  end if
end subroutine save_post_config_params


subroutine load_refine_check_species
  integer fU, i, i1, ios, n
  character(len=const_len_init_abun_file_row) str
  n = GetFileLen_comment_blank( &
      combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
        a_disk_iter_params%filename_list_check_refine), '!')
  allocate(idx_Species_check_refine(n), &
           thr_Species_check_refine(n))
  call openFileSequentialRead(fU, &
    combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
      a_disk_iter_params%filename_list_check_refine), 99, getu=1)
  i1 = 0
  do
    read(fU, FMT='(A)', IOSTAT=ios) str
    if (ios .NE. 0) then
      exit
    end if
    do i=1, chem_species%nSpecies
      if (trim(str(1:const_len_species_name)) .EQ. chem_species%names(i)) then
        i1 = i1 + 1
        idx_Species_check_refine(i1) = i
        read(str(const_len_species_name+1:const_len_init_abun_file_row), &
          '(F7.0)') thr_Species_check_refine(i1)
        exit
      end if
    end do
  end do
  close(fU)
  a_disk_iter_params%nSpecies_check_refine = i1
  write(str_disp, '("! Species used for checking refinement:")')
  call display_string_both(str_disp, a_book_keeping%fU)
  do i=1, a_disk_iter_params%nSpecies_check_refine
    write(str_disp, '("! ", A12, ES12.2)') &
      chem_species%names(idx_Species_check_refine(i)), &
      thr_Species_check_refine(i)
    call display_string_both(str_disp, a_book_keeping%fU)
  end do
end subroutine load_refine_check_species


subroutine do_refine
  integer i, n_refine
  a_disk_iter_params%ncell_refine = 0
  do i=1, leaves%nlen
    if (need_to_refine(leaves%list(i)%p, n_refine)) then
      a_disk_iter_params%ncell_refine = a_disk_iter_params%ncell_refine + 1
      write(str_disp, '("!", I4, A, 4ES12.2, " into ", I4, " parts.")') &
        a_disk_iter_params%ncell_refine, ' Refining ', &
        leaves%list(i)%p%xmin, leaves%list(i)%p%xmax, &
        leaves%list(i)%p%ymin, leaves%list(i)%p%ymax, n_refine
      call display_string_both(str_disp, a_book_keeping%fU, onlyfile=.true.)
      !
      call refine_this_cell_vertical(leaves%list(i)%p, n_refine)
      !
    end if
  end do
  if (a_disk_iter_params%ncell_refine .ge. 1) then
    call remake_index
    call vertical_pressure_gravity_balance_alt(a_disk%star_mass_in_Msun, &
      fix_dust_struct=a_disk_iter_params%vertical_structure_fix_dust, &
      disk_gas_mass_preset=a_disk%disk_mass_in_Msun)
    call remake_index
  end if
end subroutine do_refine


subroutine refine_after_vertical
  type(type_cell), pointer :: c
  integer i, n_div
  integer nleaves_now
  !
  if (leaves%nlen .gt. a_disk_iter_params%max_num_of_cells) then
    return
  end if
  !
  nleaves_now = leaves%nlen
  !
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    call get_ndiv(c, n_div)
    if (n_div .le. 1) then
      cycle
    end if
    !
    nleaves_now = nleaves_now + n_div - 1
    !
    call refine_this_cell_vertical(c, n_div)
    !
    if (nleaves_now .gt. a_disk_iter_params%max_num_of_cells) then
      exit
    end if
  end do
  ! Remake the index because new cells are created during the refinement
  call remake_index
  !
  ! Readjust the densities after refining
  call vertical_pressure_gravity_balance_alt(a_disk%star_mass_in_Msun, &
    useTdust=.true., Tdust_lowerlimit=a_disk_iter_params%minimum_Tdust, &
    ndust_lowerlimit=grid_config%min_val_considered*1D-18, &
    ngas_lowerlimit=grid_config%min_val_considered, &
    disk_gas_mass_preset = a_disk%disk_mass_in_Msun, &
    fix_dust_struct=a_disk_iter_params%vertical_structure_fix_dust)
  !
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    if (.not. c%using) then
      call deallocate_when_not_using(c)
    end if
  end do
  !
  call remake_index
end subroutine refine_after_vertical



subroutine merge_cells
  type(type_cell), pointer :: prt
  integer i, j, itmp, stat
  !
  do i=1, leaves%nlen
    if (.not. associated(leaves%list(i)%p)) then
      cycle
    end if
    !
    if (.not. leaves%list(i)%p%using) then
      call deallocate_when_not_using(leaves%list(i)%p)
      nullify(leaves%list(i)%p)
      cycle
    end if
    !
    prt => leaves%list(i)%p%parent
    !
    if (prt%using) then
      write(*,*) 'In merge_cells: prt already using!'
      call error_stop()
    end if
    !
    if (size(prt%children) .ne. prt%nChildren) then
      write(*,*) 'In merge_cells: wrong children number!'
      write(*,*) prt%nChildren, prt%order, prt%nOffspring, prt%nleaves, prt%using
      write(*,*) prt%xmin, prt%xmax, prt%ymin, prt%ymax
      call error_stop()
    end if
    !
    if (need_to_merge(prt)) then
      !
      call set_par_from_children(prt)
      !
      do j=1, prt%nChildren
        prt%children(j)%p%using = .false.
        prt%children(j)%p%id = -1
        !
        itmp = prt%children(j)%p%id
        if (prt%children(j)%p%using) then
          nullify(leaves%list(itmp)%p)
        end if
        call deallocate_when_not_using(prt%children(j)%p)
        !
        nullify(prt%children(j)%p%parent)
        !deallocate(prt%children(j)%p)
      end do
      !
      deallocate(prt%children, stat=stat)
      if (stat .ne. 0) then
        write(*,*) 'In merge_cells2: Fail to deallocate children!'
        call error_stop()
      end if
      prt%using = .true.
      prt%converged = .false.
      prt%nOffspring = 0
      prt%nChildren = 0
      prt%nleaves = 1
    end if
  end do
  !
  call remake_index
  !
end subroutine merge_cells


function need_to_merge(prt) result(shouldMerge)
  type(type_cell), intent(in) :: prt
  type(type_cell), pointer :: c
  logical shouldMerge
  integer i
  double precision min_dy_here, max_dy_here
  double precision, dimension(6) :: maxs, mins
  double precision, parameter :: maxdz_ratio = 0.02D0
  !
  if (prt%nChildren .le. 1) then
    shouldMerge = .true.
    return
  end if
  !
  min_dy_here = grid_config%small_len_frac * &
    0.5D0 * sqrt((prt%xmin+prt%xmax)**2 + (prt%ymin+prt%ymax)**2)
  if ((prt%ymax - prt%ymin) .lt. min_dy_here) then
    shouldMerge = .true.
    return
  end if
  !
  max_dy_here = (prt%xmax+prt%xmin) * 0.5D0 * maxdz_ratio
  if ((prt%ymax - prt%ymin) .gt. max_dy_here) then
    shouldMerge = .false.
    return
  end if
  !
  maxs = 0D0
  mins = 1D199
  do i=1, prt%nChildren
    if (.not. associated(prt%children(i)%p)) then
      cycle
    end if
    c => prt%children(i)%p
    if (.not. associated(c%par)) then
      cycle
    end if
    maxs(1) = max(maxs(1), c%par%n_gas)
    maxs(2) = max(maxs(2), c%par%Tdust)
    maxs(3) = max(maxs(3), c%par%Av_toStar)
    maxs(4) = max(maxs(4), c%par%Av_toISM)
    maxs(5) = max(maxs(5), c%par%flux_Xray)
    maxs(6) = max(maxs(6), c%par%flux_UV)
    mins(1) = min(mins(1), c%par%n_gas)
    mins(2) = min(mins(2), c%par%Tdust)
    mins(3) = min(mins(3), c%par%Av_toStar)
    mins(4) = min(mins(4), c%par%Av_toISM)
    mins(5) = min(mins(5), c%par%flux_Xray)
    mins(6) = min(mins(6), c%par%flux_UV)
  end do
  shouldMerge = &
    (maxs(1) / mins(1) .le. grid_config%max_ratio_to_be_uniform) .and. & ! 1
    (maxs(2) / mins(2) .le. 1.1D0) .and. & ! 2
    ((maxs(3) / (mins(3)+1D-20) .le. 1.2D0) .or. & ! 3
     ((maxs(3) / (mins(3)+1D-20) .le. 3D0) .and. (maxs(3) .le. 1D-5))) .and. &
    ((maxs(4) / (mins(4)+1D-20) .le. 1.2D0) .or. & ! 4
     ((maxs(4) / (mins(4)+1D-20) .le. 3D0) .and. (maxs(4) .le. 1D-5))) .and. &
    (maxs(5) / (mins(5)+1D-20) .le. 1.2D0) .and. & ! 5
    (maxs(6) / (mins(6)+1D-20) .le. 1.2D0) ! 6
end function need_to_merge


subroutine set_par_from_children(prt)
  type(type_cell), pointer, intent(inout) :: prt
  integer i, nsum
  call set_cell_par_preliminary(prt)
  do i=1, prt%nChildren
    if (associated(prt%children(i)%p%par)) then
      call disk_set_a_cell_params(prt, prt%children(i)%p%par, &
            asCopied=.true., only_rescal=.false.)
      exit
    end if
  end do
  prt%par%Tgas       = 0D0
  prt%par%Tdusts     = 0D0
  prt%par%Tdust      = 0D0
  prt%abundances     = 0D0
  prt%col_den_toISM  = 0D0
  prt%col_den_toStar = 0D0
  nsum = 0
  do i=1, prt%nChildren
    if (associated(prt%children(i)%p%par)) then
      nsum = nsum + 1
      prt%par%Tgas       = prt%par%Tgas      + prt%children(i)%p%par%Tgas
      prt%par%Tdusts     = prt%par%Tdusts    + prt%children(i)%p%par%Tdusts
      prt%par%Tdust      = prt%par%Tdust     + prt%children(i)%p%par%Tdust
      prt%abundances     = prt%abundances    + prt%children(i)%p%abundances
      prt%col_den_toISM  = prt%col_den_toISM + prt%children(i)%p%col_den_toISM
      prt%col_den_toStar = prt%col_den_toStar+ prt%children(i)%p%col_den_toStar
    end if
  end do
  prt%par%Tgas       =  prt%par%Tgas       / dble(nsum)
  prt%par%Tdusts     =  prt%par%Tdusts     / dble(nsum)
  prt%par%Tdust      =  prt%par%Tdust      / dble(nsum)
  prt%abundances     =  prt%abundances     / dble(nsum)
  prt%col_den_toISM  =  prt%col_den_toISM  / dble(nsum)
  prt%col_den_toStar =  prt%col_den_toStar / dble(nsum)
end subroutine set_par_from_children



subroutine remake_index
  !
  write(*, '(A)') 'Remaking the global index...'
  call get_number_of_leaves(root)
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
  call deallocate_columns
  call make_columns
  write(*, '(A, I8)') 'New number of bottom cells:', bott_cells%nlen
  write(*, '(A, I8)') 'New number of leaf cells:', root%nleaves
  !
  call load_ana_points_list ! Reload, actually
  !
  call allocate_iter_stor
  !
  call disk_calc_disk_mass
  write(*, '(A, ES12.4)') 'Disk gas mass (in Msun) = ', a_disk%disk_mass_in_Msun
end subroutine remake_index


function need_to_refine(c, n_refine)
  logical need_to_refine
  type(type_cell), target :: c
  integer, intent(out), optional :: n_refine
  integer i, i0, i1, j
  double precision val_max, val_min
  logical flag1, flag2
  integer, parameter :: n_refine_max = 10
  flag1 = .false.
  flag2 = .false.
  if (present(n_refine)) then
    n_refine = 0
  end if
  if ((c%ymax-c%ymin) .le. grid_config%smallest_cell_size) then
    need_to_refine = .false.
    return
  end if
  do i=1, c%above%n
    i0 = c%above%idx(i)
    !if (leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      if (val_max .gt. thr_Species_check_refine(j)) then
        if (val_max / val_min .gt. a_disk_iter_params%threshold_ratio_refine) then
          flag1 = .true.
          if (present(n_refine)) then
            n_refine = max(n_refine, int(log10(val_max / val_min)) * 2)
          end if
        end if
      end if
    end do
  end do
  do i=1, c%below%n
    i0 = c%below%idx(i)
    !if (leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      if (val_max .gt. thr_Species_check_refine(j)) then
        if (val_max / val_min .gt. a_disk_iter_params%threshold_ratio_refine) then
          flag2 = .true.
          if (present(n_refine)) then
            n_refine = max(n_refine, int(log10(val_max / val_min)) * 2)
          end if
        end if
      end if
    end do
  end do
  need_to_refine = flag1 .or. flag2
  n_refine = min(n_refine, n_refine_max)
  return
end function need_to_refine


subroutine refine_this_cell_vertical(c, n)
  ! c is a working cell that needs to be refined.
  type(type_cell), target :: c
  double precision dy
  integer, intent(in), optional :: n
  integer i, ndivide
  !
  if (present(n)) then
    ndivide = n
  else
    ndivide = 3
  end if
  !
  if (ndivide .lt. 2) then
    return
  end if
  !
  c%nleaves = ndivide
  c%nChildren = ndivide
  call init_children(c, ndivide)
  !
  dy = (c%ymax - c%ymin) / dble(ndivide)
  !
  do i=1, c%nChildren
    associate(cc => c%children(i)%p)
      cc%xmin = c%xmin
      cc%xmax = c%xmax
      cc%ymin = c%ymin + dble(i-1) * dy
      cc%ymax = c%ymin + dble(i)   * dy
      !
      ! Re-interpolate density from the input data.
      call set_cell_par_preliminary(cc)
      cc%using = .true.
      cc%converged = .false.
      cc%nOffspring = 0
      cc%nChildren = 0
      cc%nleaves = 1
      !
      cc%iIter = c%iIter
      !
      call disk_set_a_cell_params(cc, c%par, asCopied=.true.)
      cc%par%Tgas = c%par%Tgas
      cc%par%Tdusts = c%par%Tdusts
      cc%par%Tdust  = c%par%Tdust
      !
      cc%h_c_rates = c%h_c_rates
      cc%abundances = c%abundances
      cc%col_den_toISM = c%col_den_toISM
      cc%col_den_toStar = c%col_den_toStar
    end associate
  end do
  ! Avoid numerical roundings
  c%children(1)%p%ymin       = c%ymin
  c%children(ndivide)%p%ymax = c%ymax
  do i=1, c%nChildren-1
    c%children(i)%p%ymax = c%children(i+1)%p%ymin
  end do
  !
  ! Deactivate c
  c%using = .false.
  c%converged = .false.
  c%id = -1
end subroutine refine_this_cell_vertical


subroutine load_ana_species_list
  integer fU, ios, i, n
  integer, dimension(:), allocatable :: list_tmp
  character(len=12) str
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_species), 99, getu=1)
  allocate(list_tmp(chem_species%nSpecies))
  n = 0
  do
    read(fU, '(A12)', iostat=ios) str
    if (ios .ne. 0) then
      exit
    end if
    do i=1, chem_species%nSpecies
      if (chem_species%names(i) .eq. str) then
        !
        if (.not. is_in_list_int(i, n, list_tmp(1:n))) then
          n = n + 1
          list_tmp(n) = i
        end if
        !
        exit
        !
      end if
    end do
  end do
  close(fU)
  !
  ana_splist%nlen = n
  if (n .gt. 0) then
    if (allocated(ana_splist%vals)) then
      deallocate(ana_splist%vals)
    end if
    allocate(ana_splist%vals(n))
    ana_splist%vals = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine load_ana_species_list



subroutine load_ana_points_list
  integer fU, ios, i, n
  double precision r, z
  integer, dimension(:), allocatable :: list_tmp
#ifdef DIAGNOSIS_TRACK_FUNC_CALL
  write(*,*) 'Running load_ana_points_list'
#endif
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_points), 99, getu=1)
  allocate(list_tmp(leaves%nlen))
  n = 0
  do
    read(fU, '(2F6.0)', iostat=ios) r, z
    if (ios .gt. 0) then
      ! Input error
      cycle
    else if (ios .lt. 0) then
      ! End of file
      exit
    end if
    do i=1, leaves%nlen
      if ((leaves%list(i)%p%xmin .le. r) .and. (leaves%list(i)%p%xmax .ge. r) .and. &
          (leaves%list(i)%p%ymin .le. z) .and. (leaves%list(i)%p%ymax .ge. z)) then
        !
        if (.not. is_in_list_int(i, n, list_tmp(1:n))) then
          n = n + 1
          list_tmp(n) = i
        end if
        !
        exit
        !
      end if
    end do
  end do
  close(fU)
  !
  ana_ptlist%nlen = n
  if (n .gt. 0) then
    if (allocated(ana_ptlist%vals)) then
      deallocate(ana_ptlist%vals)
    end if
    allocate(ana_ptlist%vals(n))
    ana_ptlist%vals = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine load_ana_points_list



subroutine chem_analyse(id)
  use quick_sort
  integer, intent(in) :: id
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=128) fname_pre
  character(len=32) FMTstryHistory
  double precision dy_y, dt_t
  double precision frac
  double precision, dimension(:,:), allocatable :: reacHeats
  double precision, dimension(:), allocatable :: tmp
  frac = 0.1D0
  !
  write(*, '(A/)') 'Doing some analysis... Might be very slow.'
  !
  allocate(reacHeats(2, chem_net%nReactions), tmp(chem_net%nReactions))
  !
  write(fname_pre, &
        '(I4.4, "_rz_", F0.6, "_", F0.6, "_iter_", I3.3)') &
        id, leaves%list(id)%p%xmin, leaves%list(id)%p%ymin, &
        leaves%list(id)%p%iIter
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), 999999, getu=1)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chemsol_params%n_record_real
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i)
  end do
  close(fU3)
  !
  call openFileSequentialWrite(fU1, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'ele_'//trim(fname_pre)//'.dat'), 999, getu=1)
  !
  call openFileSequentialWrite(fU2, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'contri_'//trim(fname_pre)//'.dat'), 999, getu=1)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = 1+chemsol_params%n_record_real/20
  end if
  !
  write(fU1, '(2F10.1, 3ES12.2, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, &
    chem_params%Av_toStar, &
    chem_params%Av_toISM, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  write(fU2, '(2F10.1, 3ES12.2, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, &
    chem_params%Av_toStar, &
    chem_params%Av_toISM, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  do k=1, chemsol_params%n_record_real, a_disk_ana_params%ana_i_incr
    !
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chemsol_stor%record(:, k) - chemsol_stor%record(:, k-1))) / &
        (chemsol_stor%record(:, k) + chemsol_stor%record(:, k-1) + 1D-15))
      dt_t = (chemsol_stor%touts(k) - chemsol_stor%touts(k-1)) / &
             (chemsol_stor%touts(k) + chemsol_stor%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    chemsol_stor%y = chemsol_stor%record(:, k)
    chem_params%Tgas = chemsol_stor%y(chemsol_params%NEQ)
    !
    write(fU1, '("Time = ", ES14.4)') chemsol_stor%touts(k)
    write(fU1, '("Tgas = ", ES14.4)') chemsol_stor%y(chem_species%nSpecies+1)
    !
    call chem_elemental_residence
    !
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * &
        dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * &
        abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("Time = ", ES14.4)') chemsol_stor%touts(k)
    !
    if (ana_splist%nlen .le. 0) then
      cycle
    end if
    !
    call get_contribution_each
    !
    do i=1, chem_species%nSpecies
      if (.not. is_in_list_int(i, ana_splist%nlen, ana_splist%vals)) then
        cycle
      end if
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri) + 1D-100
      sum_dest = sum(chem_species%destr(i)%contri) + 1D-100
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, &
          chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%produ(i)%contri(j) .le. &
            chem_species%produ(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
      write(fU2, '(2X, A, 2X, ES12.2)') 'Destruction', sum_dest
      accum = 0D0
      do j=1, min(chem_species%destr(i)%nItem, 20)
        i0 = chem_species%destr(i)%list(j)
        accum = accum + chem_species%destr(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, &
          chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
    !
    write(fU2, '("Tgas = ", ES14.4)') chemsol_stor%y(chem_species%nSpecies+1)
    do i=1, chem_net%nReactions
      reacHeats(1, i) = dble(i)
    end do
    call  heating_chemical_termbyterm(chem_params%Tgas, chem_net%nReactions, tmp)
    reacHeats(2, :) = -abs(tmp)
    call quick_sort_array(reacHeats(:, :), 2, chem_net%nReactions, 1, (/2/))
    do i=1, min(chem_net%nReacWithHeat, 50)
        i0 = int(reacHeats(1, i))
        write(fU2, '(4X, I4, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          i, tmp(i0), chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
    end do
    !
  end do
  close(fU1)
  close(fU2)
  deallocate(reacHeats, tmp)
end subroutine chem_analyse



function get_H2_form_rate(c, XgH, XH, ngas) result(r)
  ! dn(H2)/dt
  double precision r
  double precision, intent(in) :: c, XgH, XH, ngas
  if (chemsol_params%H2_form_use_moeq) then
    r = c * XgH * XH * ngas
  else
    if (chem_idx_some_spe%i_gH .gt. 0) then
      r = c * XgH * XgH * ngas
    else
      r = c * XH * ngas
    end if
  end if
end function get_H2_form_rate



subroutine post_disk_iteration
  integer i
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      !TODO
      ! 1.
      !if (c%abundances(chem_idx_some_spe%i_H2O) .ge. 1D-6) then
      ! - Does not work
      !  c%using = .false.
      !  c%abundances = 0D0
      !end if
      ! 2.
      !if (c%par%Tgas .ge. 6D2) then
      ! - Does not work
      !  c%par%Tgas = 6D2
      !end if
      ! 3.
      ! - This reduces the water emission significantly.
      !c%par%Tgas = c%par%Tdust
      ! 4.
      ! Does not work
      !if (c%par%Tgas .ge. 5D2) then
      !  c%par%Tgas = c%par%Tdust
      !end if
      ! 5.
      ! It works
      !if (c%par%Tgas .ge. 2D2) then
      !  c%par%Tgas = c%par%Tdust
      !end if
      !
      !TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 2014-05-23 Fri 14:53:11
      ! ONLY FOR TESTING
      !c%par%Tgas = c%par%Tdust
      !c%par%sound_speed = sqrt(phy_kBoltzmann_CGS*c%par%Tgas / &
      !        (phy_mProton_CGS * c%par%MeanMolWeight*2D0))
      !c%par%velo_width_turb = c%par%sound_speed
      !c%abundances(chem_idx_some_spe%i_H2O) = &
      !  1D-2 * c%abundances(chem_idx_some_spe%i_H2O)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !if (c%par%sound_speed**2 .ge. &
      !  (2D0 * phy_GravitationConst_CGS * a_star%mass * phy_Msun_CGS &
      !   / ((c%xmax+c%xmin)*0.5D0 * phy_AU2cm))) then
      !  c%using = .false.
      !  if (allocated(c%abundances)) then
      !    c%abundances = 0D0
      !  end if
      !end if
      ! 2014-05-26 Mon 12:19:48
      !if (c%par%pressure_thermal .ge. &
      !    1D0 * abs(c%par%gravity_acc_z / c%par%area_T)) then
      !  c%using = .false.
      !  if (allocated(c%abundances)) then
      !    c%abundances = 0D0
      !  end if
      !end if
      ! 2014-05-26 Mon 23:50:47
      !c%par%n_gas = c%par%n_gas * &
      !  abs(c%par%gravity_acc_z / c%par%area_T) / c%par%pressure_thermal
      !c%par%Tgas = c%par%Tdust
      ! 2014-06-03 Tue 23:15:50
      !if (c%par%Ncol_toISM .le. 1D19) then
      !  c%abundances = 0D0
      !  c%using = .false.
      !end if
      ! 2014-06-09 Mon 11:03:49
      !if (c%xmax .le. 20D0) then
      !  c%par%Tgas = c%par%Tdust
      !end if
      ! 2014-06-09 Mon 12:55:54
      !if ((c%par%Av_toISM .lt. 1D0) .or.(c%par%Av_toStar .lt. 1D0)) then
      !  c%par%Tgas = c%par%Tdust + 1d2
      !else
      !  c%par%Tgas = c%par%Tdust
      !end if
      ! 2014-06-09 Mon 13:42:09
      !if (c%quality .ne. 0) then
      !  c%abundances = 0D0
      !end if
      !if (c%par%Tgas .ge. 5d2) then
      !  c%par%Tgas = c%par%tdust
      !end if
      ! 2014-06-12 Thu 10:01:28
      !if (c%xmax .le. 30D0) then
      !  c%par%n_gas = c%par%n_gas * 1D-2
      !end if
      ! 2014-07-08 Tue 01:15:28
      !if (c%xmax .le. 4D0) then
      !  c%par%n_gas = c%par%n_gas * 1D-1
      !end if
      ! 2014-07-08 Tue 11:27:38
      !if (c%xmax .le. 4D0) then
      !  c%par%n_gas = c%par%n_gas * 1D-2
      !end if
      ! 2014-07-08 Tue 12:19:16
      !if (c%xmax .le. 4D0) then
      !  c%par%Tgas = c%par%Tdust
      !end if
      ! 2014-07-11 Fri 18:10:45
      !if ((c%par%Av_toStar .ge. 1D0) .and. (c%par%Av_toISM .ge. 1D0)) then
      !  c%par%Tgas = c%par%Tdust
      !else
      !  c%par%Tgas = c%par%Tdust + 135D0
      !end if
      ! 2014-07-15 Tue 00:48:09
      !if (c%xmax .le. 4D0) then
      !  c%abundances = 0D0
      !end if
      ! 2014-07-27 Sun 01:40:02
      !if (c%xmax .le. 4D0) then
      !  c%par%Tgas = c%par%Tdust
      !end if
      ! 2014-09-05 Fri 14:32:29
      !if ((c%xmin .ge. 30D0) .and. (c%xmax .le. 60D0)) then
      !  c%abundances = c%abundances * 0.1D0
      !end if
      ! 2014-09-16 Tue 00:44:42
      !if (c%xmin .le. 10D0) then
      !  c%abundances = 0D0
      !end if
      ! 2014-11-07 Fri 23:37:25
      !if (c%xmin .le. 50D0) then
      !  c%par%n_gas = 0D0
      !  c%par%ndust_tot = 0D0
      !  c%par%n_dusts = 0D0
      !  c%par%rho_dusts = 0D0
      !  c%par%mdust_tot = 0D0
      !  c%par%mdusts_cell = 0D0
      !  c%abundances = 0D0
      !  c%optical%ext_tot = 0D0
      !  c%optical%flux = 0D0
      !end if
      ! 2015-04-11 Sat 17:48:04
      !c%par%Tgas = c%par%Tdust
      ! 2016-03-08 Tue 15:10:57
      !if (c%xmin .le. 20D0) then
      !  c%par%n_gas = c%par%n_gas * 0.01D0
      !else
      !  c%par%n_gas = c%par%n_gas * 0.3D0
      !end if
      ! 2016-08-01 Mon 13:55:40
      !if (c%xmin .le. 17D0) then
      !  c%par%n_gas = c%par%n_gas * 6D0
      !end if
      ! 2018-07-03 Tue 11:40:57
      ! c%par%Tgas = c%par%Tgas + 50D0
      !
    end associate
  end do
end subroutine post_disk_iteration



subroutine do_save_only_structure
  integer i
  ! Reload just in case the grid has changed
  !
  call disk_save_results_pre
  write(*, *) 'Saving structure to file: ', filename_save_results
  do i=1, leaves%nlen
    call update_params_above_alt(i)
    call disk_save_results_write(fU_save_results, leaves%list(i)%p)
  end do
  close(fU_save_results)
end subroutine do_save_only_structure




subroutine do_rerun_single_points
  integer i
  ! Reload just in case the grid has changed
  call load_ana_snippet
  !
  do i=1, ana_ptlist%nlen
    call calc_this_cell(ana_ptlist%vals(i))
  end do
end subroutine do_rerun_single_points


subroutine load_ana_snippet
  call load_ana_species_list
  call load_ana_points_list
  call get_species_produ_destr
  a_disk_ana_params%analyse_out_dir = &
    trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'ana/'))
  if (.not. dir_exist(a_disk_ana_params%analyse_out_dir)) then
    call my_mkdir(a_disk_ana_params%analyse_out_dir)
  end if
end subroutine load_ana_snippet



subroutine solve_a_Tdust(j)
  integer, intent(in) :: j
  double precision Td0, Ts1, Ts2, dTd, k, tmp
  integer, parameter :: nitermax_dusts=32
  integer i, i1
  ! Using Newton-Raphson iteration
  tmp = 1D0 / (4*phy_Pi*hc_params%mdusts_cell(j))
  do i=1, nitermax_dusts
    Td0 = hc_params%Tdusts(j)
    call update_en_exchange_with_dust
    Ts1 = get_Tdust_from_LUT( &
            (hc_params%en_gains(j) + hc_params%en_exchange(j)) * tmp, &
            luts%list(j), i1)
    dTd = 1D-2*hc_params%Tdusts(j) + 1D-1
    hc_params%Tdusts(j) = hc_params%Tdusts(j) + dTd
    call update_en_exchange_with_dust
    Ts2 = get_Tdust_from_LUT( &
            (hc_params%en_gains(j) + hc_params%en_exchange(j)) * tmp, &
            luts%list(j), i1)
    if (i1 .le. 0) then
      write(*, '(A, I6, 3ES16.6)') 'Error in solving Tdust:', &
        j, Ts1, Ts2, hc_params%Tgas
    end if
    k = (Ts2 - Ts1) / dTd
    if (abs(1D0 - k) .le. 1D-15) then
      hc_params%Tdusts(j) = Td0
      return
    else
      hc_params%Tdusts(j) = (Ts1 - k * Td0) / (1D0 - k)
    end if
    if (abs(Td0 - hc_params%Tdusts(j)) .le. 1D-5*Td0) then
      return
    end if
  end do
end subroutine solve_a_Tdust


function get_dEmit_dTd(j) result(dEmit_dTd)
  integer, intent(in) :: j
  double precision :: dEmit_dTd
  double precision :: Ts1, Ts2, E_Emit, dEmit, tmp
  double precision, parameter :: del_frac = 1D-2
  integer i1
  tmp = 1D0 / (4*phy_Pi*hc_params%mdusts_cell(j))
  E_Emit = hc_params%en_gains(j) + hc_params%en_exchange(j)
  dEmit = E_Emit * del_frac
  Ts1 = get_Tdust_from_LUT(E_Emit * tmp, luts%list(j), i1)
  Ts2 = get_Tdust_from_LUT((E_Emit + dEmit) * tmp, luts%list(j), i1)
  dEmit_dTd = dEmit / (Ts2 - Ts1) / hc_params%volume
end function get_dEmit_dTd


end module disk



subroutine chem_ode_f(NEQ, t, y, ydot)
  use chemistry
  use heating_cooling
  implicit none
  integer NEQ, i, j, i1
  double precision t, y(NEQ), ydot(NEQ), rtmp, tmp, tmp1
  ydot = 0D0
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    chem_params%Tgas = y(chem_species%nSpecies+1)
    call chem_cal_rates
  end if
  !
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 6, 21, 64) ! A + B -> C ! 53
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(2, i))
        if ((y(chem_net%reac(1, i)) .lt. 0D0) .and. &
            (y(chem_net%reac(2, i)) .lt. 0D0)) then
          rtmp = -rtmp
        end if
      case (1, 2, 3, 13, 61, 20) ! A -> B
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case (62)
        tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain
        if (tmp1 .le. 0D0) then
          rtmp = chem_net%rates(i)
        else
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-4) then
            rtmp = chem_net%rates(i) * tmp
          else
            rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
          end if
        end if
      case (75)
        tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain * chem_net%ABC(3, i)
        if (tmp1 .le. 0D0) then
          rtmp = chem_net%rates(i)
        else
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-4) then
            rtmp = chem_net%rates(i) * tmp
          else
            rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
          end if
        end if
      case (63) ! gA + gA -> gB
        ! dt(N(H2)) = k_HH * <H(H-1)>
        ! Moment equation:
        ! dt(N(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * N(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * X(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * X(H) * X(gH) * n_dust / D2G
        ! Rate equation:
        ! dt(X(H2)) = k_HH * X(H)**2 / D2G
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          if (chemsol_params%H2_form_use_moeq) then
            i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
            rtmp = chem_net%rates(i) * y(i1) * y(chem_net%reac(1, i))
            ydot(i1) = ydot(i1) - rtmp ! It's like H + gH -> gH2. So dt(H) -= rtmp, dt(gH) += rtmp
            ydot(chem_net%reac(1, i)) = ydot(chem_net%reac(1, i)) + rtmp
          else
            rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
          end if
        else
          rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
        end if
        if (y(chem_net%reac(1, i)) .lt. 0D0) then
          rtmp = -rtmp
        end if
      case (0)
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case default
        cycle
    end select
    !
    do j=1, chem_net%n_reac(i)
      ydot(chem_net%reac(j, i)) = ydot(chem_net%reac(j, i)) - rtmp
    end do
    do j=1, chem_net%n_prod(i)
      ydot(chem_net%prod(j, i)) = ydot(chem_net%prod(j, i)) + rtmp
    end do
  end do
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    call realtime_heating_cooling_rate(ydot(chem_species%nSpecies+1), NEQ, y)
  else
    ydot(chem_species%nSpecies+1) = 0D0
  end if
  !
end subroutine chem_ode_f




subroutine realtime_heating_cooling_rate(r, NEQ, y)
  ! Heating/cooling rate for a single average gas particle.
  use chemistry
  use heating_cooling
  use disk
  implicit none
  double precision, intent(out) :: r
  integer, intent(in) :: NEQ
  double precision, dimension(NEQ), intent(in) :: y
  integer j
  !
  hc_params%X_H2    = y(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = y(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = y(chem_idx_some_spe%i_CI)
  hc_params%X_CII   = y(chem_idx_some_spe%i_CII)
  hc_params%X_OI    = y(chem_idx_some_spe%i_OI)
  hc_params%X_NII   = y(chem_idx_some_spe%i_NII)
  hc_params%X_FeII  = y(chem_idx_some_spe%i_FeII)
  hc_params%X_SiII  = y(chem_idx_some_spe%i_SiII)
  hc_params%X_CO    = y(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = y(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = y(chem_idx_some_spe%i_OH)
  hc_params%X_E     = y(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = y(chem_idx_some_spe%i_Hplus)
  hc_params%X_Heplus= y(chem_idx_some_spe%i_Heplus)
  if (chem_idx_some_spe%i_gH .gt. 0) then
    hc_params%X_gH   = y(chem_idx_some_spe%i_gH)
  end if
  hc_params%R_H2_form_rate = &
    get_H2_form_rate( &
      hc_params%R_H2_form_rate_coeff, &
      hc_params%X_gH, &
      hc_params%X_HI, &
      hc_params%n_gas)
  !
  hc_params%Tgas = y(chem_species%nSpecies+1)
  hc_Tgas        = y(chem_species%nSpecies+1)
  !
  if (a_disk%allow_gas_dust_en_exch .and. a_disk%Tdust_iter_tandem &
      .and. (.not. isnan(hc_Tgas))) then
    !
    do j=1, dusts%n
      if ((hc_params%n_dusts(j) .le. 1D-20) .or. &
          (isnan(hc_params%en_exchange(j)))) then
        cycle
      end if
      if (heating_cooling_config%dust_gas_linear_couple) then
        hc_params%dEmit_dTd(j) = get_dEmit_dTd(j)
        !hc_params%Tdusts(j) = cooling_gas_grain_collision(Td_real)
      else
        call solve_a_Tdust(j)
      end if
    end do
    !
    hc_params%Tdust = &
        sum(hc_params%n_dusts * a_disk%dustcompo(:)%mrn%r2av * hc_params%Tdusts) / &
        sum(hc_params%n_dusts * a_disk%dustcompo(:)%mrn%r2av)
    hc_params%Tdust = max(hc_params%Tdust, a_disk_iter_params%minimum_Tdust)
  end if
  !
  hc_Tdust = hc_params%Tdust
  !
  !if (hc_Tdust .ge. mc_conf%TdustMax) then
  !  write(*, '(A, 2ES16.6)') 'Tdust too high.  Tdust,Tgas: ', hc_Tdust, hc_Tgas
  !  write(*, '(A, 4ES16.6)') 'Tdusts: ', hc_params%Tdusts
  !  write(*, '(A, 4ES16.6)') 'n_dusts: ', hc_params%n_dusts
  !  write(*, '(A, 4ES16.6)') 'en_gains: ', hc_params%en_gains
  !  write(*, '(A, 4ES16.6)') 'en_exchange: ', hc_params%en_exchange
  !end if
  !
  !hc_params%grand_gas_abundance = &
  !  sum(y(1:chem_species%nSpecies)) - sum(y(chem_species%idxGrainSpecies))
  !
  call get_alpha_viscosity_alt(hc_params, y, NEQ, a_disk%base_alpha)
  !
  r = heating_minus_cooling() * phy_SecondsPerYear / &
      (hc_params%n_gas * phy_kBoltzmann_CGS)
end subroutine realtime_heating_cooling_rate




subroutine chem_ode_jac(NEQ, t, y, j, ian, jan, pdj)
  use chemistry
  use heating_cooling
  use trivials
  implicit none
  double precision t, rtmp, tmp, tmp1, tmp2
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, i, j, k, i1
  double precision dT_dt_1, dT_dt_2, del_ratio, del_ratio_T, del_0, del_0_T, delta_y
  double precision, dimension(NEQ) :: ydot1, ydot2
  !
  del_ratio = 1D-2
  del_0 = chem_params%ratioDust2HnucNum*1D-6
  del_ratio_T = 1D-2
  del_0_T = 1D0
  !
  pdj = 0D0
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 6, 21, 64) ! A + B -> C
        if (j .EQ. chem_net%reac(1, i)) then
          if (chem_net%reac(1, i) .ne. chem_net%reac(2, i)) then
            rtmp = chem_net%rates(i) * y(chem_net%reac(2, i))
          else
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(2, i))
          end if
        else if (j .EQ. chem_net%reac(2, i)) then
          if (chem_net%reac(1, i) .ne. chem_net%reac(2, i)) then
            rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
          else
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
          end if
        else
          rtmp = 0D0
        end if
        if ((y(chem_net%reac(1, i)) .lt. 0D0) .and. &
            (y(chem_net%reac(2, i)) .lt. 0D0)) then
          rtmp = -rtmp
        end if
      case (1, 2, 3, 13, 61, 20) ! A -> B
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          rtmp = chem_net%rates(i)
        end if
      case (62)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp2 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain
          if (tmp2 .le. 0D0) then
            rtmp = 0D0
          else
            tmp1 = 1D0 / tmp2
            tmp = y(chem_net%reac(1, i)) * tmp1
            if (tmp .le. 1D-4) then
              rtmp = chem_net%rates(i) * tmp1
            else
              rtmp = chem_net%rates(i) * tmp1 * exp(-tmp)
            end if
          end if
        end if
      case (75)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp2 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain * chem_net%ABC(3, i)
          if (tmp2 .le. 0D0) then
            rtmp = 0D0
          else
            tmp1 = 1D0 / tmp2
            tmp = y(chem_net%reac(1, i)) * tmp1
            if (tmp .le. 1D-4) then
              rtmp = chem_net%rates(i) * tmp1
            else
              rtmp = chem_net%rates(i) * tmp1 * exp(-tmp)
            end if
          end if
        end if
      case (63) ! gA + gA -> gB
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          if (chemsol_params%H2_form_use_moeq) then
            i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
            if (j .eq. chem_net%reac(1, i)) then
              rtmp = chem_net%rates(i) * y(i1)
              pdj(i1) = pdj(i1) - rtmp
              pdj(chem_net%reac(1, i)) = pdj(chem_net%reac(1, i)) + rtmp
            else if (j .eq. i1) then
              rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
              pdj(i1) = pdj(i1) - rtmp
              pdj(chem_net%reac(1, i)) = pdj(chem_net%reac(1, i)) + rtmp
            else
              rtmp = 0D0
            end if
          else
            if (j .eq. chem_net%reac(1, i)) then
              rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
            else
              rtmp = 0D0
            end if
          end if
        else
          if (j .eq. chem_net%reac(1, i)) then
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
          else
            rtmp = 0D0
          end if
        end if
        if (y(chem_net%reac(1, i)) .lt. 0D0) then
          rtmp = -rtmp
        end if
      case (0)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          rtmp = chem_net%rates(i)
        end if
      case default
        cycle
    end select
    !
    if (rtmp .NE. 0D0) then
      do k=1, chem_net%n_reac(i)
        pdj(chem_net%reac(k, i)) = pdj(chem_net%reac(k, i)) - rtmp
      end do
      do k=1, chem_net%n_prod(i)
        pdj(chem_net%prod(k, i)) = pdj(chem_net%prod(k, i)) + rtmp
      end do
    end if
  end do
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (is_in_list_int(j, chem_idx_some_spe%nItem, chem_idx_some_spe%idx)) then
      if (y(j) .ge. 0D0) then
        call realtime_heating_cooling_rate(dT_dt_1, NEQ, y)
        delta_y = y(j) * del_ratio + del_0
        rtmp = y(j)
        y(j) = y(j) + delta_y
        call realtime_heating_cooling_rate(dT_dt_2, NEQ, y)
        pdj(chem_species%nSpecies+1) = (dT_dt_2 - dT_dt_1) / delta_y
        y(j) = rtmp
      else
        pdj(chem_species%nSpecies+1) = 0D0
      end if
    else if (j .eq. (chem_species%nSpecies+1)) then
      call chem_ode_f(NEQ, t, y, ydot1)
      delta_y = y(j) * del_ratio_T + del_0_T
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call chem_ode_f(NEQ, t, y, ydot2)
      pdj = (ydot2 - ydot1) / delta_y
      y(j) = rtmp
    end if
  else
    pdj(chem_species%nSpecies+1) = 0D0
  end if
end subroutine chem_ode_jac
