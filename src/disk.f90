module disk

use data_struct
use grid
use chemistry
use heating_cooling
use montecarlo

implicit none

type :: type_disk_basic_info
  double precision star_luminosity_in_Lsun
  double precision star_mass_in_Msun, star_radius_in_Rsun, star_temperature
  double precision disk_mass_in_Msun
  double precision ratio_uv2total
  double precision ratio_lyman2uv
  double precision ratio_xray2total
  double precision Lyman_phlumi_star_surface, &
                   UV_cont_phlumi_star_surface, &
                   Xray_phlumi_star_surface
  character(len=32) filename_exe
  logical            :: backup_src = .true.
  character(len=128) :: backup_src_cmd = &
    'find *.f90 *.f *.py makefile | cpio -pdm --insecure '
  double precision :: geometric_factor_UV   = 0.01D0
  double precision :: geometric_factor_Xray = 0.001D0
  !
  double precision :: dust2gas_mass_bg = 1D-5
  type(type_Andrews_disk) andrews_gas, andrews_dust, andrews_dust_bg
  !
  !double precision :: colDen2Av_coeff = 1D-21 ! Sun Kwok, eq 10.21
  !double precision :: colDen2Av_coeff = 5.3D-22 ! Draine 2011, eq 21.7
end type type_disk_basic_info


type :: type_disk_iter_params
  integer :: n_iter=128, n_iter_used = 0
  integer :: nlocal_iter = 2
  !
  double precision :: rtol_T = 0.1D0,    atol_T = 2D0
  double precision :: rtol_abun = 0.2D0, atol_abun = 1D-12
  !
  logical flag_converged
  integer n_cell_converged
  real converged_cell_percentage_stop
  !
  logical :: redo_montecarlo = .true.
  logical :: flag_save_rates = .false.
  logical :: flag_shortcut_ini = .false.
  !
  integer :: nSpecies_check_refine = 0
  integer :: ncell_refine = 0, count_refine = 0
  integer :: nMax_refine = 2
  double precision :: threshold_ratio_refine = 10D0
  character(len=128) filename_list_check_refine
  !
  character(len=128) iter_files_dir
  character(len=128) :: notes1 = ''
  character(len=128) :: notes2 = ''
  character(len=128) :: notes3 = ''
  character(len=128) :: notes4 = ''
end type type_disk_iter_params


type :: type_simple_integer_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: vals
end type type_simple_integer_list


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

! Basic config params of the disk model
type(type_disk_basic_info) a_disk

! Initialization params for the disk model
type(type_disk_basic_info) disk_params_ini

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

! Columns of cells
type(type_leaves), dimension(:), allocatable :: columns

! Index list of the cells that are being calculated
integer, dimension(:), allocatable :: calculating_cells
integer n_calculating_cells, n_calculating_cells_max

! Filename for saving the results for each iteration; the name will be changed
! for different iterations.
character(len=128) :: filename_save_results
integer fU_save_results

! Energy of a typical X-ray particle, in kev.
! For calculating the X-ray number flux.
! Should find a better way to deal with this.
double precision, parameter, private :: xray_energy_kev = 1D0

! Index of species that are to be used for check whether cell refinement is
! needed.
integer, dimension(:), allocatable, private :: idx_Species_check_refine
double precision, dimension(:), allocatable, private :: &
    thr_Species_check_refine

! For displaying some text to the screen
character(len=256) str_disp

! Field length for certain output
integer, parameter :: len_item=14

! Namelists for reading config params
namelist /disk_configure/ &
  disk_params_ini

namelist /cell_configure/ &
  cell_params_ini

namelist /iteration_configure/ &
  a_disk_iter_params

namelist /analyse_configure/ &
  a_disk_ana_params


contains



subroutine disk_iteration
  use my_timer
  type(date_time) a_date_time
  integer i, i0, i_count, l_count, ii
  !
  call disk_iteration_prepare
  !
  call montecarlo_prep
  !
  ! call montecarlo_reset_cells
  ! !
  ! call montecarlo_do(mc_conf, root)
  ! !
  ! call post_montecarlo
  ! !
  ! call openFileSequentialWrite(ii, &
  !   combine_dir_filename(mc_conf%mc_dir_out, 'cmp_radmc.dat'), 999)
  ! do i=1, leaves%nlen
  !   write(ii, '(2F9.2, 2I11, 18ES14.6)') &
  !     leaves%list(i)%p%par%Tdust1, leaves%list(i)%p%par%Tdust, &
  !     leaves%list(i)%p%optical%ab_count_dust, &
  !     leaves%list(i)%p%optical%cr_count, &
  !     leaves%list(i)%p%optical%en_gain_dust, &
  !     leaves%list(i)%p%optical%en_gain_abso, &
  !     leaves%list(i)%p%par%n_gas, &
  !     leaves%list(i)%p%par%mdust_cell, &
  !     leaves%list(i)%p%par%Ncol_toISM, &
  !     leaves%list(i)%p%par%Ncol_toStar, &
  !     leaves%list(i)%p%par%flux_UV, &
  !     leaves%list(i)%p%par%flux_Lya, &
  !     leaves%list(i)%p%par%dir_UV_r, &
  !     leaves%list(i)%p%par%dir_UV_z, &
  !     leaves%list(i)%p%par%aniso_UV, &
  !     leaves%list(i)%p%par%dir_Lya_r, &
  !     leaves%list(i)%p%par%dir_Lya_z, &
  !     leaves%list(i)%p%par%aniso_Lya, &
  !     leaves%list(i)%p%par%rmin, &
  !     leaves%list(i)%p%par%rmax, &
  !     leaves%list(i)%p%par%zmin, &
  !     leaves%list(i)%p%par%zmax
  ! end do
  ! close(ii)
  ! return
  !
  call save_post_config_params
  !
  ! Now start the major big loop.
  !
  a_disk_iter_params%count_refine = 0
  !
  write(str_disp, '("! ", A)') "Iteration begins."
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  do ii = 1, a_disk_iter_params%n_iter
    !
    a_disk_iter_params%n_iter_used = ii
    !
    write(str_disp, '("! ", A)') "Monte Carlo begins."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call montecarlo_reset_cells
    !
    call montecarlo_do(mc_conf, root)
    !
    ! Retrieve physical parameters from the monte carlo results
    call post_montecarlo
    !
    write(str_disp, '("! ", A)') "Monte Carlo finished."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! Write header to the file
    call disk_save_results_pre
    !
    ! Calculate layer by layer.
    ! Start from the surface layer.
    n_calculating_cells = surf_cells%nlen
    calculating_cells(1:surf_cells%nlen) = surf_cells%idx
    i_count = 0 ! Counter for cells
    l_count = 0 ! Counter for layers
    do
      l_count = l_count + 1
      !
      do i=1, n_calculating_cells
        i_count = i_count + 1
        i0 = calculating_cells(i)
        !
        write(*, '(3(A, I5, A, I5, ",", 2X), (A, I4, ","), 2X, A, 4F8.3)') &
          "Iter:", a_disk_iter_params%n_iter_used, "/", a_disk_iter_params%n_iter, &
          "Cell:", i_count, '/', leaves%nlen, &
          "cell:", i, '/', n_calculating_cells, &
          "Layer:", l_count, &
          'rz:', &
          leaves%list(i0)%p%par%rmin, &
          leaves%list(i0)%p%par%rmax, &
          leaves%list(i0)%p%par%zmin, &
          leaves%list(i0)%p%par%zmax
        write(*, '(2(A, ES10.3, 2X), 2X, 2A, 2X, 2A)') &
          'n_gas: ', leaves%list(i0)%p%par%n_gas, &
          'Tdust: ', leaves%list(i0)%p%par%Tdust, &
          'exe: ', trim(a_disk%filename_exe), &
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
        a_iter_stor%abundances(:, i0) = &
          leaves%list(i0)%p%abundances(chem_idx_some_spe%idx)
        !
        call disk_save_results_write(fU_save_results, leaves%list(i0)%p)
        flush(fU_save_results)
      end do
      !
      call update_calculating_cells
      !
      if (n_calculating_cells .eq. 0) then
        exit
      end if
    end do
    !
    write(str_disp, '("! ", A, I4, A)') "Iteration ", ii, " finished."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! At this point all the layers have been walked through.
    call check_convergency_whole_disk
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
    !
    else if (a_disk_iter_params%redo_montecarlo) then
      !if (ii .ge. min(4, a_disk_iter_params%n_iter/3)) then
      if (mod(ii, 7) .eq. 6) then
        ! Adjust the vertical structure
        call vertical_pressure_gravity_balance
        write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
        call display_string_both(str_disp, a_book_keeping%fU)
        !
      end if
      cycle
    !
    else if (a_disk_iter_params%count_refine .gt. &
             a_disk_iter_params%nMax_refine) then
      write(str_disp, '(A, I4, " > ", I4)') &
        '! Will not refine any more. count_refine: ', &
        a_disk_iter_params%count_refine, a_disk_iter_params%nMax_refine
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    !
    else
      write(*, '(/A)') 'Doing refinements where necessary.'
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
        call remake_index
        !
        call load_ana_points_list ! Reload, actually
        !
        if (allocated(a_iter_stor%T_s)) then
          deallocate(a_iter_stor%T_s, a_iter_stor%abundances)
        end if
        allocate(a_iter_stor%T_s(leaves%nlen), &
                 a_iter_stor%abundances(chem_idx_some_spe%nItem, &
                                                leaves%nlen))
        do i=1, leaves%nlen
          a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
          a_iter_stor%abundances(:,i) = &
            leaves%list(i)%p%abundances(chem_idx_some_spe%idx)
        end do
        !
        if (allocated(calculating_cells)) then
          deallocate(calculating_cells)
        end if
        n_calculating_cells_max = leaves%nlen
        allocate(calculating_cells(n_calculating_cells_max))
        !
        write(str_disp, '("!", A, 2X, I5)') 'New number of cells (leaf):', leaves%nlen
        call display_string_both(str_disp, a_book_keeping%fU)
        write(str_disp, '("!", A, 2X, I5)') 'New number of cells (total):', root%nOffspring
        call display_string_both(str_disp, a_book_keeping%fU)
      else
        write(str_disp, '("! ", A)') "No further refinement needed."
        call display_string_both(str_disp, a_book_keeping%fU)
      end if
    end if
  end do
  !
  if (a_disk_iter_params%flag_converged) then
    write(*, '(A/)') "Iteration has converged!"
  else
    write(*, '(A/)') "Iteration hasn't converged. :("
  end if
  !
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=iteration_configure)
  end if
  write(str_disp, '("!Final number of cells =", I4)') leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  ! call disk_iteration_postproc
  !
end subroutine disk_iteration




subroutine post_montecarlo
  integer i
  integer i1, i2
  integer, parameter :: cr_TH = 10
  double precision vx, vy, vz
  double precision RR
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      !if (c%optical%cr_count .ge. cr_TH) then
        c%par%Tdust1 = &
          get_Tdust_from_LUT(c%optical%en_gain_abso / &
                             (4*phy_Pi*c%par%mdust_cell), lut_0, i1)
        c%par%Tdust = c%par%Tdust1
      !else
      !end if
      !
      ! Flux of each cell as a function of wavelength
      c%optical%flux = c%optical%flux * (phy_AU2cm / c%par%volume)
      !
      ! Get some properties of the radiation field
      !
      ! UV
      i1 = max(1, get_idx_for_kappa(lam_range_UV(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_UV(2), dust_0))
      c%par%flux_UV = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      c%par%dir_UV_r = vx
      c%par%dir_UV_z = vz
      c%par%aniso_UV = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! Lya
      i1 = max(1, get_idx_for_kappa(lam_range_LyA(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_LyA(2), dust_0))
      c%par%flux_Lya = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      c%par%dir_Lya_r = vx
      c%par%dir_Lya_z = vz
      c%par%aniso_Lya = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! NIR
      i1 = max(1, get_idx_for_kappa(lam_range_NIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_NIR(2), dust_0))
      c%par%flux_NIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      c%par%dir_NIR_r = vx
      c%par%dir_NIR_z = vz
      c%par%aniso_NIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! MIR
      i1 = max(1, get_idx_for_kappa(lam_range_MIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_MIR(2), dust_0))
      c%par%flux_MIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      c%par%dir_MIR_r = vx
      c%par%dir_MIR_z = vz
      c%par%aniso_MIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! FIR
      i1 = max(1, get_idx_for_kappa(lam_range_FIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_FIR(2), dust_0))
      c%par%flux_FIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      c%par%dir_FIR_r = vx
      c%par%dir_FIR_z = vz
      c%par%aniso_FIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! Local number flux of Lyman alpha
      c%par%phflux_Lya = c%par%flux_Lya / phy_LyAlpha_energy_CGS
      c%par%G0_Lya_atten = c%par%flux_Lya / phy_Habing_energy_flux_CGS
      !
      ! Calculate the total column density to the star and to the ISM
      call calc_Ncol_to_ISM(leaves%list(i)%p)
      call calc_Ncol_to_Star(leaves%list(i)%p)
      !
      RR = (c%par%rcen**2 + c%par%zcen**2) * phy_AU2cm**2
      c%par%flux_UV_star_unatten = star_0%lumi_UV0 / (4D0*phy_Pi*RR)
      !
      ! Calculate the G0 factors
      ! The G0 is the unattenuated one, so a further
      ! exp(-k*Av) should be applied.
      c%par%G0_UV_toStar = c%par%flux_UV_star_unatten / phy_Habing_energy_flux_CGS
      c%par%G0_UV_toISM  = c%par%UV_G0_factor_background
      !
      c%par%Av_toStar = max(0D0, &
        -log(c%par%flux_UV / c%par%flux_UV_star_unatten) / phy_UVext2Av)
      ! The Av to ISM is a simple scaling of the dust column density
      c%par%Av_toISM = 1.086D0 * (phy_Pi * c%par%GrainRadius_CGS**2 * 2D0) * &
                       calc_Ncol_from_cell_to_point(c, c%par%rcen, root%ymax*2D0, -5)
    end associate
  end do
end subroutine post_montecarlo



subroutine montecarlo_reset_cells
  integer i
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      !
      c%par%Tdust1 = 0D0
      !
      c%par%X_HI = c%abundances(chem_idx_some_spe%i_HI)
      c%par%X_H2O = c%abundances(chem_idx_some_spe%i_H2O)
      !
      call calc_Ncol_to_ISM(leaves%list(i)%p)
      call calc_Ncol_to_Star(leaves%list(i)%p)
      !
      call allocate_local_optics(leaves%list(i)%p, gl_coll_0, dust_0)
      call reset_local_optics(leaves%list(i)%p)
    end associate
  end do
  !
end subroutine montecarlo_reset_cells



subroutine disk_iteration_prepare
  integer i
  !
  ! The density structure is needed for making the grid.
  a_andrews_4ini = disk_params_ini%andrews_gas
  a_andrews_4ini%particlemass = 1.4D0 * phy_mProton_CGS
  !
  call make_grid
  !
  n_calculating_cells_max = leaves%nlen
  allocate(calculating_cells(n_calculating_cells_max))
  !
  ! Prepare the chemical stuff
  call chem_read_reactions()
  call chem_load_reactions()
  call chem_parse_reactions()
  call chem_get_dupli_reactions()
  call chem_get_idx_for_special_species()
  call load_species_enthalpies
  call get_reaction_heat
  !
  call load_refine_check_species
  !
  chemsol_params%fU_log = a_book_keeping%fU
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !
  write(str_disp, '("!", A, 2X, I5)') 'Number of cells (leaf):', leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of cells (total):', root%nOffspring
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of reactions:', chem_net%nReactions
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of species ', chem_species%nSpecies
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  ! Set the disk and cell parameters
  call disk_set_disk_params
  call disk_set_gridcell_params
  call make_columns
  !
  if (.NOT. allocated(a_iter_stor%T_s)) then
    allocate(a_iter_stor%T_s(leaves%nlen), &
             a_iter_stor%abundances(chem_idx_some_spe%nItem, &
                                    leaves%nlen))
  end if
  !
  do i=1, leaves%nlen
    leaves%list(i)%p%abundances = chemsol_stor%y(1:chem_species%nSpecies)
    a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
    a_iter_stor%abundances(:, i) = chemsol_stor%y(chem_idx_some_spe%idx)
  end do
  !
  call disk_calc_disk_mass
  !
  call heating_cooling_prepare
  !
  if (a_disk_ana_params%do_analyse) then
    call load_ana_species_list
    call load_ana_points_list
    call get_species_produ_destr
    a_disk_ana_params%analyse_out_dir = &
      trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'ana/'))
    if (.not. dir_exist(a_disk_ana_params%analyse_out_dir)) then
      call my_mkdir(a_disk_ana_params%analyse_out_dir)
    end if
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
  star_0%mass   = a_disk%star_mass_in_Msun
  star_0%radius = a_disk%star_radius_in_Rsun
  star_0%T      = a_disk%star_temperature
  !
end subroutine disk_iteration_prepare



subroutine calc_this_cell(id)
  integer, intent(in) :: id
  integer j
  double precision tmp
  double precision R3, Z
  !
  leaves%list(id)%p%iIter = a_disk_iter_params%n_iter_used
  !
  do j=1, a_disk_iter_params%nlocal_iter
    !
    write(*, '("Local iter: ", I4, " of ", I4)') j, a_disk_iter_params%nlocal_iter
    !
    ! Set the initial condition for chemical evolution
    call set_initial_for_chem_solver(id, j)
    !
    write(*, '(4X, A, F12.3/)') 'Tgas_old: ', leaves%list(id)%p%par%Tgas
    !
    call update_params_above_alt(id)
    !
    call set_chemistry_params_from_cell(id)
    call chem_cal_rates
    !
    call set_heatingcooling_params_from_cell(id)
    !
    chemsol_params%flag_chem_evol_save = .false.
    !
    call chem_set_solver_flags_alt(j)
    if (j .eq. 1) then
      chemsol_params%evolT = .true.
      chemsol_params%maySwitchT = .false.
    else if (abs(heating_cooling_rates%hc_net_rate) .le. &
            1D-4 * min(max_heating_rate(), max_cooling_rate())) then
      chemsol_params%evolT = .false.
      chemsol_params%maySwitchT = .true.
    else
      chemsol_stor%y(chem_species%nSpecies+1) = &
        (0.5D0 + dble(j)*0.1D0) * &
        (leaves%list(id)%p%par%Tgas + leaves%list(id)%p%par%Tdust)
      chemsol_params%evolT = .true.
      chemsol_params%maySwitchT = .false.
    end if
    !
    call chem_evol_solve
    !
    leaves%list(id)%p%abundances = chemsol_stor%y(1:chem_species%nSpecies)
    leaves%list(id)%p%par%Tgas = chemsol_stor%y(chem_species%nSpecies+1)
    leaves%list(id)%p%quality = chemsol_params%quality
    leaves%list(id)%p%par%t_final = chemsol_stor%touts(chemsol_params%n_record_real)
    !
    write(*, '(4X, A, F12.3)') 'Tgas_new: ', leaves%list(id)%p%par%Tgas
    !
    !if (minval(chemsol_stor%y(1:chem_species%nSpecies)) .lt. 0D0) then
    !  if (leaves%list(id)%p%quality .eq. 0) then
    !    leaves%list(id)%p%quality = -1
    !  else
    !    leaves%list(id)%p%quality = -leaves%list(id)%p%quality
    !  end if
    !end if
    !
    call update_params_above_alt(id)
    !
    leaves%list(id)%p%par%pressure_thermal = &
      leaves%list(id)%p%par%n_gas * &
      (leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI) + &
       leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)) * &
      leaves%list(id)%p%par%Tgas * phy_kBoltzmann_CGS
    !
    ! R3 = R^3 = (sqrt(r^2 + z^2))^3
    R3 = (sqrt((leaves%list(id)%p%par%rcen)**2 + &
               (leaves%list(id)%p%par%zcen)**2))**3
    Z = leaves%list(id)%p%par%zcen
    leaves%list(id)%p%par%gravity_z = &
      phy_GravitationConst_CGS * (star_0%mass * phy_Msun_CGS) * &
      (leaves%list(id)%p%par%mgas_cell + &
       leaves%list(id)%p%par%mdust_cell) * &
      (-Z / R3 / (phy_AU2cm**2))
    !
    if (leaves%list(id)%p%quality .eq. 0) then
      exit
    end if
    !
  end do
  !
  a_iter_stor%T_s(id) = leaves%list(id)%p%par%Tgas
  !
  if (.not. chemsol_params%evolT) then
    call chem_cal_rates
    call realtime_heating_cooling_rate(tmp, chemsol_params%NEQ, chemsol_stor%y)
  end if
  leaves%list(id)%p%h_c_rates = heating_cooling_rates
  !
  !if (a_disk_iter_params%flag_save_rates) then
  !  call save_chem_rates(id)
  !end if
  !
  if (a_disk_ana_params%do_analyse) then
    if ((leaves%list(id)%p%quality .gt. 0) .or. &
        ((a_disk_iter_params%n_iter_used .gt. 1) .and. &
         (is_in_list_int(id, ana_ptlist%nlen, ana_ptlist%vals) .or. &
          need_to_refine(leaves%list(id)%p)))) then
      call chem_analyse(id)
    end if
  end if
end subroutine calc_this_cell



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
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiH2O) * const_LyAlpha_cross_H2O)))
    c%par%f_selfshielding_toStar_H2O = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiH2O) * const_LyAlpha_cross_H2O)))
    !
    c%par%f_selfshielding_toISM_OH = &
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiOH) * const_LyAlpha_cross_OH)))
    c%par%f_selfshielding_toStar_OH = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiOH) * const_LyAlpha_cross_OH)))
    !
    c%par%f_selfshielding_toISM_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toISM(chem_idx_some_spe%iiH2), &
      c%col_den_toISM(chem_idx_some_spe%iiCO))))
    !
    c%par%f_selfshielding_toStar_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toStar(chem_idx_some_spe%iiH2), &
      c%col_den_toStar(chem_idx_some_spe%iiCO))))
    !
    ! Calculate the gravitational force from above
    if (c%above%n .eq. 0) then
      c%par%gravity_acc_z = &
        phy_GravitationConst_CGS * (star_0%mass * phy_Msun_CGS) * &
          (calc_Ncol_from_cell_to_point(leaves%list(i0)%p, &
                                        c%par%rcen, root%ymax*2D0, -4) * &
           c%par%area_T * phy_mProton_CGS * c%par%MeanMolWeight) * &
          (-c%ymax / (sqrt(c%xmax**2 + c%ymax**2))**3 / (phy_AU2cm**2))
    else
      c%par%gravity_acc_z = &
        calc_Ncol_from_cell_to_point(leaves%list(i0)%p, &
                                     c%par%rcen, root%ymax*2D0, -1)
    end if
  end associate
end subroutine update_params_above_alt



!subroutine update_params_above(i0)
!  use load_Visser_CO_selfshielding
!  integer, intent(in) :: i0
!  integer j, j0
!  associate(p  => leaves%list(i0)%p, &
!            dz => leaves%list(i0)%p%par%dz * phy_AU2cm)
!    p%col_den   = p%abundances(chem_idx_some_spe%idx) * p%par%n_gas * dz
!    p%par%dNcol = p%par%n_gas * dz
!    p%col_den_acc = 0D0
!    p%par%Ncol = 0D0
!    if (leaves%list(i0)%p%above%n .gt. 0) then
!      do j=1, leaves%list(i0)%p%above%n
!        j0 = leaves%list(i0)%p%above%idx(j)
!        p%col_den_acc = p%col_den_acc + &
!                        (leaves%list(j0)%p%col_den_acc + &
!                         leaves%list(j0)%p%col_den) * p%above%fra(j)
!        p%par%Ncol    = p%par%Ncol + &
!                        (leaves%list(j0)%p%par%Ncol + &
!                         leaves%list(j0)%p%par%dNcol) * p%above%fra(j)
!      end do
!    end if
!  end associate
!  associate(p        => leaves%list(i0)%p, &
!            dz       => leaves%list(i0)%p%par%dz * phy_AU2cm, &
!            Ncol_H2  => leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiH2), &
!            dcol_H2  => leaves%list(i0)%p%col_den(chem_idx_some_spe%iiH2), &
!            Ncol_H   => leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiHI), &
!            dcol_H   => leaves%list(i0)%p%col_den(chem_idx_some_spe%iiHI), &
!            Ncol_H2O => leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiH2O), &
!            dcol_H2O => leaves%list(i0)%p%col_den(chem_idx_some_spe%iiH2O), &
!            Ncol_OH  => leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiOH), &
!            dcol_OH  => leaves%list(i0)%p%col_den(chem_idx_some_spe%iiOH), &
!            Ncol_CO  => leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiCO), &
!            dcol_CO  => leaves%list(i0)%p%col_den(chem_idx_some_spe%iiCO))
!    ! Kwok eq 10.20
!    p%par%Av = 1.086D0 * p%par%ratioDust2HnucNum * &
!      (phy_Pi * p%par%GrainRadius_CGS**2) * 2D0 * &
!      p%par%Ncol
!      !(p%par%Ncol + p%par%dNcol * 0.5D0)
!    p%par%f_selfshielding_H2  = &
!      min(1D0, get_H2_self_shielding(Ncol_H2, p%par%velo_width_turb))
!      !min(1D0, get_H2_self_shielding(Ncol_H2 + dcol_H2*0.5D0, p%par%velo_width_turb))
!      !min(1D0, ((Ncol_H2 + dcol_H2*0.5D0)/1D14)**(-0.75D0)) ! Tielens 2005, equation 8.39
!    p%par%f_selfshielding_H2O = &
!      min(1D0, exp(-(Ncol_H2O * const_LyAlpha_cross_H2O))) !* &
!      !tau2beta(dcol_H2O * const_LyAlpha_cross_H2O)
!    p%par%f_selfshielding_OH  = &
!      min(1D0, exp(-(Ncol_OH * const_LyAlpha_cross_OH))) !* &
!      !tau2beta(dcol_OH * const_LyAlpha_cross_OH)
!    p%par%f_selfshielding_CO = min(1D0, max(0D0, get_12CO_shielding(Ncol_H2, Ncol_CO)))
!  end associate
!end subroutine update_params_above


function get_H2_self_shielding(N_H2, dv_turb)
  ! Draine 1996, equation 37
  double precision get_H2_self_shielding
  double precision, intent(in) :: N_H2, dv_turb
  double precision x, b5
  x = N_H2 / 5D14
  b5 = dv_turb / 1D5
  get_H2_self_shielding = 0.965D0 / (1D0 + x/b5)**2 + &
    0.035 / sqrt(1D0 + x) * exp(-8.5D-4 * sqrt(1D0 + x))
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
  write(str_disp, '("! Iter", I4, " Number of cells converged: ", I6, "/", I6)') &
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


subroutine set_initial_for_chem_solver(id, iloc_iter)
  integer, intent(in) :: id, iloc_iter
  double precision tmp
  logical found_neighbor
  integer i, i0, ntmp
  !
  if (a_disk_iter_params%flag_shortcut_ini) then
    if (a_disk_iter_params%n_iter_used .eq. 1) then
      found_neighbor = .false.
      do i=1, leaves%list(id)%p%around%n
        i0 = leaves%list(id)%p%around%idx(i)
        if (leaves%list(i0)%p%iIter .gt. leaves%list(id)%p%iIter) then
          chemsol_stor%y(1:chem_species%nSpecies) = leaves%list(i0)%p%abundances
          found_neighbor = .true.
          exit
        end if
      end do
      if (.not. found_neighbor) then
        chemsol_stor%y(1:chem_species%nSpecies) = &
          chemsol_stor%y0(1:chem_species%nSpecies)
      end if
    else
      chemsol_stor%y(1:chem_species%nSpecies) = leaves%list(id)%p%abundances
    end if
    if (chemsol_params%neutralize) then
      tmp = sum(chemsol_stor%y(1:chem_species%nSpecies) * &
                        dble(chem_species%elements(1,:)))
      if (abs(tmp) .ge. 1D-2*chemsol_stor%y(chem_idx_some_spe%i_E)) then
        chemsol_stor%y(1:chem_species%nSpecies) = &
          chemsol_stor%y0(1:chem_species%nSpecies)
      else
        chemsol_stor%y(chem_idx_some_spe%i_E) = &
          chemsol_stor%y(chem_idx_some_spe%i_E) + tmp
        if (chemsol_stor%y(chem_idx_some_spe%i_E) .lt. 0D0) then
          ! When it is not possible to neutralize the composition by artificially
          ! changing the electron abundance, then use the general initial abundances,
          ! which should be absolutely neutral.
          chemsol_stor%y(1:chem_species%nSpecies) = &
            chemsol_stor%y0(1:chem_species%nSpecies)
          write(str_disp, '("! Cannot neutralize: X(E-) = ", ES12.4)') &
            chemsol_stor%y(chem_idx_some_spe%i_E)
          call display_string_both(str_disp, a_book_keeping%fU)
          write(str_disp, '("! Use y0 as initial abundance.")')
          call display_string_both(str_disp, a_book_keeping%fU)
          write(str_disp, '("! x, y = ", 2ES10.2, " iIter = ", I4)') &
            leaves%list(id)%p%xmin, leaves%list(id)%p%ymin, leaves%list(id)%p%iIter
          call display_string_both(str_disp, a_book_keeping%fU)
        end if
      end if
    end if
  else
    chemsol_stor%y(1:chem_species%nSpecies) = &
      chemsol_stor%y0(1:chem_species%nSpecies)
  end if
  !
  ! Initial abundance of *neutral* dust
  ! There should be no dust in the input initial abundance file
  chemsol_stor%y(chem_idx_some_spe%i_Grain0) = leaves%list(id)%p%par%ratioDust2HnucNum
  !
  ! Always use the temperature of the above cell as init
  if ((a_disk_iter_params%n_iter_used .eq. 1) .and. (iloc_iter .eq. 1)) then
    if (leaves%list(id)%p%above%n .gt. 0) then
      leaves%list(id)%p%par%Tgas = 0D0
      ntmp = 0
      do i=1, leaves%list(id)%p%above%n
        i0 = leaves%list(id)%p%above%idx(i)
        if (leaves%list(i0)%p%iIter .lt. leaves%list(id)%p%iIter) then
          cycle
        end if 
        ntmp = ntmp + 1
        leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tgas + &
                                      leaves%list(i0)%p%par%Tgas
      end do
      if (ntmp .gt. 0) then
        leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tgas / &
                                      dble(ntmp)
      else
        leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tdust
      end if
    else
      leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tdust
    end if
  end if
  !
  ! Set the initial temeprature
  chemsol_stor%y(chem_species%nSpecies+1) = leaves%list(id)%p%par%Tgas
  !
end subroutine set_initial_for_chem_solver



subroutine vertical_pressure_gravity_balance
  integer i
  double precision dznew, pnew, nnew, vnew
  if (.not. grid_config%columnwise) then
    return
  end if
  !
  write(str_disp, '(A)') '! Adjusting the vertical structure.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  ! Calculate the new density and size of each cell
  do i=1, leaves%nlen
    pnew = -leaves%list(i)%p%par%gravity_acc_z / leaves%list(i)%p%par%area_T
    ! Not change too much in one go
    pnew = sqrt(pnew * leaves%list(i)%p%par%pressure_thermal)
    !
    nnew = pnew / &
      (leaves%list(i)%p%par%Tgas * phy_kBoltzmann_CGS * &
       (leaves%list(i)%p%abundances(chem_idx_some_spe%i_HI) + &
        leaves%list(i)%p%abundances(chem_idx_some_spe%i_H2)))
    vnew = leaves%list(i)%p%par%mgas_cell / nnew / &
      (phy_mProton_CGS * leaves%list(i)%p%par%MeanMolWeight)
    dznew = vnew / leaves%list(i)%p%par%area_T / phy_AU2cm
    !
    leaves%list(i)%p%par%n_dust = leaves%list(i)%p%par%n_dust * nnew / &
                                  leaves%list(i)%p%par%n_gas
    leaves%list(i)%p%par%n_gas = nnew
    leaves%list(i)%p%val(1) = nnew
    leaves%list(i)%p%par%volume = vnew
    leaves%list(i)%p%par%dz = dznew
    !leaves%list(i)%p%par%n_dust = leaves%list(i)%p%par%n_gas * &
    !    leaves%list(i)%p%par%ratioDust2HnucNum
    leaves%list(i)%p%par%pressure_thermal = pnew
  end do
  !
  ! Move the cells column by column from bottom to top.
  call shift_and_scale_above
  !
  ! Remake neigbor information
  call grid_make_neighbors
  !
  call grid_make_surf_bott
  !
  mc_conf%maxw = get_surf_max_angle()
  call get_mc_stellar_par(mc_conf)
  !
  write(str_disp, '(A, 2F9.4//)') 'New minw,maxw: ', mc_conf%minw, mc_conf%maxw
  call display_string_both(str_disp, a_book_keeping%fU)
  !
end subroutine vertical_pressure_gravity_balance



subroutine shift_and_scale_above
  type(type_cell), pointer :: cthis
  double precision ybelow
  integer i, ic
  !
  do ic=1, bott_cells%nlen ! Loop over the columns
    ybelow = columns(ic)%list(1)%p%ymin
    do i=1, columns(ic)%nlen
      cthis => columns(ic)%list(i)%p
      if (cthis%using) then
        cthis%ymin = ybelow
        cthis%ymax = ybelow + cthis%par%dz
        !
        cthis%par%zmin = cthis%ymin
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        ! dz is already set
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
        !
      else
        cthis%ymax = ybelow + (cthis%ymax - cthis%ymin)
        cthis%ymin = ybelow
      end if
      !
      root%ymax = max(root%ymax, cthis%ymax)
      !
      ybelow = cthis%ymax
      !
    end do
  end do
  !
  ! Align the top edge to the domain upper boundary
  do ic=1, bott_cells%nlen
    i = columns(ic)%nlen
    cthis => columns(ic)%list(i)%p
    if (cthis%ymax .le. root%ymax) then
      ! Rescale the density
      cthis%val(1) = cthis%val(1) * &
        (cthis%ymax - cthis%ymin) / (root%ymax - cthis%ymin)
      cthis%ymax = root%ymax
      !
      if (cthis%using) then
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        cthis%par%dz = cthis%par%zmax - cthis%par%zmin
        !
        cthis%par%n_dust = cthis%par%n_dust * cthis%val(1) / cthis%par%n_gas
        cthis%par%n_gas  = cthis%val(1)
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
      end if
    else
      write(*,'(A)') 'Should not have this case:'
      write(*,'(A/)') 'in shift_and_scale_above.'
    end if
  end do
end subroutine shift_and_scale_above


subroutine make_columns
  integer i, j
  type(type_cell), pointer :: cthis, cnext
  double precision length, r, z, eps
  integer dirtype
  logical found
  type(type_ray) ray
  !
  allocate(columns(bott_cells%nlen))
  do i=1, bott_cells%nlen
    cthis => leaves%list(bott_cells%idx(i))%p
    r = cthis%par%rcen
    z = root%ymax * 2D0
    columns(i)%nlen = int(calc_Ncol_from_cell_to_point(cthis, r, z, -2))
    if (.not. allocated(columns(i)%list)) then
      allocate(columns(i)%list(columns(i)%nlen))
    end if
    !
    ray%x = cthis%par%rcen
    ray%y = 0D0
    ray%z = cthis%par%zcen
    ray%vx = 0D0
    ray%vy = 0D0
    ray%vz = 1D0
    !
    j = 0
    ! First make a list of all the cells (including null ones) above the bottom
    ! (mid-plane) one.  They are ordered from bottom to top.
    do
      call calc_intersection_ray_cell(ray, cthis, &
        length, r, z, eps, found, dirtype)
      if (.not. found) then
        write(str_disp,'(A)') 'In make_columns:'
        call display_string_both(str_disp, a_book_keeping%fU)
        write(str_disp,'(A, 9ES16.6/)') 'ray does not intersect cthis: ', &
          sqrt(ray%x**2+ray%y**2), ray%z, ray%vx, ray%vy, ray%vz, &
            cthis%xmin, cthis%xmax, cthis%ymin, cthis%ymax
        call display_string_both(str_disp, a_book_keeping%fU)
        return
      end if
      !
      j = j + 1
      columns(i)%list(j)%p => cthis
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
  end do
end subroutine make_columns



subroutine calc_Ncol_to_ISM(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toISM(iSp) = calc_Ncol_from_cell_to_point( &
      c, c%par%rcen, root%ymax * 2D0, &
      chem_idx_some_spe%idx(iSp))
  else
    c%par%Ncol_toISM = calc_Ncol_from_cell_to_point( &
      c, c%par%rcen, root%ymax * 2D0)
  end if
end subroutine calc_Ncol_to_ISM



subroutine calc_Ncol_to_Star(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toStar(iSp) = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0, chem_idx_some_spe%idx(iSp))
  else
    c%par%Ncol_toStar = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0)
  end if
end subroutine calc_Ncol_to_Star



function calc_Ncol_from_cell_to_point(c, r, z, iSpe) result(N)
  double precision N
  type(type_cell), intent(in), target :: c ! Todo: pointer or not?
  double precision, intent(in) :: r, z
  integer, intent(in), optional :: iSpe
  type(type_ray) ray
  type(type_cell), pointer :: cthis, cnext
  double precision t, length, r1, z1, eps
  logical found
  integer dirtype
  double precision, parameter :: small_dist = 1D-50
  !
  N = 0D0
  !
  ray%x = c%par%rcen
  ray%y = 0D0
  ray%z = c%par%zcen
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
      write(str_disp,'(A, 9ES13.4/)') &
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
            N = N + cthis%par%n_dust * length * phy_AU2cm
          end if
        end if
      else
        write(str_disp,'(A/)') 'I do not know what to do!'
        call display_string_both(str_disp, a_book_keeping%fU)
        stop
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
  if (.NOT. getFileUnit(fU_save_results)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  write(filename_save_results, '("iter_", I4.4, ".dat")') &
    a_disk_iter_params%n_iter_used
  filename_save_results = trim(combine_dir_filename( &
    a_disk_iter_params%iter_files_dir, filename_save_results))
  call openFileSequentialWrite(fU_save_results, filename_save_results, 99999)
  !
  call write_header(fU_save_results)
end subroutine disk_save_results_pre


subroutine write_header(fU)
  integer, intent(in) :: fU
  character(len=64) fmt_str
  character(len=8192) tmp_str
  write(fmt_str, '("(", I4, "A14)")') chem_species%nSpecies
  write(tmp_str, fmt_str) chem_species%names
  write(fU, '(A)') &
    '!' // &
    str_pad_to_len('cvg', 4) // &
    str_pad_to_len('qual', 5) // &
    str_pad_to_len('arnd', 5) // &
    str_pad_to_len('abov', 5) // &
    str_pad_to_len('belo', 5) // &
    str_pad_to_len('innr', 5) // &
    str_pad_to_len('outr', 5) // &
    str_pad_to_len('cr_count',len_item) // &
    str_pad_to_len('abc_dus', len_item) // &
    str_pad_to_len('scc_HI',  len_item) // &
    str_pad_to_len('abc_wat', len_item) // &
    str_pad_to_len('t_final', len_item) // &
    str_pad_to_len('rmin',    len_item) // &
    str_pad_to_len('rmax',    len_item) // &
    str_pad_to_len('zmin',    len_item) // &
    str_pad_to_len('zmax',    len_item) // &
    str_pad_to_len('Tgas',    len_item) // &
    str_pad_to_len('Tdust',   len_item) // &
    str_pad_to_len('n_gas',   len_item) // &
    str_pad_to_len('n_dust',  len_item) // &
    str_pad_to_len('d2gmas',  len_item) // &
    str_pad_to_len('d2gnum',  len_item) // &
    str_pad_to_len('deplet',  len_item) // &
    str_pad_to_len('rd_min',  len_item) // &
    str_pad_to_len('rd_max',  len_item) // &
    str_pad_to_len('mrn_n',   len_item) // &
    str_pad_to_len('rd_av',   len_item) // &
    str_pad_to_len('r2d_av',  len_item) // &
    str_pad_to_len('r3d_av',  len_item) // &
    str_pad_to_len('mg_cell', len_item) // &
    str_pad_to_len('md_cell', len_item) // &
    str_pad_to_len('presr_t', len_item) // &
    str_pad_to_len('grav_z',  len_item) // &
    str_pad_to_len('egain_d', len_item) // &
    str_pad_to_len('egain_ab',len_item) // &
    str_pad_to_len('flx_UV',  len_item) // &
    str_pad_to_len('flx_Lya', len_item) // &
    str_pad_to_len('flx_NIR', len_item) // &
    str_pad_to_len('flx_MIR', len_item) // &
    str_pad_to_len('flx_FIR', len_item) // &
    str_pad_to_len('vr_UV',   len_item) // &
    str_pad_to_len('vz_UV',   len_item) // &
    str_pad_to_len('ani_UV',  len_item) // &
    str_pad_to_len('vr_Lya',  len_item) // &
    str_pad_to_len('vz_Lya',  len_item) // &
    str_pad_to_len('ani_Lya', len_item) // &
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
    str_pad_to_len('XRay0',   len_item) // &
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
    str_pad_to_len('c_wa_ro', len_item) // &
    str_pad_to_len('c_wa_vi', len_item) // &
    str_pad_to_len('c_CO_ro', len_item) // &
    str_pad_to_len('c_CO_vi', len_item) // &
    str_pad_to_len('c_H2_ro', len_item) // &
    str_pad_to_len('c_LyAlp', len_item) // &
    str_pad_to_len('c_fb   ', len_item) // &
    str_pad_to_len('c_ff   ', len_item) // &
    trim(tmp_str)
end subroutine write_header


subroutine disk_save_results_write(fU, c)
  character(len=64) fmt_str
  integer, intent(in) :: fU
  type(type_cell), pointer, intent(in) :: c
  integer converged
  !
  write(fmt_str, '(", ", I4, "ES14.4E4)")') chem_species%nSpecies
  if (c%converged) then
    converged = 1
  else
    converged = 0
  end if
  write(fU, '(7I5, 4I14, 95ES14.5E3' // trim(fmt_str)) &
  converged                                              , &
  c%quality                                              , &
  c%around%n                                             , &
  c%above%n                                              , &
  c%below%n                                              , &
  c%inner%n                                              , &
  c%outer%n                                              , &
  c%optical%cr_count                                     , &
  c%optical%ab_count_dust                                , &
  c%optical%sc_count_HI                                  , &
  c%optical%ab_count_water                               , &
  c%par%t_final                                          , &
  c%par%rmin                                             , &
  c%par%rmax                                             , &
  c%par%zmin                                             , &
  c%par%zmax                                             , &
  c%par%Tgas                                             , &
  c%par%Tdust                                            , &
  c%par%n_gas                                            , &
  c%par%n_dust                                           , &
  c%par%ratioDust2GasMass                                , &
  c%par%ratioDust2HnucNum                                , &
  c%par%dust_depletion                                   , &
  c%mrn%rmin                                             , &
  c%mrn%rmax                                             , &
  c%mrn%n                                                , &
  c%mrn%rav                                              , &
  c%mrn%r2av                                             , &
  c%mrn%r3av                                             , &
  c%par%mgas_cell                                        , &
  c%par%mdust_cell                                       , &
  c%par%pressure_thermal                                 , &
  c%par%gravity_acc_z                                    , &
  c%optical%en_gain_dust                                 , &
  c%optical%en_gain_abso                                 , &
  c%par%flux_UV                                          , &
  c%par%flux_Lya                                         , &
  c%par%flux_NIR                                         , &
  c%par%flux_MIR                                         , &
  c%par%flux_FIR                                         , &
  c%par%dir_UV_r                                         , &
  c%par%dir_UV_z                                         , &
  c%par%aniso_UV                                         , &
  c%par%dir_Lya_r                                        , &
  c%par%dir_Lya_z                                        , &
  c%par%aniso_Lya                                        , &
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
  c%par%Xray_flux_0                                      , &
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
  c%h_c_rates%cooling_Neufeld_H2O_rate_rot               , &
  c%h_c_rates%cooling_Neufeld_H2O_rate_vib               , &
  c%h_c_rates%cooling_Neufeld_CO_rate_rot                , &
  c%h_c_rates%cooling_Neufeld_CO_rate_vib                , &
  c%h_c_rates%cooling_Neufeld_H2_rot_rate                , &
  c%h_c_rates%cooling_LymanAlpha_rate                    , &
  c%h_c_rates%cooling_free_bound_rate                    , &
  c%h_c_rates%cooling_free_free_rate                     , &
  c%abundances
end subroutine disk_save_results_write


subroutine disk_calc_disk_mass
  integer i
  a_disk%disk_mass_in_Msun = 0D0
  do i=1, leaves%nlen
    associate(p => leaves%list(i)%p%par)
      a_disk%disk_mass_in_Msun = &
        a_disk%disk_mass_in_Msun + &
        p%n_gas * p%MeanMolWeight * (phy_2Pi * p%rcen * p%dr * p%dz)
    end associate
  end do
  a_disk%disk_mass_in_Msun = a_disk%disk_mass_in_Msun * &
    phy_AU2cm**3 * phy_mProton_CGS / phy_Msun_CGS
end subroutine disk_calc_disk_mass


subroutine set_heatingcooling_params_from_cell(id)
  integer id
  hc_params%type_cell_rz_phy_basic = leaves%list(id)%p%par
  hc_params%Neufeld_dv_dz = leaves%list(id)%p%par%velo_gradient * 1D-5 ! cm s-1 to km s-1
  hc_params%Neufeld_G     = 1D0
  hc_params%X_H2    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CI)
  hc_params%X_Cplus = leaves%list(id)%p%abundances(chem_idx_some_spe%i_Cplus)
  hc_params%X_OI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OI)
  hc_params%X_CO    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OH)
  hc_params%X_E     = leaves%list(id)%p%abundances(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = leaves%list(id)%p%abundances(chem_idx_some_spe%i_Hplus)
  hc_params%X_gH    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_gH)
  if (chemsol_params%H2_form_use_moeq) then
    hc_params%R_H2_form_rate = &
      hc_params%R_H2_form_rate_coeff * &
      hc_params%X_gH * &
      hc_params%X_HI * &
      hc_params%n_gas
  else
    hc_params%R_H2_form_rate = &
      hc_params%R_H2_form_rate_coeff * &
      hc_params%X_gH * &
      hc_params%X_gH * &
      hc_params%n_gas
  end if
  !
  leaves%list(id)%p%par%R_H2_form_rate = hc_params%R_H2_form_rate
end subroutine set_heatingcooling_params_from_cell


subroutine set_chemistry_params_from_cell(id)
  integer id
  chem_params => leaves%list(id)%p%par
end subroutine set_chemistry_params_from_cell


subroutine disk_set_a_cell_params(c, cell_params_copy)
  type(type_cell), target :: c
  type(type_cell_rz_phy_basic), intent(in) :: cell_params_copy
  if (.not. associated(c%par)) then
    allocate(c%par)
    allocate(c%mrn)
  end if
  if (.not. allocated(c%h_c_rates)) then
    allocate(c%h_c_rates)
  end if
  !
  if (.not. allocated(c%abundances)) then
    allocate(c%abundances(chem_species%nSpecies), &
             c%col_den_toISM(chem_idx_some_spe%nItem), &
             c%col_den_toStar(chem_idx_some_spe%nItem))
  end if
  !
  c%iIter = 0
  c%quality = 0
  !
  c%par = cell_params_copy
  !
  c%par%rmin = c%xmin
  c%par%rmax = c%xmax
  c%par%rcen = (c%xmax + c%xmin) * 0.5D0
  c%par%dr   = c%xmax - c%xmin
  !
  c%par%zmin = c%ymin
  c%par%zmax = c%ymax
  c%par%zcen = (c%ymax + c%ymin) * 0.5D0
  c%par%dz   = c%ymax - c%ymin
  !
  c%par%daz  = 0D0
  !
  c%par%volume = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * c%par%dz * phy_AU2cm**3
  c%par%area_T = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * phy_AU2cm**2
  c%par%area_B = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * phy_AU2cm**2
  c%par%area_I = phy_2Pi * c%par%rmin * c%par%dz * phy_AU2cm**2
  c%par%area_O = phy_2Pi * c%par%rmax * c%par%dz * phy_AU2cm**2
  c%par%surf_area = c%par%area_T + c%par%area_B + c%par%area_I + c%par%area_O
  !
  ! Get gas number density and gas mass in each cell
  a_disk%andrews_gas%particlemass = c%par%MeanMolWeight * phy_mProton_CGS
  c%par%n_gas = Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_gas)
  c%par%mgas_cell = c%par%n_gas * c%par%volume * (phy_mProton_CGS * c%par%MeanMolWeight)
  !
  ! Local dust size distribution; preliminary
  c%mrn%rmin = c%par%aGrainMin_micron ! micron
  c%mrn%rmax = c%par%aGrainMax_micron
  c%mrn%n    = c%par%mrn_ind
  call calc_dust_MRN_par(c%mrn)
  c%par%GrainRadius_CGS = sqrt(c%mrn%r2av) * 1D-4 ! = sqrt(<r**2>)
  !
  c%par%SitesPerGrain = 4D0 * phy_Pi * c%mrn%r2av * (1D-4)**2 * SitesDensity_CGS
  !
  ! Dust particle mass
  c%par%mdust = 4.0D0*phy_Pi/3.0D0 * (1D-4)**3 * c%mrn%r3av * &
                c%par%GrainMaterialDensity_CGS
  ! Here n_dust is actually the mass density
  ! Two dust component: a background one, and a settled one.
  c%par%n_dust = Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_dust) + &
                 Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_dust_bg)
  c%par%mdust_cell = c%par%n_dust * c%par%volume
  ! Now convert to particle density for dust
  c%par%n_dust = c%par%n_dust / c%par%mdust
  !
  c%par%ratioDust2GasMass = c%par%mdust_cell / c%par%mgas_cell
  !
  c%par%ratioDust2HnucNum = c%par%n_dust / c%par%n_gas
  c%par%dust_depletion = c%par%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
  !
  c%val(1) = c%par%n_gas
  !
  if (grid_config%use_data_file_input) then
    c%par%Tgas    = c%val(2)
    c%par%Tdust   = c%val(2)
  else
    c%par%Tgas    = 400D0 / (1D0 + c%par%rcen) * (1D0 + c%par%zcen)
    c%par%Tdust   = 0D0 ! instead of c%par%Tgas
  end if
  !
  c%par%pressure_thermal = 0D0
  c%par%gravity_z = 0D0
  c%par%gravity_acc_z = 0D0
  !
  ! c%par%UV_G0_factor = c%par%UV_G0_factor_background + &
  !   a_disk%UV_cont_phlumi_star_surface &
  !      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
  !      / phy_Habing_photon_flux_CGS &
  !      * a_disk%geometric_factor_UV
  ! c%par%LymanAlpha_number_flux_0 = &
  !   a_disk%Lyman_phlumi_star_surface &
  !      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2)
  ! c%par%LymanAlpha_energy_flux_0 = c%par%LymanAlpha_number_flux_0 * phy_LyAlpha_energy_CGS
  ! c%par%LymanAlpha_G0_factor = c%par%LymanAlpha_energy_flux_0 / phy_Habing_energy_flux_CGS
  c%par%Xray_flux_0 = &
    a_disk%Xray_phlumi_star_surface &
      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
      * a_disk%geometric_factor_Xray
  associate( &
          G     => phy_GravitationConst_CGS, &
          M     => a_disk%star_mass_in_Msun * phy_Msun_CGS, &
          r     => c%par%rcen * phy_AU2cm, &
          v     => c%par%velo_Kepler, &
          w     => c%par%omega_Kepler, &
          dv_dr => c%par%velo_gradient, &
          delv  => c%par%velo_width_turb, &
          l     => c%par%coherent_length)
    v = sqrt(G * M / r)
    w = v / r
    dv_dr = 0.5D0 * v / r
    delv = v ! Todo
    l = delv / dv_dr
  end associate
end subroutine disk_set_a_cell_params


subroutine disk_set_gridcell_params
  integer i
  do i=1, leaves%nlen
    call disk_set_a_cell_params(leaves%list(i)%p, cell_params_ini)
  end do
end subroutine disk_set_gridcell_params


subroutine disk_set_disk_params
  a_disk = disk_params_ini
  !
  ! Background dust
  a_disk%andrews_dust_bg = a_disk%andrews_gas
  a_disk%andrews_dust_bg%useNumDens = .false.
  a_disk%andrews_dust_bg%Md = &
    a_disk%andrews_gas%Md * a_disk%dust2gas_mass_bg
  !
  associate( &
    Lstar => a_disk%star_luminosity_in_Lsun * phy_Lsun_CGS, &
    uv2total => a_disk%ratio_uv2total, &
    lyman2uv => a_disk%ratio_lyman2uv, &
    xray2total => a_disk%ratio_xray2total)
    a_disk%UV_cont_phlumi_star_surface = &
      Lstar * uv2total * (1D0 - lyman2uv) / phy_UV_cont_energy_CGS
    a_disk%Lyman_phlumi_star_surface = &
      Lstar * uv2total * lyman2uv         / phy_LyAlpha_energy_CGS
    a_disk%Xray_phlumi_star_surface  = &
      Lstar * xray2total / (xray_energy_kev*1D3*phy_eV2erg)
    write(str_disp, '("!Stellar total luminosity = ", ES12.4, " erg s-1")') Lstar
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar UV cont luminosity = ", ES12.4, " erg s-1")') &
      a_disk%UV_cont_phlumi_star_surface * phy_UV_cont_energy_CGS
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar UV cont photon count rate = ", ES12.4, " s-1")') &
      a_disk%UV_cont_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar LyA luminosity = ", ES12.4, " erg s-1")') &
      a_disk%Lyman_phlumi_star_surface * phy_LyAlpha_energy_CGS
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar LyA photon count rate = ", ES12.4, " s-1")') &
      a_disk%Lyman_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar X-ray luminosity = ", ES12.4, " erg s-1")') &
      a_disk%Xray_phlumi_star_surface * (xray_energy_kev*1D3*phy_eV2erg)
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar X-ray photon count rate = ", ES12.4, " s-1")') &
      a_disk%Xray_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
  end associate
end subroutine disk_set_disk_params


function get_local_doppler_kepler_scale(M, r, dv, factor)
  double precision :: get_local_doppler_kepler_scale
  double precision, intent(in) :: M, r, dv
  double precision, optional :: factor
  if (.NOT. present(factor)) then
    factor = 1D0
  end if
  get_local_doppler_kepler_scale = factor * 2D0 * r * dv / &
    sqrt(phy_GravitationConst_CGS * M * phy_Msun_CGS / (r * phy_AU2cm))
end function get_local_doppler_kepler_scale


function get_local_dv_microturb(M, r, T)
  double precision :: get_local_dv_microturb
  double precision, intent(in), optional :: M, r, T
  if (.not. present(M)) then
    get_local_dv_microturb = 1D5 ! = 1 km s-1
  else
    get_local_dv_microturb = 1D5 ! = 1 km s-1
    write(*,'(A)') 'In get_local_dv_microturb:'
    write(*,'(A/)') 'Not implemented yet!'
  end if
end function get_local_dv_microturb


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
  type(type_heating_cooling_parameters) heat_cool_log
  ! Use namelist for output some logging infomation.
  ! Not very readable, but easy to implement.
  namelist /cell_par_log/ heat_cool_log
  !
  write(filename, '("reac_rates_cell_", I4.4, ".dat")') i0
  if(.NOT. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  dir = trim(combine_dir_filename(a_book_keeping%dir, 'rates_log/'))
  if (.NOT. dir_exist(dir)) then
    call my_mkdir(dir)
  end if
  call openFileSequentialWrite(fU, combine_dir_filename(dir, filename), 99999)
  !
  heat_cool_log = hc_params
  write(fU, nml=cell_par_log)
  !
  do k=1, chem_net%nReactions
    write(fU, '(A135, ES16.4E4)') chem_reac_str%list(k), chem_net%rates(k)
  end do
  close(fU)
end subroutine save_chem_rates


subroutine save_post_config_params
  type(type_disk_basic_info) disk_params_tmp
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
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialRead(fU, &
    combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
      a_disk_iter_params%filename_list_check_refine), 99)
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
          '(F7.1)') thr_Species_check_refine(i1)
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
      call display_string_both(str_disp, a_book_keeping%fU)
      call refine_this_cell_vertical(leaves%list(i)%p, n_refine)
    end if
  end do
end subroutine do_refine


subroutine remake_index
  call get_number_of_leaves(root)
  leaves%nlen = root%nleaves
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
end subroutine remake_index


function need_to_refine(c, n_refine)
  logical need_to_refine
  type(type_cell), target :: c
  integer, intent(out), optional :: n_refine
  integer i, i0, i1, j
  double precision val_max, val_min
  logical flag1, flag2
  flag1 = .false.
  flag2 = .false.
  if (present(n_refine)) then
    n_refine = 0
  end if
  if (c%par%dz .le. grid_config%smallest_cell_size) then
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
      call disk_set_a_cell_params(cc, c%par)
      cc%par%Tgas = c%par%Tgas
      !
      cc%h_c_rates = c%h_c_rates
      cc%abundances = c%abundances
      cc%col_den = c%col_den
      cc%col_den_acc = c%col_den_acc
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
  !deallocate(c%par, c%h_c_rates, c%abundances, c%col_den, c%col_den_acc)
  !deallocate(c%inner%idx, c%inner%fra)
  !deallocate(c%outer%idx, c%outer%fra)
  !deallocate(c%above%idx, c%above%fra)
  !deallocate(c%below%idx, c%below%fra)
  !deallocate(c%around%idx, c%around%fra)
  !deallocate(c%inner, c%outer, c%above, c%below, c%around)
end subroutine refine_this_cell_vertical


subroutine disk_iteration_postproc
  integer fU, fU1, fU2, ios, i, i0, j, idx, idx_diff
  double precision r, z
  double precision sum_prod, sum_dest, accum
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  call get_species_produ_destr
  !
  write(*,*) 'Trying to find out where are the elements.'
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_points), 99)
  if (.not. getFileUnit(fU1)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialWrite(fU1, &
       combine_dir_filename(a_book_keeping%dir, &
         a_disk_ana_params%file_analyse_res_ele), 999)
  if (.not. getFileUnit(fU2)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialWrite(fU2, &
       combine_dir_filename(a_book_keeping%dir, &
         a_disk_ana_params%file_analyse_res_contri), 999)
  do
    read(fU, '(2F6.2)', iostat=ios) r, z
    if (ios .ne. 0) then
      exit
    end if
    idx = 0
    do i=1, leaves%nlen
      if ((leaves%list(i)%p%par%rmin .le. r) .and. (leaves%list(i)%p%par%rmax .ge. r) .and. &
          (leaves%list(i)%p%par%zmin .le. z) .and. (leaves%list(i)%p%par%zmax .ge. z)) then
        idx = i
        exit
      end if
    end do
    if (idx .eq. 0) then
      write(*, '("Point (", 2F6.2, ")", A)') r, z, ' not in any cells!'
      cycle
    end if
    chemsol_stor%y(1:chem_species%nSpecies) = &
        leaves%list(idx)%p%abundances(1:chem_species%nSpecies)
    call chem_elemental_residence
    write(fU1, '("(", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    if (leaves%list(idx)%p%above%n .gt. 0) then
      idx_diff = leaves%list(idx)%p%above%idx(1)
    else if (leaves%list(idx)%p%below%n .gt. 0) then
      idx_diff = leaves%list(idx)%p%below%idx(1)
    else
      idx_diff = idx
    end if
    !
    call set_chemistry_params_from_cell(idx_diff)
    call chem_cal_rates
    call get_contribution_each
    !
    write(fU2, '("This (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av
    write(fU2, '("Diff (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') &
      leaves%list(idx_diff)%p%par%rcen, &
      leaves%list(idx_diff)%p%par%zcen, &
      leaves%list(idx_diff)%p%par%Tgas,  leaves%list(idx_diff)%p%par%Tdust, &
      leaves%list(idx_diff)%p%par%n_gas, leaves%list(idx_diff)%p%par%Ncol, &
      leaves%list(idx_diff)%p%par%Av
    do i=1, chem_species%nSpecies
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(A12, ": ", ES12.2, " Diff: ", ES12.2, " Rate: ", ES12.2)') chem_species%names(i), &
        chemsol_stor%y(i), leaves%list(idx_diff)%p%abundances(i), &
        sum_prod - sum_dest
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
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
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU)
  close(fU1)
  close(fU2)
end subroutine disk_iteration_postproc




subroutine load_ana_species_list
  integer fU, ios, i, n
  integer, dimension(:), allocatable :: list_tmp
  character(len=12) str
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in load_ana_species_list.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_species), 99)
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
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in load_ana_points_list.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_points), 99)
  allocate(list_tmp(leaves%nlen))
  n = 0
  do
    read(fU, '(2F6.2)', iostat=ios) r, z
    if (ios .ne. 0) then
      exit
    end if
    do i=1, leaves%nlen
      if ((leaves%list(i)%p%par%rmin .le. r) .and. (leaves%list(i)%p%par%rmax .ge. r) .and. &
          (leaves%list(i)%p%par%zmin .le. z) .and. (leaves%list(i)%p%par%zmax .ge. z)) then
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



subroutine load_species_enthalpies
  integer j, i1, fU, ios
  double precision dblTmp
  character(Len=32) FMTstr, strTMP
  character commentChar
  character(LEN=const_len_species_name) nameSpecies_tmp
  ! The output enthalpies are in K.
  if (allocated(chem_species%enthalpies)) then
    deallocate(chem_species%enthalpies)
  end if
  allocate(chem_species%enthalpies(chem_species%nSpecies))
  chem_species%enthalpies = dblNaN()
  i1 = 0
  commentChar = '!'
  if (IsWordChar(chemsol_params%filename_species_enthalpy(1:1))) then
    if (.NOT. getFileUnit(fU)) then
      write (*,*) 'In subroutine ImportSpeciesEnthalpy:'
      write (*,*) 'Cannot allocate an output file unit!'
      stop
    end if
    write (FMTstr, FMT= '("(", "A", I2, ", F", I1, ".0)")') &
      const_len_species_name, 9
    CALL openFileSequentialRead(fU, &
        combine_dir_filename(chemsol_params%chem_files_dir, &
        chemsol_params%filename_species_enthalpy), 999999)
    do
      read (UNIT=fU, FMT='(A32)', IOSTAT=ios) strTMP
      if (ios .NE. 0) then
        exit
      end if
      if ((strTMP(1:1) .EQ. commentChar) .OR. &
          (strTMP(1:1) .EQ. ' ')) then
        cycle
      end if
      read (strTMP, FMT=FMTstr, IOSTAT=ios) nameSpecies_tmp, dblTmp
      if (ios .NE. 0) then
        write (*, *) 'Error in importing enthalpies: ios = ', ios
        stop
      end if
      do j=1, chem_species%nSpecies
        if (trim(chem_species%names(j)) .EQ. trim(nameSpecies_tmp)) then
          ! Convert from kJ/mol to K to erg.
          chem_species%enthalpies(j) = dblTmp * 1D3 / phy_IdealGasConst_SI * phy_kBoltzmann_CGS
          i1 = i1 + 1
          exit
        end if
      end do
    end do
    close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
    write(str_disp, '("! ", I5, A, I5, A)') &
      i1, ' of ', chem_species%nSpecies, ' species have enthalpy.'
    call display_string_both(str_disp, a_book_keeping%fU)
  end if
end subroutine load_species_enthalpies



subroutine get_reaction_heat
  integer i, j
  if (allocated(chem_net%heat)) then
    deallocate(chem_net%heat)
  end if
  allocate(chem_net%heat(chem_net%nReactions))
  chem_net%heat = dblNaN()
  do i=1, chem_net%nReactions
    if (chem_net%itype(i) .ne. 5) then
      cycle
    end if
    if ((chem_net%ctype(i) .eq. 'RA') .or. (chem_net%ctype(i) .eq. 'RR')) then
      cycle
    end if
    chem_net%heat(i) = 0D0
    do j=1, chem_net%n_reac(i)
      chem_net%heat(i) = chem_net%heat(i) + chem_species%enthalpies(chem_net%reac(j,i))
    end do
    do j=1, chem_net%n_prod(i)
      chem_net%heat(i) = chem_net%heat(i) - chem_species%enthalpies(chem_net%prod(j,i))
    end do
  end do
end subroutine get_reaction_heat



subroutine chem_analyse(id)
  integer, intent(in) :: id
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=128) fname_pre
  character(len=32) FMTstryHistory
  double precision dy_y, dt_t
  double precision frac
  frac = 0.1D0
  !
  write(*, '(/A/)') 'Doing some analysis... Might be very slow.'
  !
  if (.not. getFileUnit(fU3)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  write(fname_pre, &
        '(I4.4, "_rz_", F0.6, "_", F0.6, "_iter_", I3.3)') &
        id, leaves%list(id)%p%xmin, leaves%list(id)%p%ymin, &
        leaves%list(id)%p%iIter
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), 999999)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chemsol_params%n_record
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i)
  end do
  close(fU3)
  !
  if (.not. getFileUnit(fU1)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  !
  call openFileSequentialWrite(fU1, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'ele_'//trim(fname_pre)//'.dat'), 999)
  !
  if (.not. getFileUnit(fU2)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  call openFileSequentialWrite(fU2, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'contri_'//trim(fname_pre)//'.dat'), 999)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = 1+chemsol_params%n_record/20
  end if
  !
  write(fU1, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  write(fU2, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
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
    write(fU1, '("time = ", ES14.4)') chemsol_stor%touts(k)
    !
    chemsol_stor%y = chemsol_stor%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chemsol_stor%touts(k)
    if (ana_splist%nlen .le. 0) then
      cycle
    end if
    call get_contribution_each
    do i=1, chem_species%nSpecies
      if (.not. is_in_list_int(i, ana_splist%nlen, ana_splist%vals)) then
        cycle
      end if
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
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
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU1)
  close(fU2)
end subroutine chem_analyse




subroutine a_test_case
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=64) FMTstryHistory, fname_pre
  !
  allocate(chem_params)
  chem_params = a_disk_ana_params%chempar
  !
  call chem_read_reactions
  call chem_load_reactions
  call chem_parse_reactions
  call chem_get_dupli_reactions
  call chem_get_idx_for_special_species
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !chemsol_stor%y(chem_idx_some_spe%i_E) = &
  !  chemsol_stor%y(chem_idx_some_spe%i_E) + &
  !    sum(chemsol_stor%y(1:chem_species%nSpecies) * &
  !        dble(chem_species%elements(1,:)))
  write(*, '(A)') 'Initial abundances:'
  do i=1, chem_species%nSpecies
    if (chemsol_stor%y(i) .ge. 1D-20) then
      write(*, '(I5, 2X, A, ES12.5)') i, chem_species%names(i), &
        chemsol_stor%y(i)
    end if
  end do
  write(*, '(//)')
  !
  associate(ch => chem_params)
    !ch%Tgas = 10D0
    !ch%Tdust = 10D0
    !ch%n_gas = 1D6
    !ch%UV_G0_factor = 0D0
    !ch%UV_G0_factor_background = 1D0
    !ch%Av = 10D0
    !ch%LymanAlpha_number_flux_0 = 0D0
    !ch%Xray_flux_0 = 0D0
    !ch%Ncol = 1D22
    !ch%dNcol = 1D21
    !ch%f_selfshielding_H2 = 0D0
    !ch%f_selfshielding_CO = 0D0
    !ch%f_selfshielding_H2O = 0D0
    !ch%f_selfshielding_OH = 0D0
    ch%GrainMaterialDensity_CGS = 2D0
    ch%ratioDust2GasMass = 0.01D0
    ch%MeanMolWeight = 1.4D0
    ch%ratioDust2HnucNum = &
          ch%ratioDust2GasMass * (phy_mProton_CGS * ch%MeanMolWeight) &
          / (4.0D0*phy_Pi/3.0D0 * (ch%GrainRadius_CGS)**3 * &
             ch%GrainMaterialDensity_CGS)
    ch%dust_depletion = ch%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
    ch%n_dust = ch%n_gas * ch%ratioDust2HnucNum
    chemsol_stor%y(chem_species%nSpecies+1) = ch%Tgas
  end associate
  !
  call chem_cal_rates
  call chem_set_solver_flags
  chemsol_params%evolT = .false.
  call chem_evol_solve
  !
  write(*,*) 'Doing some analysis... Might be very slow.'
  !
  call get_species_produ_destr
  !
  if (.not. getFileUnit(fU3)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_iter_params%iter_files_dir, 'func_of_time.dat'), 999999)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, &
    '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chemsol_params%n_record
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i)
  end do
  close(fU3)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = chemsol_params%n_record / 4
  end if
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
    write(fname_pre, '(I4.4, "_")') k
    !
    if (.not. getFileUnit(fU1)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
    call openFileSequentialWrite(fU1, &
      combine_dir_filename(a_disk_iter_params%iter_files_dir, &
        trim(fname_pre)//'elemental_residence.dat'), 999)
    !
    if (.not. getFileUnit(fU2)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
    call openFileSequentialWrite(fU2, &
      combine_dir_filename(a_disk_iter_params%iter_files_dir, &
        trim(fname_pre)//'contribution_reactions.dat'), 999)
    !
    chemsol_stor%y = chemsol_stor%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chemsol_stor%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    call get_contribution_each
    !
    write(fU2, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chemsol_stor%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    do i=1, chem_species%nSpecies
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, 3ES12.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
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
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, 3ES11.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
    close(fU1)
    close(fU2)
  end do
end subroutine a_test_case



subroutine b_test_case
  integer i, j, fU
  type(type_cell_rz_phy_basic), pointer :: ch => null()
  double precision Tmin, Tmax, dT, ratio
  double precision n_gas_min, n_gas_max, dn
  double precision h_c_net_rate
  character(len=128) filename, fname_pre, header
  type(type_cell), pointer :: c => null()
  !
  filename = 'Tgas_hc_abundances.dat'
  !
  allocate(ch)
  ch = a_disk_ana_params%chempar
  chem_params => ch
  !
  call chem_read_reactions
  call chem_load_reactions
  call chem_parse_reactions
  call chem_get_dupli_reactions
  call chem_get_idx_for_special_species
  call load_species_enthalpies
  call get_reaction_heat
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !
  call heating_cooling_prepare
  !
  call load_ana_species_list
  call get_species_produ_destr
  !
  Tmin = 1D2
  Tmax = 1200D0
  dT = 10D0
  n_gas_min = 3.5D5
  n_gas_max = 3.6D5
  dn = 1D3
  ratio = 1D0
  !
  ch%Tgas = Tmin
  ch%n_gas = n_gas_min
  !
  allocate(c)
  allocate(c%col_den_acc(chem_idx_some_spe%nItem), &
           c%col_den(chem_idx_some_spe%nItem), &
           c%abundances(chem_species%nSpecies))
  allocate(c%around, c%above, c%below, c%inner, c%outer)
  allocate(c%h_c_rates)
  !
  c%par => ch
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit!'
    stop
  end if
  call openFileSequentialWrite(fU, &
    combine_dir_filename(a_disk_iter_params%iter_files_dir, filename), 99999)
  !
  call write_header(fU)
  !
  a_disk_ana_params%analyse_out_dir = &
    trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'ana/'))
  if (.not. dir_exist(a_disk_ana_params%analyse_out_dir)) then
    call my_mkdir(a_disk_ana_params%analyse_out_dir)
  end if
  !
  do i=1, 1
    do j=1, 299
      !
      chemsol_params%evolT = .true.
      !
      ch%Tgas = ch%Tdust
      !
      ch%GrainMaterialDensity_CGS = 2D0
      ch%ratioDust2GasMass = 0.01D0
      ch%MeanMolWeight = 1.4D0
      ch%ratioDust2HnucNum = &
            ch%ratioDust2GasMass * (phy_mProton_CGS * ch%MeanMolWeight) &
            / (4.0D0*phy_Pi/3.0D0 * (ch%GrainRadius_CGS)**3 * &
               ch%GrainMaterialDensity_CGS)
      ch%dust_depletion = ch%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
      ch%n_dust = ch%n_gas * ch%ratioDust2HnucNum
      write(*,*) 'Dust density ', ch%n_dust
      write(*,*) ch%ratioDust2HnucNum
      write(*,*) ch%n_gas
      !
      ch%velo_Kepler = 30D5
      ch%omega_Kepler = ch%velo_Kepler / phy_AU2cm
      ch%velo_gradient = 0.5D0 * ch%velo_Kepler / phy_AU2cm
      ch%velo_width_turb = ch%velo_Kepler
      ch%coherent_length = ch%velo_width_turb / ch%velo_gradient
      !
      write(*,'(I4, F9.1, ES12.4, F9.1/)') i, ch%Tgas, ch%n_gas, ch%Tdust
      !
      chemsol_stor%y(1:chem_species%nSpecies) = chemsol_stor%y0(1:chem_species%nSpecies)
      chemsol_stor%y(chem_species%nSpecies+1) = ch%Tgas
      !
      call chem_cal_rates
      write(*,'(2ES12.2/)') chem_params%f_selfshielding_H2, chem_params%Av
      call chem_set_solver_flags_alt(1)
      !
      hc_params%type_cell_rz_phy_basic = ch
      !
      hc_params%Neufeld_dv_dz = 10D0/phy_AU2cm
      hc_params%Neufeld_G     = 1D0
      !
      hc_params%X_H2    = chemsol_stor%y(chem_idx_some_spe%i_H2)
      hc_params%X_HI    = chemsol_stor%y(chem_idx_some_spe%i_HI)
      hc_params%X_CI    = chemsol_stor%y(chem_idx_some_spe%i_CI)
      hc_params%X_Cplus = chemsol_stor%y(chem_idx_some_spe%i_Cplus)
      hc_params%X_OI    = chemsol_stor%y(chem_idx_some_spe%i_OI)
      hc_params%X_CO    = chemsol_stor%y(chem_idx_some_spe%i_CO)
      hc_params%X_H2O   = chemsol_stor%y(chem_idx_some_spe%i_H2O)
      hc_params%X_OH    = chemsol_stor%y(chem_idx_some_spe%i_OH)
      hc_params%X_E     = chemsol_stor%y(chem_idx_some_spe%i_E)
      hc_params%X_Hplus = chemsol_stor%y(chem_idx_some_spe%i_Hplus)
      hc_params%X_gH    = chemsol_stor%y(chem_idx_some_spe%i_gH)
      !
      if (chemsol_params%H2_form_use_moeq) then
        hc_params%R_H2_form_rate = &
          hc_params%R_H2_form_rate_coeff * &
          hc_params%X_gH * &
          hc_params%X_HI * &
          hc_params%n_gas
      else
        hc_params%R_H2_form_rate = &
          hc_params%R_H2_form_rate_coeff * &
          hc_params%X_gH * &
          hc_params%X_gH * &
          hc_params%n_gas
      end if
      ch%R_H2_form_rate = hc_params%R_H2_form_rate
      !
      !call realtime_heating_cooling_rate(tmp, chemsol_params%NEQ, chemsol_stor%y)
      !write(*,'(2ES16.6/)') tmp, heating_minus_cooling()
      !call disp_h_c_rates
      call chem_evol_solve
      !
      c%abundances  = chemsol_stor%y(1:chem_species%nSpecies)
      c%col_den     = c%abundances(chem_idx_some_spe%idx) * c%par%dNcol
      c%col_den_acc = c%abundances(chem_idx_some_spe%idx) * c%par%Ncol
      !
      hc_Tgas = ch%Tgas
      hc_Tdust = ch%Tdust
      h_c_net_rate = heating_minus_cooling()
      !
      c%h_c_rates = heating_cooling_rates
      c%par%t_final = chemsol_stor%touts(chemsol_params%n_record_real)
      !
      call disk_save_results_write(fU, c)
      !
      write(fname_pre, '(I4.4, "_", I4.4)') i, j
      write(header, '("n_gas = ", ES13.6)') ch%n_gas
      !
      a_disk_ana_params%ana_i_incr = 1
      call do_a_analysis(fname_pre, header)
      !
      ch%n_gas = ch%n_gas + dn
      dn = dn * ratio
      if (ch%n_gas .gt. n_gas_max) then
        exit
      end if
    end do
    ch%Tgas = ch%Tgas + dT
    dT = dT * ratio
    if (ch%Tgas .gt. Tmax) then
      exit
    end if
  end do
  close(fU)
  !
end subroutine b_test_case




subroutine do_a_analysis(fname_pre, header)
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=128), intent(in) :: fname_pre, header
  character(len=32) FMTstryHistory
  double precision dy_y, dt_t
  double precision frac
  double precision r
  frac = 0.1D0
  !
  write(*, '(/A/)') 'Doing some analysis... Might be slow.'
  !
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), maxRowLen=999999, getu=1)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 3
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, &
    '  Tgas        ', &
    '  hc          '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 3
  do i=1, chemsol_params%n_record
    call realtime_heating_cooling_rate(r, chemsol_params%NEQ, chemsol_stor%record(:, i))
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i), r
  end do
  close(fU3)
  !
  call openFileSequentialWrite(fU1, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'ele_'//trim(fname_pre)//'.dat'), maxRowLen=999, getu=1)
  call openFileSequentialWrite(fU2, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'contri_'//trim(fname_pre)//'.dat'), maxRowLen=999, getu=1)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = 1+chemsol_params%n_record/20
  end if
  !
  write(fU1, '(A)') trim(header)
  write(fU2, '(A)') trim(header)
  !
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
    !+++
    if ((chemsol_stor%touts(k) .le. 1D2) .or. (chemsol_stor%touts(k) .ge. 1D5)) then
      cycle
    end if
    if (chemsol_stor%touts(k) .ge. 4D2) then
      if (mod(k, 50) .ne. 0) then
        cycle
      end if
    end if
    !---
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chemsol_stor%record(1:chem_species%nSpecies, k) - &
             chemsol_stor%record(1:chem_species%nSpecies, k-1))) / &
        (chemsol_stor%record(1:chem_species%nSpecies, k) + &
         chemsol_stor%record(1:chem_species%nSpecies, k-1) + 1D-15))
      dt_t = (chemsol_stor%touts(k) - chemsol_stor%touts(k-1)) / &
             (chemsol_stor%touts(k) + chemsol_stor%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    write(fU1, '("time = ", ES14.4)') chemsol_stor%touts(k)
    !
    chemsol_stor%y(1:chem_species%nSpecies) = chemsol_stor%record(1:chem_species%nSpecies, k)
    chem_params%Tgas = chemsol_stor%record(chem_species%nSpecies+1, k)
    call chem_cal_rates
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chemsol_stor%touts(k)
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
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
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
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU1)
  close(fU2)
end subroutine do_a_analysis


end module disk



subroutine chem_ode_f(NEQ, t, y, ydot)
  use chemistry
  use heating_cooling
  implicit none
  integer NEQ, i, j, i1
  double precision t, y(NEQ), ydot(NEQ), rtmp, tmp
  ydot = 0D0
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (y(chem_species%nSpecies+1) .ne. chem_params%Tgas) then
      chem_params%Tgas = y(chem_species%nSpecies+1)
      call chem_cal_rates
    end if
  end if
  !
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 21, 64) ! A + B -> C ! 53
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(2, i))
      case (1, 2, 3, 13, 61, 0, 20) ! A -> B
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case (62)
        tmp = y(chem_net%reac(1, i)) / &
          (chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain)
        if (tmp .le. 1D-9) then
          rtmp = chem_net%rates(i) * tmp
        else
          rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
        end if
      case (75)
        tmp = y(chem_net%reac(1, i)) / &
          (chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain &
           * chem_net%ABC(3, i))
        if (tmp .le. 1D-9) then
          rtmp = chem_net%rates(i) * tmp
        else
          rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
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
  use chemistry
  use heating_cooling
  double precision, intent(out) :: r
  integer, intent(in) :: NEQ
  double precision, dimension(NEQ), intent(in) :: y
  hc_params%Tgas    = y(chem_species%nSpecies+1)
  hc_params%X_H2    = y(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = y(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = y(chem_idx_some_spe%i_CI)
  hc_params%X_Cplus = y(chem_idx_some_spe%i_Cplus)
  hc_params%X_OI    = y(chem_idx_some_spe%i_OI)
  hc_params%X_CO    = y(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = y(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = y(chem_idx_some_spe%i_OH)
  hc_params%X_E     = y(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = y(chem_idx_some_spe%i_Hplus)
  hc_params%X_gH    = y(chem_idx_some_spe%i_gH)
  hc_params%R_H2_form_rate_coeff = chem_params%R_H2_form_rate_coeff
  if (chemsol_params%H2_form_use_moeq) then
    hc_params%R_H2_form_rate = &
      hc_params%R_H2_form_rate_coeff * &
      hc_params%X_gH * &
      hc_params%X_HI * &
      hc_params%n_gas
  else
    hc_params%R_H2_form_rate = &
      hc_params%R_H2_form_rate_coeff * &
      hc_params%X_gH * &
      hc_params%X_gH * &
      hc_params%n_gas
  end if
  hc_Tgas = y(chem_species%nSpecies+1)
  hc_Tdust = hc_params%Tdust
  r = &
    heating_minus_cooling() * phy_SecondsPerYear / (chem_params%n_gas * phy_kBoltzmann_CGS)
end subroutine realtime_heating_cooling_rate




subroutine chem_ode_jac(NEQ, t, y, j, ian, jan, pdj)
  use chemistry
  use heating_cooling
  use trivials
  implicit none
  double precision t, rtmp, tmp, tmp1
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, i, j, k, i1
  double precision dT_dt_1, dT_dt_2, del_ratio, del_0, delta_y
  double precision, dimension(NEQ) :: ydot1, ydot2
  del_ratio = 1D-3
  del_0 = 1D-16
  pdj = 0D0
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 21, 64) ! A + B -> C
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
      case (1, 2, 3, 13, 61, 0, 20) ! A -> B
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          rtmp = chem_net%rates(i)
        end if
      case (62)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-9) then
            rtmp = chem_net%rates(i) / tmp1
          else
            rtmp = chem_net%rates(i) / tmp1 * exp(-tmp)
          end if
        end if
      case (75)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain &
                 * chem_net%ABC(3, i)
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-9) then
            rtmp = chem_net%rates(i) / tmp1
          else
            rtmp = chem_net%rates(i) / tmp1 * exp(-tmp)
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
      call realtime_heating_cooling_rate(dT_dt_1, NEQ, y)
      delta_y = y(j) * del_ratio + del_0
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call realtime_heating_cooling_rate(dT_dt_2, NEQ, y)
      pdj(chem_species%nSpecies+1) = (dT_dt_2 - dT_dt_1) / delta_y
      y(j) = rtmp
    else if (j .eq. (chem_species%nSpecies+1)) then
      call chem_ode_f(NEQ, t, y, ydot1)
      delta_y = y(j) * del_ratio * 1D1 + del_0
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
