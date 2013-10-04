module disk

use data_struct
use grid
use chemistry
use heating_cooling

implicit none

type :: phy_chem_rad_disk_params
  double precision star_luminosity_in_Lsun
  double precision star_mass_in_Msun, disk_mass_in_Msun
  double precision ratio_uv2total
  double precision ratio_lyman2uv
  double precision ratio_xray2total
  double precision Lyman_alpha_freq, Lyman_alpha_Lorentz_width, Lyman_alpha_crosssec_peak
  double precision Lyman_phlumi_star_surface, UV_cont_phlumi_star_surface, Xray_phlumi_star_surface
  character(len=32)   filename_exe
  logical          :: backup_src = .true.
  character(len=128) :: backup_src_cmd = 'find *.f90 *.f *.py makefile | cpio -pdm --insecure '
  !double precision :: colDen2Av_coeff = 1D-21 ! Sun Kwok, eq 10.21
  !double precision :: colDen2Av_coeff = 5.3D-22 ! Draine 2011, eq 21.7
  double precision :: geometric_factor_UV   = 0.01D0
  double precision :: geometric_factor_Xray = 0.001D0
end type phy_chem_rad_disk_params


type :: a__disk
  type(phy_chem_rad_disk_params) :: params
end type a__disk


type :: disk_iteration_params
  integer :: n_iter=128, n_iter_used = 0, ncell_refine = 0
  integer :: nlocal_iter = 2
  double precision :: rtol_T = 0.1D0, atol_T = 2D0
  double precision :: rtol_abun = 0.2D0, atol_abun = 1D-12
  logical flag_converged
  integer n_cell_converged
  real converged_cell_percentage_stop
  character(len=128) iter_files_dir
  character(len=128) :: notes1 = ''
  character(len=128) :: notes2 = ''
  character(len=128) :: notes3 = ''
  character(len=128) :: notes4 = ''
  logical :: flag_save_rates = .FALSE.
  logical :: flag_shortcut_ini = .FALSE.
  logical :: redo_couple_every_column = .FALSE.
  logical :: redo_couple_every_cell = .FALSE.
  logical :: iter_cell_outwards = .TRUE.
  logical :: iter_cell_upwards = .TRUE.
  integer :: nSpecies_check_refine = 0
  integer :: count_refine = 0
  integer :: nMax_refine = 2
  double precision :: threshold_ratio_refine = 10D0
  character(len=128) filename_list_check_refine
end type disk_iteration_params


type :: type_simple_integer_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: vals
end type type_simple_integer_list


type :: disk_analyse_params
  logical :: do_analyse = .false.
  integer ana_i_incr
  character(len=128) analyse_points_inp_dir, analyse_out_dir
  character(len=128) file_list_analyse_points, file_list_analyse_species, &
    file_analyse_res_ele, file_analyse_res_contri
  type(type_cell_rz_phy_basic) chempar
end type disk_analyse_params


type :: disk_iteration_storage
  double precision, dimension(:), allocatable :: T_s
  double precision, dimension(:,:), allocatable :: abundances
end type disk_iteration_storage


type :: book_keeping
  integer fU
  character(len=128) dir, filename_log
end type book_keeping


type(book_keeping) a_book_keeping

type(a__disk) a_disk

type(disk_iteration_params) a_disk_iter_params
type(disk_iteration_storage) a_disk_iter_storage

type(phy_chem_rad_disk_params) disk_params_ini

type(type_cell_rz_phy_basic) cell_params_ini

type(disk_analyse_params) a_disk_ana_params
type(type_simple_integer_list) :: ana_ptlist, ana_splist


integer, dimension(:), allocatable :: calculating_cells_list
integer n_calculating_cells, n_calculating_cells_max

character(len=128) :: filename_save_results
integer fU_save_results

double precision, parameter, private :: xray_energy_kev = 1D0

integer, dimension(:), allocatable, private :: idx_Species_check_refine
double precision, dimension(:), allocatable, private :: thr_Species_check_refine

character(len=256) str_display

integer, parameter :: len_item=14

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
  call save_post_config_params
  !
  ! Now start the major big loop.
  !
  a_disk_iter_params%count_refine = 0
  !
  do ii = 1, a_disk_iter_params%n_iter
    !
    a_disk_iter_params%n_iter_used = ii
    !
    call disk_save_results_pre
    !
    ! Calculate layer by layer.
    ! Start from the surface layer.
    n_calculating_cells = surf_cells%nlen
    calculating_cells_list(1:surf_cells%nlen) = surf_cells%idx
    i_count = 0 ! Counter for cells
    l_count = 0 ! Counter for layers
    do
      l_count = l_count + 1
      !
      do i=1, n_calculating_cells
        i_count = i_count + 1
        i0 = calculating_cells_list(i)
        !
        write(*, '(3(A, I5, A, I5, ",", 2X), (A, I4, ","), 2X, A, 4F8.3)') &
          "Iter:", a_disk_iter_params%n_iter_used, "/", a_disk_iter_params%n_iter, &
          "Cell:", i_count, '/', cell_leaves%nlen, &
          "cell:", i, '/', n_calculating_cells, &
          "Layer:", l_count, &
          'rz:', &
          cell_leaves%list(i0)%p%par%rmin, &
          cell_leaves%list(i0)%p%par%rmax, &
          cell_leaves%list(i0)%p%par%zmin, &
          cell_leaves%list(i0)%p%par%zmax
        write(*, '(A, ES10.3, ",", 4X, 2A, ",", 2X, 2A)') &
          'n_gas: ', cell_leaves%list(i0)%p%par%n_gas, &
          'exe: ', trim(a_disk%params%filename_exe), &
          'dir: ', trim(a_disk_iter_params%iter_files_dir)
        !
        call calc_this_cell(i0)
        !
        call check_convergency_cell(i0)
        !
        write(*, '(12X, 10A10)') chem_idx_some_spe%names(1:10)
        write(*, '(A, 2X, 10ES10.3, L3/)') &
          'Abundances:', cell_leaves%list(i0)%p%abundances(chem_idx_some_spe%idx(1:10)), &
          cell_leaves%list(i0)%p%converged
        !
        a_disk_iter_storage%abundances(:, i0) = &
          cell_leaves%list(i0)%p%abundances(chem_idx_some_spe%idx)
        !
        call disk_save_results_write(fU_save_results, cell_leaves%list(i0)%p)
        flush(fU_save_results)
      end do
      !
      call update_calculating_cells_list
      !
      if (n_calculating_cells .eq. 0) then
        exit
      end if
    end do
    !
    write(str_display, '("! ", A, I4, A)') "Iteration ", ii, " finished."
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_display, a_book_keeping%fU)
    !
    ! At this point all the layers have been walked through.
    call check_convergency_whole_disk
    !
    write(fU_save_results, '(A, L)') '! flag_converged = ', a_disk_iter_params%flag_converged
    write(fU_save_results, '(A)') '! Finish saving ' // trim(filename_save_results)
    write(fU_save_results, '(A)') '! at ' // trim(a_date_time%date_time_str())
    flush(fU_save_results)
    close(fU_save_results)
    !
    if (a_disk_iter_params%flag_converged) then
      exit
    else if (a_disk_iter_params%count_refine .gt. a_disk_iter_params%nMax_refine) then
      write(str_display, &
        '("! Will not refine any more. count_refine: ", I4, " > ", I4)') &
        a_disk_iter_params%count_refine+1, a_disk_iter_params%nMax_refine
      call display_string_both(str_display, a_book_keeping%fU)
      exit
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
          a_disk_iter_params%ncell_refine, cell_leaves%nlen
        !
        call remake_index
        !
        call load_ana_points_list ! Reload, actually
        !
        if (allocated(a_disk_iter_storage%T_s)) then
          deallocate(a_disk_iter_storage%T_s, a_disk_iter_storage%abundances)
        end if
        allocate(a_disk_iter_storage%T_s(cell_leaves%nlen), &
                 a_disk_iter_storage%abundances(chem_idx_some_spe%nItem, &
                                                cell_leaves%nlen))
        do i=1, cell_leaves%nlen
          a_disk_iter_storage%T_s(i) = cell_leaves%list(i)%p%par%Tgas
          a_disk_iter_storage%abundances(:,i) = cell_leaves%list(i)%p%abundances(chem_idx_some_spe%idx)
        end do
        !
        if (allocated(calculating_cells_list)) then
          deallocate(calculating_cells_list)
        end if
        n_calculating_cells_max = cell_leaves%nlen
        allocate(calculating_cells_list(n_calculating_cells_max))
        !
        write(str_display, '("!", A, 2X, I5)') 'New number of cells (leaf):', cell_leaves%nlen
        call display_string_both(str_display, a_book_keeping%fU)
        write(str_display, '("!", A, 2X, I5)') 'New number of cells (total):', root%nOffspring
        call display_string_both(str_display, a_book_keeping%fU)
      else
        write(str_display, '("! ", A)') "No further refinement needed."
        call display_string_both(str_display, a_book_keeping%fU)
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
  write(str_display, '("!Final number of cells =", I4)') cell_leaves%nlen
  call display_string_both(str_display, a_book_keeping%fU)
  !
  ! call disk_iteration_postproc
  !
end subroutine disk_iteration


subroutine disk_iteration_prepare
  integer i
  !
  call make_grid
  n_calculating_cells_max = cell_leaves%nlen
  allocate(calculating_cells_list(n_calculating_cells_max))
  !
  call chem_read_reactions()
  call chem_load_reactions()
  call chem_parse_reactions()
  call chem_get_dupli_reactions()
  call chem_get_idx_for_special_species()
  call load_species_enthalpies
  call get_reaction_heat
  !
  write(str_display, '("!", A, 2X, I5)') 'Number of cells (leaf):', cell_leaves%nlen
  call display_string_both(str_display, a_book_keeping%fU)
  write(str_display, '("!", A, 2X, I5)') 'Number of cells (total):', root%nOffspring
  call display_string_both(str_display, a_book_keeping%fU)
  write(str_display, '("!", A, 2X, I5)') 'Number of reactions:', chem_net%nReactions
  call display_string_both(str_display, a_book_keeping%fU)
  write(str_display, '("!", A, 2X, I5)') 'Number of species ', chem_species%nSpecies
  call display_string_both(str_display, a_book_keeping%fU)
  !
  chem_solver_params%fU_log = a_book_keeping%fU
  !
  call load_refine_check_species
  !
  call disk_set_disk_params
  call disk_set_gridcell_params
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  if (.NOT. allocated(a_disk_iter_storage%T_s)) then
    allocate(a_disk_iter_storage%T_s(cell_leaves%nlen), &
             a_disk_iter_storage%abundances(chem_idx_some_spe%nItem, &
                                            cell_leaves%nlen))
  end if
  !
  call chem_load_initial_abundances
  !
  do i=1, cell_leaves%nlen
    cell_leaves%list(i)%p%abundances = chem_solver_storage%y(1:chem_species%nSpecies)
    a_disk_iter_storage%T_s(i) = cell_leaves%list(i)%p%par%Tgas
    a_disk_iter_storage%abundances(:, i) = chem_solver_storage%y(chem_idx_some_spe%idx)
  end do
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
end subroutine disk_iteration_prepare


subroutine update_params_above(i0)
  use load_Visser_CO_selfshielding
  integer, intent(in) :: i0
  integer j, j0
  associate(p  => cell_leaves%list(i0)%p, &
            dz => cell_leaves%list(i0)%p%par%dz * phy_AU2cm)
    p%col_den   = p%abundances(chem_idx_some_spe%idx) * p%par%n_gas * dz
    p%par%dNcol = p%par%n_gas * dz
    p%col_den_acc = 0D0
    p%par%Ncol = 0D0
    if (cell_leaves%list(i0)%p%above%n .gt. 0) then
      do j=1, cell_leaves%list(i0)%p%above%n
        j0 = cell_leaves%list(i0)%p%above%idx(j)
        p%col_den_acc = p%col_den_acc + &
                        (cell_leaves%list(j0)%p%col_den_acc + &
                         cell_leaves%list(j0)%p%col_den) * p%above%fra(j)
        p%par%Ncol    = p%par%Ncol + &
                        (cell_leaves%list(j0)%p%par%Ncol + &
                         cell_leaves%list(j0)%p%par%dNcol) * p%above%fra(j)
      end do
    end if
  end associate
  associate(p        => cell_leaves%list(i0)%p, &
            dz       => cell_leaves%list(i0)%p%par%dz * phy_AU2cm, &
            Ncol_H2  => cell_leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiH2), &
            dcol_H2  => cell_leaves%list(i0)%p%col_den(chem_idx_some_spe%iiH2), &
            Ncol_H   => cell_leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiHI), &
            dcol_H   => cell_leaves%list(i0)%p%col_den(chem_idx_some_spe%iiHI), &
            Ncol_H2O => cell_leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiH2O), &
            dcol_H2O => cell_leaves%list(i0)%p%col_den(chem_idx_some_spe%iiH2O), &
            Ncol_OH  => cell_leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiOH), &
            dcol_OH  => cell_leaves%list(i0)%p%col_den(chem_idx_some_spe%iiOH), &
            Ncol_CO  => cell_leaves%list(i0)%p%col_den_acc(chem_idx_some_spe%iiCO), &
            dcol_CO  => cell_leaves%list(i0)%p%col_den(chem_idx_some_spe%iiCO))
    ! Kwok eq 10.20
    p%par%Av = 1.086D0 * p%par%ratioDust2HnucNum * &
      (phy_Pi * p%par%GrainRadius_CGS**2) * 2D0 * &
      p%par%Ncol
      !(p%par%Ncol + p%par%dNcol * 0.5D0)
    p%par%f_selfshielding_H2  = &
      min(1D0, get_H2_self_shielding(Ncol_H2, p%par%velo_width_turb))
      !min(1D0, get_H2_self_shielding(Ncol_H2 + dcol_H2*0.5D0, p%par%velo_width_turb))
      !min(1D0, ((Ncol_H2 + dcol_H2*0.5D0)/1D14)**(-0.75D0)) ! Tielens 2005, equation 8.39
    p%par%f_selfshielding_H2O = &
      min(1D0, exp(-(Ncol_H2O * const_LyAlpha_cross_H2O))) !* &
      !tau2beta(dcol_H2O * const_LyAlpha_cross_H2O)
    p%par%f_selfshielding_OH  = &
      min(1D0, exp(-(Ncol_OH * const_LyAlpha_cross_OH))) !* &
      !tau2beta(dcol_OH * const_LyAlpha_cross_OH)
    p%par%f_selfshielding_CO = min(1D0, max(0D0, get_12CO_shielding(Ncol_H2, Ncol_CO)))
  end associate
end subroutine update_params_above


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
  if (maxval(abs(cell_leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                 - a_disk_iter_storage%abundances(:, i0)) &
             - (a_disk_iter_params%atol_abun + &
                a_disk_iter_params%rtol_abun * &
                abs(cell_leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                  + a_disk_iter_storage%abundances(:, i0))) &
            ) .le. 0D0) then
    cell_leaves%list(i0)%p%converged = .true.
  else
    cell_leaves%list(i0)%p%converged = .false.
  end if
end subroutine check_convergency_cell


subroutine check_convergency_whole_disk
  integer i
  a_disk_iter_params%n_cell_converged = 0
  do i=1, cell_leaves%nlen
    if (cell_leaves%list(i)%p%converged) then
      a_disk_iter_params%n_cell_converged = a_disk_iter_params%n_cell_converged + 1
    end if
  end do
  a_disk_iter_params%flag_converged = &
    a_disk_iter_params%n_cell_converged .ge. &
    int(a_disk_iter_params%converged_cell_percentage_stop * real(cell_leaves%nlen))
  write(str_display, '("! Iter", I4, " Number of cells converged: ", I6, "/", I6)') &
    a_disk_iter_params%n_iter_used, a_disk_iter_params%n_cell_converged, cell_leaves%nlen
  call display_string_both(str_display, a_book_keeping%fU)
end subroutine check_convergency_whole_disk


subroutine update_calculating_cells_list
  integer, dimension(:), allocatable :: list_tmp
  integer i, i0, itmp, j, k, n
  logical flag_notyet
  allocate(list_tmp(n_calculating_cells_max))
  n = 0
  do i=1, n_calculating_cells
    i0 = calculating_cells_list(i)
    do j=1, cell_leaves%list(i0)%p%below%n
      itmp = cell_leaves%list(i0)%p%below%idx(j)
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
    calculating_cells_list(1:n) = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine update_calculating_cells_list


subroutine calc_this_cell(id)
  integer, intent(in) :: id
  integer i, j, i0, ntmp
  double precision tmp
  logical found_neighbor
  !
  cell_leaves%list(id)%p%iIter = a_disk_iter_params%n_iter_used
  !
  do j=1, min(4,a_disk_iter_params%nlocal_iter)
    !
    write(*, '("Local iter: ", I4, " of ", I4)') j, min(4,a_disk_iter_params%nlocal_iter)
    !
    ! Set the initial condition for chemical evolution
    if (a_disk_iter_params%flag_shortcut_ini) then
      if (a_disk_iter_params%n_iter_used .eq. 1) then
        found_neighbor = .false.
        do i=1, cell_leaves%list(id)%p%around%n
          i0 = cell_leaves%list(id)%p%around%idx(i)
          if (cell_leaves%list(i0)%p%iIter .gt. cell_leaves%list(id)%p%iIter) then
            chem_solver_storage%y(1:chem_species%nSpecies) = cell_leaves%list(i0)%p%abundances
            found_neighbor = .true.
            exit
          end if
        end do
        if (.not. found_neighbor) then
          chem_solver_storage%y(1:chem_species%nSpecies) = chem_solver_storage%y0(1:chem_species%nSpecies)
        end if
      else
        chem_solver_storage%y(1:chem_species%nSpecies) = cell_leaves%list(id)%p%abundances
      end if
      if (chem_solver_params%neutralize) then
        tmp = sum(chem_solver_storage%y(1:chem_species%nSpecies) * &
                          dble(chem_species%elements(1,:)))
        if (abs(tmp) .ge. 1D-2*chem_solver_storage%y(chem_idx_some_spe%i_E)) then
          chem_solver_storage%y(1:chem_species%nSpecies) = chem_solver_storage%y0(1:chem_species%nSpecies)
        else
          chem_solver_storage%y(chem_idx_some_spe%i_E) = &
            chem_solver_storage%y(chem_idx_some_spe%i_E) + tmp
          if (chem_solver_storage%y(chem_idx_some_spe%i_E) .lt. 0D0) then
            ! When it is not possible to neutralize the composition by artificially
            ! changing the electron abundance, then use the general initial abundances,
            ! which should be absolutely neutral.
            chem_solver_storage%y(1:chem_species%nSpecies) = &
              chem_solver_storage%y0(1:chem_species%nSpecies)
            write(str_display, '("! Cannot neutralize: X(E-) = ", ES12.4)') &
              chem_solver_storage%y(chem_idx_some_spe%i_E)
            call display_string_both(str_display, a_book_keeping%fU)
            write(str_display, '("! Use y0 as initial abundance.")')
            call display_string_both(str_display, a_book_keeping%fU)
            write(str_display, '("! x, y = ", 2ES10.2, " iIter = ", I4)') &
              cell_leaves%list(id)%p%xmin, cell_leaves%list(id)%p%ymin, cell_leaves%list(id)%p%iIter
            call display_string_both(str_display, a_book_keeping%fU)
          end if
        end if
      end if
    else
      chem_solver_storage%y(1:chem_species%nSpecies) = chem_solver_storage%y0(1:chem_species%nSpecies)
    end if
    !
    ! Always use the temperature of the above cell as init
    if ((a_disk_iter_params%n_iter_used .eq. 1) .and. (j .eq. 1)) then
      if (cell_leaves%list(id)%p%above%n .gt. 0) then
        cell_leaves%list(id)%p%par%Tgas = 0D0
        ntmp = 0
        do i=1, cell_leaves%list(id)%p%above%n
          i0 = cell_leaves%list(id)%p%above%idx(i)
          if (cell_leaves%list(i0)%p%iIter .lt. cell_leaves%list(id)%p%iIter) then
            cycle
          end if 
          ntmp = ntmp + 1
          cell_leaves%list(id)%p%par%Tgas = cell_leaves%list(id)%p%par%Tgas + &
                                            cell_leaves%list(i0)%p%par%Tgas
        end do
        if (ntmp .gt. 0) then
          cell_leaves%list(id)%p%par%Tgas = cell_leaves%list(id)%p%par%Tgas / &
                                            dble(ntmp)
        else
          cell_leaves%list(id)%p%par%Tgas = cell_leaves%list(id)%p%par%Tdust
        end if
      else
        cell_leaves%list(id)%p%par%Tgas = cell_leaves%list(id)%p%par%Tdust
      end if
    end if
    !
    chem_solver_storage%y(chem_species%nSpecies+1) = cell_leaves%list(id)%p%par%Tgas
    !
    write(*, '(4X, A, F12.3/)') 'Tgas_old: ', cell_leaves%list(id)%p%par%Tgas
    !
    call update_params_above(id)
    !
    call set_chemistry_params_from_cell(id)
    call chem_cal_rates
    !
    call set_heatingcooling_params_from_cell(id)
    !
    chem_solver_params%flag_chem_evol_save = .false.
    !
    call chem_set_solver_flags_alt(j)
    if (j .eq. 1) then
      chem_solver_params%evolT = .true.
    else if (abs(heating_cooling_rates%hc_net_rate) .le. &
            1D-2 * min(max_heating_rate(), max_cooling_rate())) then
      chem_solver_params%evolT = .false.
    else
      chem_solver_storage%y(chem_species%nSpecies+1) = &
        (0.5D0 + dble(j)*0.1D0) * (cell_leaves%list(id)%p%par%Tgas + cell_leaves%list(id)%p%par%Tdust)
      chem_solver_params%evolT = .true.
    end if
    !
    call chem_evol_solve
    !
    cell_leaves%list(id)%p%abundances = chem_solver_storage%y(1:chem_species%nSpecies)
    cell_leaves%list(id)%p%par%Tgas = chem_solver_storage%y(chem_species%nSpecies+1)
    cell_leaves%list(id)%p%quality = chem_solver_params%quality
    cell_leaves%list(id)%p%par%t_final = chem_solver_storage%touts(chem_solver_params%n_record_real)
    !
    !if (minval(chem_solver_storage%y(1:chem_species%nSpecies)) .lt. 0D0) then
    !  if (cell_leaves%list(id)%p%quality .eq. 0) then
    !    cell_leaves%list(id)%p%quality = -1
    !  else
    !    cell_leaves%list(id)%p%quality = -cell_leaves%list(id)%p%quality
    !  end if
    !end if
    !
    call update_params_above(id)
    !
    write(*, '(4X, A, F12.3)') 'Tgas_new: ', cell_leaves%list(id)%p%par%Tgas
    !
    if (cell_leaves%list(id)%p%quality .eq. 0) then
      exit
    end if
    !
  end do
  !
  a_disk_iter_storage%T_s(id) = cell_leaves%list(id)%p%par%Tgas
  !
  if (.not. chem_solver_params%evolT) then
    call chem_cal_rates
    call realtime_heating_cooling_rate(tmp, chem_solver_params%NEQ, chem_solver_storage%y)
  end if
  cell_leaves%list(id)%p%h_c_rates = heating_cooling_rates
  !
  !if (a_disk_iter_params%flag_save_rates) then
  !  call save_chem_rates(id)
  !end if
  !
  if (a_disk_ana_params%do_analyse) then
    if ((cell_leaves%list(id)%p%quality .gt. 0) .or. &
        ((a_disk_iter_params%n_iter_used .gt. 1) .and. &
         (is_in_list_int(id, ana_ptlist%nlen, ana_ptlist%vals) .or. &
          need_to_refine(cell_leaves%list(id)%p)))) then
      call chem_analyse(id)
    end if
  end if
end subroutine calc_this_cell


subroutine disk_save_results_pre
  if (.NOT. getFileUnit(fU_save_results)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  write(filename_save_results, '("iter_", I4.4, ".dat")') a_disk_iter_params%n_iter_used
  filename_save_results = combine_dir_filename(a_disk_iter_params%iter_files_dir, filename_save_results)
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
    str_pad_to_len('t_final', len_item) // &
    str_pad_to_len('rmin',    len_item) // &
    str_pad_to_len('rmax',    len_item) // &
    str_pad_to_len('zmin',    len_item) // &
    str_pad_to_len('zmax',    len_item) // &
    str_pad_to_len('Tgas',    len_item) // &
    str_pad_to_len('Tdust',   len_item) // &
    str_pad_to_len('n_gas',   len_item) // &
    str_pad_to_len('Av',      len_item) // &
    str_pad_to_len('UV_G0',   len_item) // &
    str_pad_to_len('LyA_G0',  len_item) // &
    str_pad_to_len('LyANF0',  len_item) // &
    str_pad_to_len('XRay0',   len_item) // &
    str_pad_to_len('Ncol',    len_item) // &
    str_pad_to_len('dNcol',   len_item) // &
    str_pad_to_len('N_H2',    len_item) // &
    str_pad_to_len('N_H2O',   len_item) // &
    str_pad_to_len('N_OH',    len_item) // &
    str_pad_to_len('N_CO',    len_item) // &
    str_pad_to_len('dN_H2',   len_item) // &
    str_pad_to_len('dN_H2O',  len_item) // &
    str_pad_to_len('dN_OH',   len_item) // &
    str_pad_to_len('dN_CO',   len_item) // &
    str_pad_to_len('f_H2',    len_item) // &
    str_pad_to_len('f_H2O',   len_item) // &
    str_pad_to_len('f_OH',    len_item) // &
    str_pad_to_len('f_CO',    len_item) // &
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
  write(fU, '(7I5, 54ES14.6E2' // trim(fmt_str)) &
  converged                                              , &
  c%quality                                              , &
  c%around%n                                             , &
  c%above%n                                              , &
  c%below%n                                              , &
  c%inner%n                                              , &
  c%outer%n                                              , &
  c%par%t_final                                          , &
  c%par%rmin                                             , &
  c%par%rmax                                             , &
  c%par%zmin                                             , &
  c%par%zmax                                             , &
  c%par%Tgas                                             , &
  c%par%Tdust                                            , &
  c%par%n_gas                                            , &
  c%par%Av                                               , &
  c%par%UV_G0_factor                                     , &
  c%par%LymanAlpha_G0_factor                             , &
  c%par%LymanAlpha_number_flux_0                         , &
  c%par%Xray_flux_0                                      , &
  c%par%Ncol                                             , &
  c%par%dNcol                                            , &
  c%col_den_acc(chem_idx_some_spe%iiH2)                  , &
  c%col_den_acc(chem_idx_some_spe%iiH2O)                 , &
  c%col_den_acc(chem_idx_some_spe%iiOH)                  , &
  c%col_den_acc(chem_idx_some_spe%iiCO)                  , &
  c%col_den(chem_idx_some_spe%iiH2)                      , &
  c%col_den(chem_idx_some_spe%iiH2O)                     , &
  c%col_den(chem_idx_some_spe%iiOH)                      , &
  c%col_den(chem_idx_some_spe%iiCO)                      , &
  c%par%f_selfshielding_H2                               , &
  c%par%f_selfshielding_H2O                              , &
  c%par%f_selfshielding_OH                               , &
  c%par%f_selfshielding_CO                               , &
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
  a_disk%params%disk_mass_in_Msun = 0D0
  do i=1, cell_leaves%nlen
    associate(p => cell_leaves%list(i)%p%par)
      a_disk%params%disk_mass_in_Msun = &
        a_disk%params%disk_mass_in_Msun + &
        p%n_gas * p%MeanMolWeight * (phy_2Pi * p%rcen * p%dr * p%dz)
    end associate
  end do
  a_disk%params%disk_mass_in_Msun = a_disk%params%disk_mass_in_Msun * &
    phy_AU2cm**3 * phy_mProton_CGS / phy_Msun_CGS
end subroutine disk_calc_disk_mass


subroutine set_heatingcooling_params_from_cell(id)
  integer id
  heating_cooling_params%type_cell_rz_phy_basic = cell_leaves%list(id)%p%par
  heating_cooling_params%Neufeld_dv_dz = cell_leaves%list(id)%p%par%velo_gradient * 1D-5 ! cm s-1 to km s-1
  heating_cooling_params%Neufeld_G     = 1D0
  heating_cooling_params%X_H2    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)
  heating_cooling_params%X_HI    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI)
  heating_cooling_params%X_CI    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_CI)
  heating_cooling_params%X_Cplus = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_Cplus)
  heating_cooling_params%X_OI    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_OI)
  heating_cooling_params%X_CO    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_CO)
  heating_cooling_params%X_H2O   = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2O)
  heating_cooling_params%X_OH    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_OH)
  heating_cooling_params%X_E     = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_E)
  heating_cooling_params%X_Hplus = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_Hplus)
  heating_cooling_params%X_gH    = cell_leaves%list(id)%p%abundances(chem_idx_some_spe%i_gH)
  if (chem_solver_params%H2_form_use_moeq) then
    heating_cooling_params%R_H2_form_rate = &
      heating_cooling_params%R_H2_form_rate_coeff * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%X_HI * &
      heating_cooling_params%n_gas
  else
    heating_cooling_params%R_H2_form_rate = &
      heating_cooling_params%R_H2_form_rate_coeff * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%n_gas
  end if
  !
  cell_leaves%list(id)%p%par%R_H2_form_rate = heating_cooling_params%R_H2_form_rate
end subroutine set_heatingcooling_params_from_cell


subroutine set_chemistry_params_from_cell(id)
  integer id
  chem_params => cell_leaves%list(id)%p%par
end subroutine set_chemistry_params_from_cell


subroutine disk_set_a_cell_params(c, cell_params_copy)
  type(type_cell), target :: c
  type(type_cell_rz_phy_basic), intent(in) :: cell_params_copy
  if (.not. associated(c%par)) then
    allocate(c%par)
  end if
  if (.not. allocated(c%h_c_rates)) then
    allocate(c%h_c_rates)
  end if
  !
  if (.not. allocated(c%abundances)) then
    allocate(c%abundances(chem_species%nSpecies), &
             c%col_den(chem_idx_some_spe%nItem), &
             c%col_den_acc(chem_idx_some_spe%nItem))
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
  c%par%n_gas  = c%val(1) ! Already set
  if (grid_config%use_data_file_input) then
    c%par%Tgas    = c%val(2)
    c%par%Tdust   = c%val(2)
  else
    c%par%Tgas    = 400D0 / (1D0 + c%par%rcen) * (1D0 + c%par%zcen)
    c%par%Tdust   = c%par%Tgas
  end if
  !
  c%par%ratioDust2HnucNum = &
      c%par%ratioDust2GasMass * (phy_mProton_CGS * c%par%MeanMolWeight) &
      / (4.0D0*phy_Pi/3.0D0 * (c%par%GrainRadius_CGS)**3 * &
         c%par%GrainMaterialDensity_CGS)
  c%par%dust_depletion = c%par%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
  c%par%n_dust = c%par%n_gas * c%par%ratioDust2HnucNum
  !
  c%par%UV_G0_factor = c%par%UV_G0_factor_background + &
    a_disk%params%UV_cont_phlumi_star_surface &
       / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
       / phy_Habing_photon_flux_CGS &
       * a_disk%params%geometric_factor_UV
  c%par%LymanAlpha_number_flux_0 = &
    a_disk%params%Lyman_phlumi_star_surface &
       / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2)
  c%par%LymanAlpha_energy_flux_0 = c%par%LymanAlpha_number_flux_0 * phy_LyAlpha_energy_CGS
  c%par%LymanAlpha_G0_factor = c%par%LymanAlpha_energy_flux_0 / phy_Habing_energy_flux_CGS
  c%par%Xray_flux_0 = &
    a_disk%params%Xray_phlumi_star_surface &
       / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
       * a_disk%params%geometric_factor_Xray
  associate( &
          G     => phy_GravitationConst_CGS, &
          M     => a_disk%params%star_mass_in_Msun * phy_Msun_CGS, &
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
  do i=1, cell_leaves%nlen
    call disk_set_a_cell_params(cell_leaves%list(i)%p, cell_params_ini)
  end do
end subroutine disk_set_gridcell_params


subroutine disk_set_disk_params
  a_disk%params = disk_params_ini
  associate( &
    Lstar => a_disk%params%star_luminosity_in_Lsun * phy_Lsun_CGS, &
    uv2total => a_disk%params%ratio_uv2total, &
    lyman2uv => a_disk%params%ratio_lyman2uv, &
    xray2total => a_disk%params%ratio_xray2total)
    a_disk%params%UV_cont_phlumi_star_surface = &
      Lstar * uv2total * (1D0 - lyman2uv) / phy_UV_cont_energy_CGS
    a_disk%params%Lyman_phlumi_star_surface = &
      Lstar * uv2total * lyman2uv         / phy_LyAlpha_energy_CGS
    a_disk%params%Xray_phlumi_star_surface  = &
      Lstar * xray2total / (xray_energy_kev*1D3*phy_eV2erg)
    write(str_display, '("!Stellar total luminosity = ", ES12.4, " erg s-1")'), Lstar
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar UV cont luminosity = ", ES12.4, " erg s-1")'), &
      a_disk%params%UV_cont_phlumi_star_surface * phy_UV_cont_energy_CGS
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar UV cont photon count rate = ", ES12.4, " s-1")'), &
      a_disk%params%UV_cont_phlumi_star_surface
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar LyA luminosity = ", ES12.4, " erg s-1")'), &
      a_disk%params%Lyman_phlumi_star_surface * phy_LyAlpha_energy_CGS
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar LyA photon count rate = ", ES12.4, " s-1")'), &
      a_disk%params%Lyman_phlumi_star_surface
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar X-ray luminosity = ", ES12.4, " erg s-1")'), &
      a_disk%params%Xray_phlumi_star_surface * (xray_energy_kev*1D3*phy_eV2erg)
    call display_string_both(str_display, a_book_keeping%fU)
    write(str_display, '("!Stellar X-ray photon count rate = ", ES12.4, " s-1")'), &
      a_disk%params%Xray_phlumi_star_surface
    call display_string_both(str_display, a_book_keeping%fU)
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
  get_local_dv_microturb = 1D5 ! = 1 km s-1
end function get_local_dv_microturb


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
  dir = combine_dir_filename(a_book_keeping%dir, 'rates_log/')
  if (.NOT. dir_exist(dir)) then
    call my_mkdir(dir)
  end if
  call openFileSequentialWrite(fU, combine_dir_filename(dir, filename), 99999)
  !
  heat_cool_log = heating_cooling_params
  write(fU, nml=cell_par_log)
  !
  do k=1, chem_net%nReactions
    write(fU, '(A135, ES16.4E4)') chem_reac_str%list(k), chem_net%rates(k)
  end do
  close(fU)
end subroutine save_chem_rates


subroutine save_post_config_params
  type(phy_chem_rad_disk_params) disk_params_tmp
  namelist /disk_params_log/ disk_params_tmp
  disk_params_tmp = a_disk%params
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
  write(str_display, '("! Species used for checking refinement:")')
  call display_string_both(str_display, a_book_keeping%fU)
  do i=1, a_disk_iter_params%nSpecies_check_refine
    write(str_display, '("! ", A12, ES12.2)') &
      chem_species%names(idx_Species_check_refine(i)), &
      thr_Species_check_refine(i)
    call display_string_both(str_display, a_book_keeping%fU)
  end do
end subroutine load_refine_check_species


subroutine do_refine
  integer i, n_refine
  a_disk_iter_params%ncell_refine = 0
  do i=1, cell_leaves%nlen
    if (need_to_refine(cell_leaves%list(i)%p, n_refine)) then
      a_disk_iter_params%ncell_refine = a_disk_iter_params%ncell_refine + 1
      write(str_display, '("!", I4, A, 4ES12.2, " into ", I4, " parts.")') &
        a_disk_iter_params%ncell_refine, ' Refining ', &
        cell_leaves%list(i)%p%xmin, cell_leaves%list(i)%p%xmax, &
        cell_leaves%list(i)%p%ymin, cell_leaves%list(i)%p%ymax, n_refine
      call display_string_both(str_display, a_book_keeping%fU)
      call refine_this_cell_vertical(cell_leaves%list(i)%p, n_refine)
    end if
  end do
end subroutine do_refine


subroutine remake_index
  call get_number_of_leaves(root)
  cell_leaves%nlen = root%nleaves
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
  if (c%par%dz .le. grid_config%very_small_len) then
    need_to_refine = .false.
    return
  end if
  do i=1, c%above%n
    i0 = c%above%idx(i)
    !if (cell_leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (cell_leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(cell_leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(cell_leaves%list(i0)%p%abundances(i1), c%abundances(i1))
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
    !if (cell_leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (cell_leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(cell_leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(cell_leaves%list(i0)%p%abundances(i1), c%abundances(i1))
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
    do i=1, cell_leaves%nlen
      if ((cell_leaves%list(i)%p%par%rmin .le. r) .and. (cell_leaves%list(i)%p%par%rmax .ge. r) .and. &
          (cell_leaves%list(i)%p%par%zmin .le. z) .and. (cell_leaves%list(i)%p%par%zmax .ge. z)) then
        idx = i
        exit
      end if
    end do
    if (idx .eq. 0) then
      write(*, '("Point (", 2F6.2, ")", A)') r, z, ' not in any cells!'
      cycle
    end if
    chem_solver_storage%y(1:chem_species%nSpecies) = &
        cell_leaves%list(idx)%p%abundances(1:chem_species%nSpecies)
    call chem_elemental_residence
    write(fU1, '("(", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      cell_leaves%list(idx)%p%par%Tgas, cell_leaves%list(idx)%p%par%Tdust, &
      cell_leaves%list(idx)%p%par%n_gas, cell_leaves%list(idx)%p%par%Ncol, &
      cell_leaves%list(idx)%p%par%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chem_solver_storage%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    if (cell_leaves%list(idx)%p%above%n .gt. 0) then
      idx_diff = cell_leaves%list(idx)%p%above%idx(1)
    else if (cell_leaves%list(idx)%p%below%n .gt. 0) then
      idx_diff = cell_leaves%list(idx)%p%below%idx(1)
    else
      idx_diff = idx
    end if
    !
    call set_chemistry_params_from_cell(idx_diff)
    call chem_cal_rates
    call get_contribution_each
    !
    write(fU2, '("This (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      cell_leaves%list(idx)%p%par%Tgas, cell_leaves%list(idx)%p%par%Tdust, &
      cell_leaves%list(idx)%p%par%n_gas, cell_leaves%list(idx)%p%par%Ncol, &
      cell_leaves%list(idx)%p%par%Av
    write(fU2, '("Diff (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') &
      cell_leaves%list(idx_diff)%p%par%rcen, &
      cell_leaves%list(idx_diff)%p%par%zcen, &
      cell_leaves%list(idx_diff)%p%par%Tgas,  cell_leaves%list(idx_diff)%p%par%Tdust, &
      cell_leaves%list(idx_diff)%p%par%n_gas, cell_leaves%list(idx_diff)%p%par%Ncol, &
      cell_leaves%list(idx_diff)%p%par%Av
    do i=1, chem_species%nSpecies
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(A12, ": ", ES12.2, " Diff: ", ES12.2, " Rate: ", ES12.2)') chem_species%names(i), &
        chem_solver_storage%y(i), cell_leaves%list(idx_diff)%p%abundances(i), &
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
  allocate(list_tmp(cell_leaves%nlen))
  n = 0
  do
    read(fU, '(2F6.2)', iostat=ios) r, z
    if (ios .ne. 0) then
      exit
    end if
    do i=1, cell_leaves%nlen
      if ((cell_leaves%list(i)%p%par%rmin .le. r) .and. (cell_leaves%list(i)%p%par%rmax .ge. r) .and. &
          (cell_leaves%list(i)%p%par%zmin .le. z) .and. (cell_leaves%list(i)%p%par%zmax .ge. z)) then
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
  if (IsWordChar(chem_solver_params%filename_species_enthalpy(1:1))) then
    if (.NOT. getFileUnit(fU)) then
      write (*,*) 'In subroutine ImportSpeciesEnthalpy:'
      write (*,*) 'Cannot allocate an output file unit!'
      stop
    end if
    write (FMTstr, FMT= '("(", "A", I2, ", F", I1, ".0)")') &
      const_len_species_name, 9
    CALL openFileSequentialRead(fU, &
        combine_dir_filename(chem_solver_params%chem_files_dir, &
        chem_solver_params%filename_species_enthalpy), 999999)
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
    write(str_display, '("! ", I5, A, I5, A)') &
      i1, ' of ', chem_species%nSpecies, ' species have enthalpy.'
    call display_string_both(str_display, a_book_keeping%fU)
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
        id, cell_leaves%list(id)%p%xmin, cell_leaves%list(id)%p%ymin, &
        cell_leaves%list(id)%p%iIter
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), 999999)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chem_solver_params%n_record
    write(fU3, FMTstryHistory) chem_solver_storage%touts(i), chem_solver_storage%record(:, i)
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
    a_disk_ana_params%ana_i_incr = 1+chem_solver_params%n_record/20
  end if
  !
  write(fU1, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, cell_leaves%list(id)%p%iIter, &
    cell_leaves%list(id)%p%xmin, cell_leaves%list(id)%p%xmax, &
    cell_leaves%list(id)%p%ymin, cell_leaves%list(id)%p%ymax
  write(fU2, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, cell_leaves%list(id)%p%iIter, &
    cell_leaves%list(id)%p%xmin, cell_leaves%list(id)%p%xmax, &
    cell_leaves%list(id)%p%ymin, cell_leaves%list(id)%p%ymax
  do k=1, chem_solver_params%n_record, a_disk_ana_params%ana_i_incr
    !
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chem_solver_storage%record(:, k) - chem_solver_storage%record(:, k-1))) / &
        (chem_solver_storage%record(:, k) + chem_solver_storage%record(:, k-1) + 1D-15))
      dt_t = (chem_solver_storage%touts(k) - chem_solver_storage%touts(k-1)) / &
             (chem_solver_storage%touts(k) + chem_solver_storage%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    write(fU1, '("time = ", ES14.4)') chem_solver_storage%touts(k)
    !
    chem_solver_storage%y = chem_solver_storage%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chem_solver_storage%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chem_solver_storage%touts(k)
    if (ana_splist%nlen .le. 0) then
      cycle
    end if
    call get_contribution_each
    do i=1, chem_species%nSpecies
      if (.not. is_in_list_int(i, ana_splist%nlen, ana_splist%vals)) then
        cycle
      end if
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chem_solver_storage%y(i)
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
    !ch%stickCoeffH = 1D0
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
  end associate
  !
  call chem_cal_rates
  call chem_set_solver_flags
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
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 1
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 1
  do i=1, chem_solver_params%n_record
    write(fU3, FMTstryHistory) chem_solver_storage%touts(i), chem_solver_storage%record(:, i)
  end do
  close(fU3)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = chem_solver_params%n_record / 4
  end if
  do k=1, chem_solver_params%n_record, a_disk_ana_params%ana_i_incr
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
    chem_solver_storage%y = chem_solver_storage%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chem_solver_storage%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chem_solver_storage%y(:) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chem_solver_storage%y(:) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chem_solver_storage%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    call get_contribution_each
    !
    write(fU2, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chem_solver_storage%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    do i=1, chem_species%nSpecies
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chem_solver_storage%y(i)
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
    close(fU1)
    close(fU2)
  end do
end subroutine a_test_case



subroutine b_test_case
  integer i, j, k, fU, fU3
  type(type_cell_rz_phy_basic), pointer :: ch => null()
  double precision Tmin, Tmax, dT, ratio
  double precision n_gas_min, n_gas_max, dn
  double precision h_c_net_rate, tmp
  character(len=128) filename, fname_pre, header
  character(len=32) FMTstryHistory
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
      chem_solver_params%evolT = .true.
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
      chem_solver_storage%y(1:chem_species%nSpecies) = chem_solver_storage%y0(1:chem_species%nSpecies)
      chem_solver_storage%y(chem_species%nSpecies+1) = ch%Tgas
      !
      call chem_cal_rates
      write(*,'(2ES12.2/)') chem_params%f_selfshielding_H2, chem_params%Av
      call chem_set_solver_flags_alt(1)
      !
      heating_cooling_params%type_cell_rz_phy_basic = ch
      !
      heating_cooling_params%Neufeld_dv_dz = 10D0/phy_AU2cm
      heating_cooling_params%Neufeld_G     = 1D0
      !
      heating_cooling_params%X_H2    = chem_solver_storage%y(chem_idx_some_spe%i_H2)
      heating_cooling_params%X_HI    = chem_solver_storage%y(chem_idx_some_spe%i_HI)
      heating_cooling_params%X_CI    = chem_solver_storage%y(chem_idx_some_spe%i_CI)
      heating_cooling_params%X_Cplus = chem_solver_storage%y(chem_idx_some_spe%i_Cplus)
      heating_cooling_params%X_OI    = chem_solver_storage%y(chem_idx_some_spe%i_OI)
      heating_cooling_params%X_CO    = chem_solver_storage%y(chem_idx_some_spe%i_CO)
      heating_cooling_params%X_H2O   = chem_solver_storage%y(chem_idx_some_spe%i_H2O)
      heating_cooling_params%X_OH    = chem_solver_storage%y(chem_idx_some_spe%i_OH)
      heating_cooling_params%X_E     = chem_solver_storage%y(chem_idx_some_spe%i_E)
      heating_cooling_params%X_Hplus = chem_solver_storage%y(chem_idx_some_spe%i_Hplus)
      heating_cooling_params%X_gH    = chem_solver_storage%y(chem_idx_some_spe%i_gH)
      !
      if (chem_solver_params%H2_form_use_moeq) then
        heating_cooling_params%R_H2_form_rate = &
          heating_cooling_params%R_H2_form_rate_coeff * &
          heating_cooling_params%X_gH * &
          heating_cooling_params%X_HI * &
          heating_cooling_params%n_gas
      else
        heating_cooling_params%R_H2_form_rate = &
          heating_cooling_params%R_H2_form_rate_coeff * &
          heating_cooling_params%X_gH * &
          heating_cooling_params%X_gH * &
          heating_cooling_params%n_gas
      end if
      ch%R_H2_form_rate = heating_cooling_params%R_H2_form_rate
      !
      !call realtime_heating_cooling_rate(tmp, chem_solver_params%NEQ, chem_solver_storage%y)
      !write(*,'(2ES16.6/)') tmp, heating_minus_cooling()
      !call disp_h_c_rates
      call chem_evol_solve
      !
      c%abundances  = chem_solver_storage%y(1:chem_species%nSpecies)
      c%col_den     = c%abundances(chem_idx_some_spe%idx) * c%par%dNcol
      c%col_den_acc = c%abundances(chem_idx_some_spe%idx) * c%par%Ncol
      !
      hc_Tgas = ch%Tgas
      hc_Tdust = ch%Tdust
      h_c_net_rate = heating_minus_cooling()
      !
      c%h_c_rates = heating_cooling_rates
      c%par%t_final = chem_solver_storage%touts(chem_solver_params%n_record_real)
      !
      call disk_save_results_write(fU, c)
      !
      write(fname_pre, '(I4.4, "_", I4.4)') i, j
      write(header, '("n_gas = ", ES12.6)') ch%n_gas
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
  do i=1, chem_solver_params%n_record
    call realtime_heating_cooling_rate(r, chem_solver_params%NEQ, chem_solver_storage%record(:, i))
    write(fU3, FMTstryHistory) chem_solver_storage%touts(i), chem_solver_storage%record(:, i), r
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
    a_disk_ana_params%ana_i_incr = 1+chem_solver_params%n_record/20
  end if
  !
  write(fU1, '(A)') trim(header)
  write(fU2, '(A)') trim(header)
  !
  do k=1, chem_solver_params%n_record, a_disk_ana_params%ana_i_incr
    !+++
    if ((chem_solver_storage%touts(k) .le. 1D2) .or. (chem_solver_storage%touts(k) .ge. 1D5)) then
      cycle
    end if
    if (chem_solver_storage%touts(k) .ge. 4D2) then
      if (mod(k, 50) .ne. 0) then
        cycle
      end if
    end if
    !---
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chem_solver_storage%record(1:chem_species%nSpecies, k) - &
             chem_solver_storage%record(1:chem_species%nSpecies, k-1))) / &
        (chem_solver_storage%record(1:chem_species%nSpecies, k) + &
         chem_solver_storage%record(1:chem_species%nSpecies, k-1) + 1D-15))
      dt_t = (chem_solver_storage%touts(k) - chem_solver_storage%touts(k-1)) / &
             (chem_solver_storage%touts(k) + chem_solver_storage%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    write(fU1, '("time = ", ES14.4)') chem_solver_storage%touts(k)
    !
    chem_solver_storage%y(1:chem_species%nSpecies) = chem_solver_storage%record(1:chem_species%nSpecies, k)
    chem_params%Tgas = chem_solver_storage%record(chem_species%nSpecies+1, k)
    call chem_cal_rates
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chem_solver_storage%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chem_solver_storage%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chem_solver_storage%touts(k)
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
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chem_solver_storage%y(i)
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
