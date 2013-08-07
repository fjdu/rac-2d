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
  character(len=32):: backup_src_pattern = '*.f90 *.py'
  !double precision :: colDen2Av_coeff = 1D-21 ! Sun Kwok, eq 10.21
  !double precision :: colDen2Av_coeff = 5.3D-22 ! Draine 2011, eq 21.7
end type phy_chem_rad_disk_params


type :: a__disk
  type(phy_chem_rad_disk_params) :: params
end type a__disk


type :: disk_iteration_params
  integer :: n_iter=128, n_iter_used
  double precision :: rtol_T = 0.1D0, atol_T = 2D0
  double precision :: rtol_abun = 0.2D0, atol_abun = 1D-12
  logical flag_converged
  integer n_cell_converged
  real converged_cell_percentage_stop
  character(len=128) iter_files_dir
  character(len=256) notes
  logical :: flag_log_rates = .FALSE.
  logical :: flag_shortcut_ini = .FALSE.
  logical :: redo_couple_every_column = .FALSE.
  logical :: redo_couple_every_cell = .FALSE.
  logical :: iter_cell_outwards = .TRUE.
  logical :: iter_cell_upwards = .TRUE.
end type disk_iteration_params


type :: disk_iteration_storage
  double precision, dimension(:,:), allocatable :: T_s
  double precision, dimension(:,:,:), allocatable :: abundances
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

integer, dimension(:), allocatable :: calculating_cells_list
integer n_calculating_cells, n_calculating_cells_max

character(len=128) :: filename_save_results
integer fU_save_results

double precision, parameter, private :: const_geometric_factor_UV   = 0.01D0
double precision, parameter, private :: const_geometric_factor_Xray = 0.001D0
double precision, parameter, private :: ratioDust2GasMass_ISM = 0.01D0
double precision, parameter, private :: xray_energy_kev = 1D0

integer, parameter :: len_item=14

namelist /disk_configure/ &
  disk_params_ini

namelist /cell_configure/ &
  cell_params_ini

namelist /iteration_configure/ &
  a_disk_iter_params


contains


subroutine disk_iteration
  use my_timer
  type(date_time) a_date_time
  integer i, i0, i_count, l_count, ii
  double precision Told, Tnew
  logical isTgood
  !
  call disk_iteration_prepare
  !
  call save_post_config_params
  !
  ! Now start the major big loop.
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
    i_count = 0
    l_count = 0
    do
      l_count = l_count + 1
      call update_params_layers_this
      !
      do i=1, n_calculating_cells
        i_count = i_count + 1
        i0 = calculating_cells_list(i)
        !
        write(*, '(3(A, I4, A, I6, 2X), (A, I4), 2X, A, 4F10.3)') &
          "Iter", a_disk_iter_params%n_iter_used, "/", a_disk_iter_params%n_iter, &
          "Cell", i_count, '/', cell_leaves%nlen, &
          "cell", i, '/', n_calculating_cells, &
          "Layer", l_count, &
          'r-z', &
          cell_leaves%list(i0)%p%par%rmin, &
          cell_leaves%list(i0)%p%par%rmax, &
          cell_leaves%list(i0)%p%par%zmin, &
          cell_leaves%list(i0)%p%par%zmax
        write(*, '(A, F12.3, 2X, A, ES12.4, 4X, 2A, X, A/)') &
          'Tgas_old: ', cell_leaves%list(i0)%p%par%Tgas, &
          'n_gas: ', cell_leaves%list(i0)%p%par%n_gas, &
          'CMD: ', trim(a_disk%params%filename_exe), &
          trim(a_disk_iter_params%iter_files_dir)
        !
        call calc_this_cell(i0)
        !
        write(*, '(A, 2X, 10ES12.4)') &
          'Abundances:', chem_solver_storage%y(chem_idx_some_spe%idx(1:10))
        !
        write(*, '(A, F12.3/)') &
          'Tgas_new: ', cell_leaves%list(i0)%p%par%Tgas
        !
        call disk_save_results_write(i0)
        flush(fU_save_results)
        !
        if (a_disk_iter_params%flag_log_rates) then
          call save_chem_rates(i0)
        end if
      end do
      !
      call update_calculating_cells_list
      !
      if (n_calculating_cells .eq. 0) then
        exit
      end if
    end do
    !
    call check_convergency
    !
    write(fU_save_results, '(A, L)') '! flag_converged = ', a_disk_iter_params%flag_converged
    write(fU_save_results, '(A)') '! Finish saving ' // trim(filename_save_results)
    write(fU_save_results, '(A)') '! at ' // trim(a_date_time%date_time_str())
    flush(fU_save_results)
    close(fU_save_results)
    !
    if (a_disk_iter_params%flag_converged) then
      exit
    end if
  end do
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=iteration_configure)
    write(a_book_keeping%fU, '("Number of cells =", I4)') cell_leaves%nlen
    flush(a_book_keeping%fU)
  end if
end subroutine disk_iteration


subroutine disk_iteration_prepare
  integer i
  !
  call make_grid
  n_calculating_cells_max = cell_leaves%nlen / 2
  allocate(calculating_cells_list(n_calculating_cells_max))
  !
  call chem_read_reactions()
  call chem_load_reactions()
  call chem_parse_reactions()
  call chem_get_dupli_reactions()
  call chem_get_idx_for_special_species()
  !
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, '("!", A, 2X, I5)') 'Number of cells (leaf):', cell_leaves%nlen
    write(a_book_keeping%fU, '("!", A, 2X, I5)') 'Number of cells (total):', root%nOffspring
    write(a_book_keeping%fU, '("!", A, 2X, I5)') 'Number of reactions:', chem_net%nReactions
    write(a_book_keeping%fU, '("!", A, 2X, I5)') 'Number of species ', chem_species%nSpecies
    flush(a_book_keeping%fU)
  end if
  !
  call disk_set_disk_params
  call disk_set_gridcell_params
  !
  ! Load cell parameters from output of RADMC.
  ! Tgas, Tdust and n_gas are obtained from the RADMC output file.
  ! n_gas is calculated with a dust to gas ratio.
  ! Tgas is set to Tdust.
  ! It overrides the values set by disk_set_cell_params.
  !call load_dust_dens_T_from_RADMC
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  if (.NOT. allocated(a_disk_iter_storage%T_s)) then
    allocate(a_disk_iter_storage%T_s(cell_leaves%nlen, &
                                     0:a_disk_iter_params%n_iter), &
             a_disk_iter_storage%abundances(chem_idx_some_spe%nItem, &
                                            cell_leaves%nlen, &
                                            0:a_disk_iter_params%n_iter))
  end if
  !
  call chem_load_initial_abundances
  call disk_set_cell_init_abundances
  !
  do i=1, cell_leaves%nlen
    a_disk_iter_storage%T_s(i, 0) = cell_leaves%list(i)%p%par%Tgas
    a_disk_iter_storage%abundances(:, i, 0) = chem_solver_storage%y(chem_idx_some_spe%idx)
  end do
  !call disk_calc_column_densities
  call disk_calc_disk_mass
end subroutine disk_iteration_prepare


subroutine update_params_layers_this
  use load_Visser_CO_selfshielding
  integer i, i0, j, j0
  do i=1, n_calculating_cells
    i0 = calculating_cells_list(i)
    associate(p  => cell_leaves%list(i0)%p, &
              dz => cell_leaves%list(i0)%p%par%dz * phy_AU2cm)
      p%col_den   = p%abundances(chem_idx_some_spe%idx) * p%par%n_gas * dz
      p%par%dNcol = p%par%n_gas * dz
      if (cell_leaves%list(i0)%p%above%n .gt. 0) then
        do j=1, cell_leaves%list(i0)%p%above%n
          j0 = cell_leaves%list(i0)%p%above%idx(j)
          p%col_den_acc = (cell_leaves%list(j0)%p%col_den_acc + &
                           cell_leaves%list(j0)%p%col_den) * p%above%fra(j)
          p%par%Ncol    = (cell_leaves%list(j0)%p%par%Ncol + &
                           cell_leaves%list(j0)%p%par%dNcol) * p%above%fra(j)
        end do
      else
        p%col_den_acc = 0D0
        p%par%Ncol    = 0D0
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
        (p%par%Ncol + p%par%dNcol * 0.5D0)
      p%par%f_selfshielding_H2  = &
        min(1D0, ((Ncol_H2 + dcol_H2*0.5D0)/1D14)**(-0.75D0)) ! Tielens 2005, equation 8.39
      p%par%f_selfshielding_H2O = &
        min(1D0, exp(-(Ncol_H2O * const_LyAlpha_cross_H2O))) * &
        tau2beta(dcol_H2O * const_LyAlpha_cross_H2O)
      p%par%f_selfshielding_OH  = &
        min(1D0, exp(-(Ncol_OH * const_LyAlpha_cross_OH))) * &
        tau2beta(dcol_OH * const_LyAlpha_cross_OH)
      p%par%f_selfshielding_CO = min(1D0, get_12CO_shielding(Ncol_H2, Ncol_CO))
    end associate
  end do
end subroutine update_params_layers_this


subroutine update_params_above(i0)
  use load_Visser_CO_selfshielding
  integer, intent(in) :: i0
  integer j, j0
  associate(p  => cell_leaves%list(i0)%p, &
            dz => cell_leaves%list(i0)%p%par%dz * phy_AU2cm)
    p%col_den   = p%abundances(chem_idx_some_spe%idx) * p%par%n_gas * dz
    p%par%dNcol = p%par%n_gas * dz
    if (cell_leaves%list(i0)%p%above%n .gt. 0) then
      do j=1, cell_leaves%list(i0)%p%above%n
        j0 = cell_leaves%list(i0)%p%above%idx(j)
        p%col_den_acc = (cell_leaves%list(j0)%p%col_den_acc + &
                         cell_leaves%list(j0)%p%col_den) * p%above%fra(j)
        p%par%Ncol    = (cell_leaves%list(j0)%p%par%Ncol + &
                         cell_leaves%list(j0)%p%par%dNcol) * p%above%fra(j)
      end do
    else
      p%col_den_acc = 0D0
      p%par%Ncol    = 0D0
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
      (p%par%Ncol + p%par%dNcol * 0.5D0)
    p%par%f_selfshielding_H2  = &
      min(1D0, ((Ncol_H2 + dcol_H2*0.5D0)/1D14)**(-0.75D0)) ! Tielens 2005, equation 8.39
    p%par%f_selfshielding_H2O = &
      min(1D0, exp(-(Ncol_H2O * const_LyAlpha_cross_H2O))) * &
      tau2beta(dcol_H2O * const_LyAlpha_cross_H2O)
    p%par%f_selfshielding_OH  = &
      min(1D0, exp(-(Ncol_OH * const_LyAlpha_cross_OH))) * &
      tau2beta(dcol_OH * const_LyAlpha_cross_OH)
    p%par%f_selfshielding_CO = min(1D0, get_12CO_shielding(Ncol_H2, Ncol_CO))
  end associate
end subroutine update_params_above


subroutine check_convergency
  integer ii
  ii = a_disk_iter_params%n_iter_used
  a_disk_iter_params%n_cell_converged = &
    count( &
      (abs(a_disk_iter_storage%abundances(:, :, ii) &
          - a_disk_iter_storage%abundances(:, :, ii-1)) &
       -  (a_disk_iter_params%atol_abun + &
           a_disk_iter_params%rtol_abun * &
             abs(a_disk_iter_storage%abundances(:, :, ii) &
               + a_disk_iter_storage%abundances(:, :, ii-1)) &
           )) .le. 0D0)
  a_disk_iter_params%flag_converged = &
    a_disk_iter_params%n_cell_converged .le. &
    a_disk_iter_params%converged_cell_percentage_stop * real(cell_leaves%nlen)
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, '("! Number of cells converged: ", I6, "/", I6)') &
      a_disk_iter_params%n_cell_converged, cell_leaves%nlen
  end if
end subroutine check_convergency


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
  calculating_cells_list(1:n) = list_tmp(1:n)
end subroutine update_calculating_cells_list


subroutine calc_this_cell(id)
  integer, intent(in) :: id
  integer i, i0, ntmp
  double precision Tnew, Told
  logical found_neighbor, isTgood
    ! Try to set the initial condition for chemical evolution
  if (a_disk_iter_params%n_iter_used .eq. 1) then
    found_neighbor = .false.
    do i=1, cell_leaves%list(id)%p%around%n
      i0 = cell_leaves%list(id)%p%around%idx(i)
      if (cell_leaves%list(i0)%p%iIter .gt. cell_leaves%list(id)%p%iIter) then
        chem_solver_storage%y = cell_leaves%list(i0)%p%abundances
        cell_leaves%list(id)%p%par%Tgas = cell_leaves%list(i0)%p%par%Tgas
        found_neighbor = .true.
        exit
      end if
    end do
    if (.not. found_neighbor) then
      call chem_load_initial_abundances
    end if
  else
    chem_solver_storage%y = cell_leaves%list(id)%p%abundances
  end if
  cell_leaves%list(id)%p%iIter = a_disk_iter_params%n_iter_used
  !
  call update_params_above(id)
  !
  call set_chemistry_params_from_cell(id)
  call chem_cal_rates
  call chem_set_solver_flags
  call chem_evol_solve
  !
  cell_leaves%list(id)%p%abundances = chem_solver_storage%y
  !
  a_disk_iter_storage%abundances(:, id, a_disk_iter_params%n_iter_used) = &
    chem_solver_storage%y(chem_idx_some_spe%idx)
  !
  call set_heatingcooling_params_from_cell(id)
  Told = cell_leaves%list(id)%p%par%Tgas
  Tnew = solve_bisect_T(Told, ntmp, isTgood)
  !
  if (.not. isTgood) then
    write(*, '(/A/)') 'Heating/cooling not converged! Reset to previous temperature.'
  else
    cell_leaves%list(id)%p%par%Tgas = Tnew
  end if
  !
  a_disk_iter_storage%T_s(id, a_disk_iter_params%n_iter_used) = cell_leaves%list(id)%p%par%Tgas
  !
  cell_leaves%list(id)%p%h_c_rates = heating_cooling_rates
end subroutine calc_this_cell


subroutine disk_save_results_pre
  integer i
  character(len=64) fmt_str
  character(len=1024) tmp_str
  if (.NOT. getFileUnit(fU_save_results)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  write(filename_save_results, '("iter_", I4.4, ".dat")') a_disk_iter_params%n_iter_used
  filename_save_results = combine_dir_filename(a_disk_iter_params%iter_files_dir, filename_save_results)
  call openFileSequentialWrite(fU_save_results, filename_save_results, 99999)
  !
  write(fmt_str, '("(", I4, "A14)")') chem_idx_some_spe%nItem
  write(tmp_str, fmt_str) chem_idx_some_spe%names
  write(fU_save_results, '(A)') &
    '!' // &
    str_pad_to_len('i', 3) // &
    str_pad_to_len('rmin',    len_item) // &
    str_pad_to_len('rmax',    len_item) // &
    str_pad_to_len('zmin',    len_item) // &
    str_pad_to_len('zmax',    len_item) // &
    str_pad_to_len('Tgas',    len_item) // &
    str_pad_to_len('Tdust',   len_item) // &
    str_pad_to_len('n_gas',   len_item) // &
    str_pad_to_len('Av',      len_item) // &
    str_pad_to_len('UV_G0',   len_item) // &
    str_pad_to_len('LyA0 ',   len_item) // &
    str_pad_to_len('XRay0',   len_item) // &
    str_pad_to_len('Ncol ',   len_item) // &
    str_pad_to_len('dNcol',   len_item) // &
    str_pad_to_len('f_H2',    len_item) // &
    str_pad_to_len('f_H2O',   len_item) // &
    str_pad_to_len('f_OH',    len_item) // &
    str_pad_to_len('f_CO',    len_item) // &
    str_pad_to_len('R_H2_fo', len_item) // &
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
end subroutine disk_save_results_pre


subroutine disk_save_results_write(i0)
  character(len=64) fmt_str
  integer, intent(in) :: i0
  write(fmt_str, '(", ", I4, "ES14.4E4)")') chem_idx_some_spe%nItem
  associate(c => cell_leaves%list(i0)%p)
    write(fU_save_results, '(I4, 41ES14.4E4' // trim(fmt_str)) &
    i0, &
    c%par%rmin                                             , &
    c%par%rmax                                             , &
    c%par%zmin                                             , &
    c%par%zmax                                             , &
    c%par%Tgas                                             , &
    c%par%Tdust                                            , &
    c%par%n_gas                                            , &
    c%par%Av                                               , &
    c%par%UV_G0_factor                                     , &
    c%par%LymanAlpha_flux_0                                , &
    c%par%Xray_flux_0                                      , &
    c%par%Ncol                                             , &
    c%par%dNcol                                            , &
    c%par%f_selfshielding_H2                               , &
    c%par%f_selfshielding_H2O                              , &
    c%par%f_selfshielding_OH                               , &
    c%par%f_selfshielding_CO                               , &
    c%par%R_H2_form_rate                                   , &
    c%h_c_rates%heating_photoelectric_small_grain_rate , &
    c%h_c_rates%heating_formation_H2_rate              , &
    c%h_c_rates%heating_cosmic_ray_rate                , &
    c%h_c_rates%heating_vibrational_H2_rate            , &
    c%h_c_rates%heating_ionization_CI_rate             , &
    c%h_c_rates%heating_photodissociation_H2_rate      , &
    c%h_c_rates%heating_photodissociation_H2O_rate     , &
    c%h_c_rates%heating_photodissociation_OH_rate      , &
    c%h_c_rates%heating_Xray_Bethell_rate              , &
    c%h_c_rates%heating_viscosity_rate                 , &
    c%h_c_rates%cooling_photoelectric_small_grain_rate , &
    c%h_c_rates%cooling_vibrational_H2_rate            , &
    c%h_c_rates%cooling_gas_grain_collision_rate       , &
    c%h_c_rates%cooling_OI_rate                        , &
    c%h_c_rates%cooling_CII_rate                       , &
    c%h_c_rates%cooling_Neufeld_H2O_rate_rot           , &
    c%h_c_rates%cooling_Neufeld_H2O_rate_vib           , &
    c%h_c_rates%cooling_Neufeld_CO_rate_rot            , &
    c%h_c_rates%cooling_Neufeld_CO_rate_vib            , &
    c%h_c_rates%cooling_Neufeld_H2_rot_rate            , &
    c%h_c_rates%cooling_LymanAlpha_rate                , &
    c%h_c_rates%cooling_free_bound_rate                , &
    c%h_c_rates%cooling_free_free_rate                 , &
    c%abundances(chem_idx_some_spe%idx)
  end associate
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
end subroutine set_heatingcooling_params_from_cell


subroutine set_chemistry_params_from_cell(id)
  integer id
  chem_params => cell_leaves%list(id)%p%par
end subroutine set_chemistry_params_from_cell


subroutine disk_set_cell_init_abundances
  integer i
  do i=1, cell_leaves%nlen
    cell_leaves%list(i)%p%abundances = chem_solver_storage%y
  end do
end subroutine disk_set_cell_init_abundances


subroutine disk_set_gridcell_params
  integer i
  do i=1, cell_leaves%nlen
    associate(c => cell_leaves%list(i)%p)
      if (.not. associated(cell_leaves%list(i)%p%par)) then
        allocate(cell_leaves%list(i)%p%par)
      end if
      if (.not. allocated(cell_leaves%list(i)%p%h_c_rates)) then
        allocate(cell_leaves%list(i)%p%h_c_rates)
      end if
      if (.not. allocated(cell_leaves%list(i)%p%iIter)) then
        allocate(cell_leaves%list(i)%p%iIter)
        cell_leaves%list(i)%p%iIter = 0
      end if
      !
      allocate(c%abundances(chem_species%nSpecies), &
               c%col_den(chem_idx_some_spe%nItem), &
               c%col_den_acc(chem_idx_some_spe%nItem))
      !
      cell_leaves%list(i)%p%par = cell_params_ini
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
      c%par%dust_depletion = c%par%ratioDust2GasMass / ratioDust2GasMass_ISM
      c%par%n_dust = c%par%n_gas * c%par%ratioDust2HnucNum
      c%par%velo_width_turb = 1D5 ! Todo
      !
      c%par%UV_G0_factor = 1D0 + &
        a_disk%params%UV_cont_phlumi_star_surface &
           / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
           / phy_Habing_photon_flux_CGS &
           * const_geometric_factor_UV
      c%par%LymanAlpha_flux_0 = &
        a_disk%params%Lyman_phlumi_star_surface &
           / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2)
      c%par%Xray_flux_0 = &
        a_disk%params%Xray_phlumi_star_surface &
           / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
           * const_geometric_factor_Xray
      !c%par%cosmicray_flux_top = 1D0
    end associate
    associate( &
      G     => phy_GravitationConst_CGS, &
      M     => a_disk%params%star_mass_in_Msun * phy_Msun_CGS, &
      r     => cell_leaves%list(i)%p%par%rcen * phy_AU2cm, &
      v     => cell_leaves%list(i)%p%par%velo_Kepler, &
      w     => cell_leaves%list(i)%p%par%omega_Kepler, &
      dv_dr => cell_leaves%list(i)%p%par%velo_gradient)
      v = sqrt(G * M / r)
      w = v / r
      dv_dr = 0.5D0 * v / r
    end associate
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


function tau2beta(tau)
  double precision tau, tau2beta
  if (tau .LE. 1D-4) then
    tau2beta = 1D0
  else
    tau2beta = (1D0 - exp(-tau)) / tau
  end if
end function tau2beta


subroutine save_post_config_params
  type(phy_chem_rad_disk_params) disk_params_tmp
  namelist /disk_params_log/ disk_params_tmp
  disk_params_tmp = a_disk%params
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=disk_params_log)
    flush(a_book_keeping%fU)
  end if
end subroutine save_post_config_params

end module disk
