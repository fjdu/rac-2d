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
    ch%ndust_tot = ch%n_gas * ch%ratioDust2HnucNum
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
      chem_params%Av_toISM
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
      chem_params%Av_toISM
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
  allocate(&!c%col_den_acc(chem_idx_some_spe%nItem), &
           !c%col_den(chem_idx_some_spe%nItem), &
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
      ch%ndust_tot = ch%n_gas * ch%ratioDust2HnucNum
      write(*,*) 'Dust density ', ch%ndust_tot
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
      !write(*,'(2ES12.2/)') chem_params%f_selfshielding_H2, chem_params%Av
      call chem_set_solver_flags_alt(1)
      !
      hc_params = ch
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
      hc_params%R_H2_form_rate = &
        get_H2_form_rate( &
          hc_params%R_H2_form_rate_coeff, &
          hc_params%X_gH, &
          hc_params%X_HI, &
          hc_params%n_gas)
      ch%R_H2_form_rate = hc_params%R_H2_form_rate
      !
      call chem_evol_solve
      !
      c%abundances  = chemsol_stor%y(1:chem_species%nSpecies)
      !c%col_den     = c%abundances(chem_idx_some_spe%idx) * c%par%dNcol
      !c%col_den_acc = c%abundances(chem_idx_some_spe%idx) * c%par%Ncol
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
  write(*, '(A/)') 'Doing some analysis... Might be slow.'
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



!subroutine save_fits_cube(filename, im)
!  character(len=*), intent(in) :: filename
!  type(type_image), intent(in) :: im
!  integer stat, fU, blocksize, bitpix, naxis
!  integer, dimension(3) :: naxes
!  integer i, j, group, fpixel, nelements, decimals
!  integer pcount, gcount
!  logical simple, extend
!  !
!  stat=0
!  blocksize = 1
!  pcount = 0
!  gcount = 1
!  group=1
!  fpixel=1
!  decimals = 1
!  author_info_fits = 'fdu@umich.edu'
!  !
!  call ftgiou(fU, stat)
!  !
!  call ftinit(fU, filename, blocksize, stat)
!  !
!  simple=.true.
!  bitpix=-64 ! double
!  naxis=3
!  naxes(1)=im%nx
!  naxes(2)=im%ny
!  naxes(3)=im%nz
!  extend=.true.
!  !
!  call ftphpr(fU, simple, bitpix, naxis, naxes, pcount, gcount, extend, stat)
!  !
!  nelements=naxes(1)*naxes(2)*naxes(3)
!  !
!  call ftpprd(fU, group, fpixel, nelements, im%val, stat)
!  !
!  call ftpkyd(fU, 'BZERO',  0.0D0,  decimals, 'Zero point', stat)
!  call ftpkyd(fU, 'BSCALE', 1.0D0,  decimals, 'Scaling factor', stat)
!  call ftpkyd(fU, 'CDELT1', 1.0D0,  decimals, 'dx', stat)
!  call ftpkyd(fU, 'CDELT2', 1.0D0,  decimals, 'dy', stat)
!  call ftpkyd(fU, 'CDELT3', 1.0D0,  decimals, 'dz', stat)
!  call ftpkyf(fU, 'CRPIX1', 51.0,   decimals, 'i0', stat)
!  call ftpkyf(fU, 'CRPIX2', 51.0,   decimals, 'j0', stat)
!  call ftpkyf(fU, 'CRPIX3', 51.0,   decimals, 'k0', stat)
!  call ftpkyf(fU, 'CRVAL1',  0.0,   decimals, 'x0', stat)
!  call ftpkyf(fU, 'CRVAL2',  0.0,   decimals, 'y0', stat)
!  call ftpkyf(fU, 'CRVAL3',  0.0,   decimals, 'z0', stat)
!  call ftpkys(fU, 'CTYPE1', 'X', '', stat)
!  call ftpkys(fU, 'CTYPE2', 'Y', '', stat)
!  call ftpkys(fU, 'CTYPE3', 'Z', '', stat)
!  call ftpkys(fU, 'AUTHOR', author_info_fits, '', stat)
!  !
!  call ftclos(fU, stat)
!  call ftfiou(fU, stat)
!end subroutine save_fits_cube



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
    write(fU1, '("(", 2F6.2, ")", 2F7.1, 4ES12.2)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av_toStar, &
      leaves%list(idx)%p%par%Av_toISM
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
    write(fU2, '("This (", 2F6.2, ")", 2F7.1, 4ES12.2)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av_toStar, &
      leaves%list(idx)%p%par%Av_toISM
    write(fU2, '("Diff (", 2F6.2, ")", 2F7.1, 4ES12.2)') &
      leaves%list(idx_diff)%p%par%rcen, &
      leaves%list(idx_diff)%p%par%zcen, &
      leaves%list(idx_diff)%p%par%Tgas,  leaves%list(idx_diff)%p%par%Tdust, &
      leaves%list(idx_diff)%p%par%n_gas, leaves%list(idx_diff)%p%par%Ncol, &
      leaves%list(idx_diff)%p%par%Av_toStar, &
      leaves%list(idx_diff)%p%par%Av_toISM
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




subroutine disk_set_disk_params
  ! Background dust
  !a_disk%andrews_dust_bg = a_disk%andrews_gas
  !a_disk%andrews_dust_bg%useNumDens = .false.
  !a_disk%andrews_dust_bg%Md = &
  !  a_disk%andrews_gas%Md * a_disk%dust2gas_mass_bg
  !
  associate( &
    Lstar => a_disk%star_luminosity_in_Lsun * phy_Lsun_CGS, &
    uv2total => a_disk%ratio_uv2total, &
    lyman2uv => a_disk%ratio_lyman2uv, &
    xray2total => a_disk%ratio_xray2total)
    !a_disk%UV_cont_phlumi_star_surface = &
    !  Lstar * uv2total * (1D0 - lyman2uv) / phy_UV_cont_energy_CGS
    !a_disk%Lyman_phlumi_star_surface = &
    !  Lstar * uv2total * lyman2uv         / phy_LyAlpha_energy_CGS
    !a_disk%Xray_phlumi_star_surface  = &
    !  Lstar * xray2total / (xray_energy_kev*1D3*phy_eV2erg)
    !write(str_disp, '("!Stellar total luminosity = ", ES12.4, " erg s-1")') Lstar
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar UV cont luminosity = ", ES12.4, " erg s-1")') &
    !  a_disk%UV_cont_phlumi_star_surface * phy_UV_cont_energy_CGS
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar UV cont photon count rate = ", ES12.4, " s-1")') &
    !  a_disk%UV_cont_phlumi_star_surface
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar LyA luminosity = ", ES12.4, " erg s-1")') &
    !  a_disk%Lyman_phlumi_star_surface * phy_LyAlpha_energy_CGS
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar LyA photon count rate = ", ES12.4, " s-1")') &
    !  a_disk%Lyman_phlumi_star_surface
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar X-ray luminosity = ", ES12.4, " erg s-1")') &
    !  a_disk%Xray_phlumi_star_surface * (xray_energy_kev*1D3*phy_eV2erg)
    !call display_string_both(str_disp, a_book_keeping%fU)
    !write(str_disp, '("!Stellar X-ray photon count rate = ", ES12.4, " s-1")') &
    !  a_disk%Xray_phlumi_star_surface
    !call display_string_both(str_disp, a_book_keeping%fU)
  end associate
end subroutine disk_set_disk_params


