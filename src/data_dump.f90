module data_dump

use trivials
use data_struct
use grid

implicit none

contains

subroutine back_cells_optical_data(dir_name, dump)
  character(len=*), intent(in) :: dir_name
  logical, intent(in), optional :: dump
  character(len=256) filename
  logical isdump
  integer i, fU, ios
  integer record_len
  !
  if (present(dump)) then
    isdump = dump
  else
    isdump = .true.
  end if
  !
  associate(c => leaves%list(1)%p, &
            lenDouble => kind(0D0), &
            lenInt => kind(0))
    ! Get the size info
    record_len = &
      lenDouble * ( &
          size(c%optical%X) + &
          size(c%optical%summed_ab) + &
          size(c%optical%summed_sc) + &
          size(c%optical%summed) + &
          size(c%optical%acc) + &
          size(c%optical%flux) + &
          size(c%optical%dir_wei%u)*3 + &
          size(c%cont_lut%lam) + &
          size(c%cont_lut%alpha) + &
          size(c%cont_lut%J) &
          ) + &
      lenInt * ( &
          size(c%optical%phc))
  end associate
  !
  filename = combine_dir_filename(dir_name, 'optical_data.bin')
  !
  if (isdump) then
    call my_mkdir(dir_name)
    call openFileBinary(fU, filename, rw='w', &
           record_len=record_len, getu=1)
  else
    call openFileBinary(fU, filename, rw='r', &
           record_len=record_len, getu=1)
  end if
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (isdump) then
        write(fU, rec=i) &
          c%optical%X, &
          c%optical%summed_ab, &
          c%optical%summed_sc, &
          c%optical%summed, &
          c%optical%acc, &
          c%optical%flux, &
          c%optical%phc, &
          c%optical%dir_wei%u, &
          c%optical%dir_wei%v, &
          c%optical%dir_wei%w, &
          c%cont_lut%lam, &
          c%cont_lut%alpha, &
          c%cont_lut%J
      else
        read(fU, rec=i) &
          c%optical%X, &
          c%optical%summed_ab, &
          c%optical%summed_sc, &
          c%optical%summed, &
          c%optical%acc, &
          c%optical%flux, &
          c%optical%phc, &
          c%optical%dir_wei%u, &
          c%optical%dir_wei%v, &
          c%optical%dir_wei%w, &
          c%cont_lut%lam, &
          c%cont_lut%alpha, &
          c%cont_lut%J
      end if
    end associate
  end do
  !
  flush(fU)
  close(fU, iostat=ios, status='KEEP')
end subroutine back_cells_optical_data



subroutine back_cells_chemical_data(dir_name, dump)
  character(len=*), intent(in) :: dir_name
  logical, intent(in), optional :: dump
  character(len=256) filename
  logical isdump
  integer i, fU, ios
  integer record_len
  !
  if (present(dump)) then
    isdump = dump
  else
    isdump = .true.
  end if
  !
  associate(c => leaves%list(1)%p, &
            lenDouble => kind(0D0))
    ! Get the size info
    record_len = &
      lenDouble * ( &
          size(c%abundances) + &
          size(c%col_den_toStar) + &
          size(c%col_den_toISM))
  end associate
  !
  filename = combine_dir_filename(dir_name, 'chemical_data.bin')
  !
  if (isdump) then
    call my_mkdir(dir_name)
    call openFileBinary(fU, filename, rw='w', &
           record_len=record_len, getu=1)
  else
    call openFileBinary(fU, filename, rw='r', &
           record_len=record_len, getu=1)
  end if
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (isdump) then
        write(fU, rec=i) &
          c%abundances, &
          c%col_den_toStar, &
          c%col_den_toISM
      else
        read(fU, rec=i) &
          c%abundances, &
          c%col_den_toStar, &
          c%col_den_toISM
      end if
    end associate
  end do
  !
  flush(fU)
  close(fU, iostat=ios, status='KEEP')
end subroutine back_cells_chemical_data



subroutine back_cells_physical_data(dir_name, dump)
  character(len=*), intent(in) :: dir_name
  logical, intent(in), optional :: dump
  character(len=256) filename
  logical isdump
  integer i, fU, ios
  !
  if (present(dump)) then
    isdump = dump
  else
    isdump = .true.
  end if
  !
  filename = combine_dir_filename(dir_name, 'physical_data.bin')
  !
  if (isdump) then
    call my_mkdir(dir_name)
    call openFileBinary(fU, filename, rw='w', getu=1)
  else
    call openFileBinary(fU, filename, rw='r', getu=1)
  end if
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (isdump) then
        write(fU) &
          c%par%ab_count_dust, &
          c%par%ab_count_water, &
          c%par%sc_count_dust, &
          c%par%sc_count_HI, &
          !
          c%par%ndustcompo, &
          !
          c%par%rmin, &
          c%par%rmax, &
          c%par%rcen, &
          c%par%dr, &
          c%par%zmin, &
          c%par%zmax, &
          c%par%zcen, &
          c%par%dz, &
          c%par%volume, &
          c%par%surf_area, &
          c%par%area_T, &
          c%par%area_B, &
          c%par%area_I, &
          c%par%area_O, &
          c%par%Tgas, &
          c%par%Tdust, &
          !
          c%par%n_gas, &
          !
          c%par%mgas_cell, &
          !
          c%par%Tdusts, &
          c%par%rho_dusts, &
          c%par%n_dusts, &
          c%par%mp_dusts, &
          c%par%mdusts_cell, &
          !
          c%par%en_exchange_per_vol, &
          c%par%en_exchange, &
          c%par%en_exchange_tot, &
          !
          c%par%abso_wei, &
          !
          c%par%sig_dusts, &
          c%par%en_gains, &
          c%par%en_gains_abso, &
          c%par%en_prevs, &
          c%par%kphs, &
          !
          c%par%en_gain_tot, &
          c%par%en_gain_abso_tot, &
          !
          c%par%sigdust_ave, &
          c%par%ndust_tot, &
          c%par%mdust_tot, &
          !
          c%par%UV_G0_factor_background, &
          !
          c%par%Ncol, &
          c%par%dNcol, &
          !
          c%par%ab_en_water, &
          !
          c%par%phflux_Lya, &
          !
          c%par%flux_UV_star_unatten, &
          c%par%flux_Lya_star_unatten, &
          c%par%flux_Vis_star_unatten, &
          !
          c%par%G0_UV_toISM, &
          c%par%G0_UV_toStar, &
          c%par%G0_Lya_atten, &
          !
          c%par%Av_toISM, &
          c%par%Av_toStar, &
          !
          c%par%Ncol_toISM, &
          c%par%Ncol_toStar, &
          !
          c%par%omega_albedo, &
          c%par%zeta_cosmicray_H2, &
          !
          c%par%zeta_Xray_H2, &
          !
          c%par%R_H2_form_rate_coeff, &
          c%par%R_H2_form_rate, &
          !
          c%par%f_selfshielding_toISM_H2, &
          c%par%f_selfshielding_toISM_CO, &
          c%par%f_selfshielding_toISM_H2O, &
          c%par%f_selfshielding_toISM_OH, &
          !
          c%par%f_selfshielding_toStar_H2, &
          c%par%f_selfshielding_toStar_CO, &
          c%par%f_selfshielding_toStar_H2O, &
          c%par%f_selfshielding_toStar_OH, &
          !
          c%par%SitesPerGrain, &
          c%par%GrainMaterialDensity_CGS, &
          c%par%GrainRadius_CGS, &
          !
          c%par%ratioDust2GasMass, &
          c%par%ratioDust2HnucNum, &
          c%par%dust_depletion, &
          c%par%MeanMolWeight, &
          !
          c%par%omega_Kepler, &
          c%par%velo_Kepler, &
          c%par%velo_gradient, &
          c%par%velo_width_turb, &
          c%par%coherent_length, &
          c%par%sound_speed, &
          !
          c%par%alpha_viscosity, &
          c%par%ambipolar_f, &
          c%par%ion_charge, &
          c%par%Neufeld_G, &
          c%par%Neufeld_dv_dz, &
          !
          c%par%t_final, &
          !
          c%par%X_H2, &
          c%par%X_HI, &
          c%par%X_CI, &
          c%par%X_Cplus, &
          c%par%X_OI, &
          c%par%X_CO, &
          c%par%X_H2O, &
          c%par%X_OH, &
          c%par%X_E, &
          c%par%X_Hplus, &
          c%par%X_gH, &
          c%par%flux_tot, &
          c%par%flux_Xray, &
          c%par%flux_UV, &
          c%par%flux_Lya, &
          c%par%flux_Vis, &
          c%par%flux_NIR, &
          c%par%flux_MIR, &
          c%par%flux_FIR, &
          c%par%dir_tot_r, &
          c%par%dir_tot_z, &
          c%par%dir_Xray_r, &
          c%par%dir_Xray_z, &
          c%par%dir_UV_r,  &
          c%par%dir_UV_z,  &
          c%par%dir_Lya_r, &
          c%par%dir_Lya_z, &
          c%par%dir_Vis_r, &
          c%par%dir_Vis_z, &
          c%par%dir_NIR_r, &
          c%par%dir_NIR_z, &
          c%par%dir_MIR_r, &
          c%par%dir_MIR_z, &
          c%par%dir_FIR_r, &
          c%par%dir_FIR_z, &
          c%par%aniso_tot, &
          c%par%aniso_Xray, &
          c%par%aniso_UV, &
          c%par%aniso_Lya, &
          c%par%aniso_Vis, &
          c%par%aniso_NIR, &
          c%par%aniso_MIR, &
          c%par%aniso_FIR, &
          c%par%pressure_thermal, &
          c%par%gravity_z, &
          c%par%gravity_acc_z
      else
        read(fU) &
          c%par%ab_count_dust, &
          c%par%ab_count_water, &
          c%par%sc_count_dust, &
          c%par%sc_count_HI, &
          !
          c%par%ndustcompo, &
          !
          c%par%rmin, &
          c%par%rmax, &
          c%par%rcen, &
          c%par%dr, &
          c%par%zmin, &
          c%par%zmax, &
          c%par%zcen, &
          c%par%dz, &
          c%par%volume, &
          c%par%surf_area, &
          c%par%area_T, &
          c%par%area_B, &
          c%par%area_I, &
          c%par%area_O, &
          c%par%Tgas, &
          c%par%Tdust, &
          !
          c%par%n_gas, &
          !
          c%par%mgas_cell, &
          !
          c%par%Tdusts, &
          c%par%rho_dusts, &
          c%par%n_dusts, &
          c%par%mp_dusts, &
          c%par%mdusts_cell, &
          !
          c%par%en_exchange_per_vol, &
          c%par%en_exchange, &
          c%par%en_exchange_tot, &
          !
          c%par%abso_wei, &
          !
          c%par%sig_dusts, &
          c%par%en_gains, &
          c%par%en_gains_abso, &
          c%par%en_prevs, &
          c%par%kphs, &
          !
          c%par%en_gain_tot, &
          c%par%en_gain_abso_tot, &
          !
          c%par%sigdust_ave, &
          c%par%ndust_tot, &
          c%par%mdust_tot, &
          !
          c%par%UV_G0_factor_background, &
          !
          c%par%Ncol, &
          c%par%dNcol, &
          !
          c%par%ab_en_water, &
          !
          c%par%phflux_Lya, &
          !
          c%par%flux_UV_star_unatten, &
          c%par%flux_Lya_star_unatten, &
          c%par%flux_Vis_star_unatten, &
          !
          c%par%G0_UV_toISM, &
          c%par%G0_UV_toStar, &
          c%par%G0_Lya_atten, &
          !
          c%par%Av_toISM, &
          c%par%Av_toStar, &
          !
          c%par%Ncol_toISM, &
          c%par%Ncol_toStar, &
          !
          c%par%omega_albedo, &
          c%par%zeta_cosmicray_H2, &
          !
          c%par%zeta_Xray_H2, &
          !
          c%par%R_H2_form_rate_coeff, &
          c%par%R_H2_form_rate, &
          !
          c%par%f_selfshielding_toISM_H2, &
          c%par%f_selfshielding_toISM_CO, &
          c%par%f_selfshielding_toISM_H2O, &
          c%par%f_selfshielding_toISM_OH, &
          !
          c%par%f_selfshielding_toStar_H2, &
          c%par%f_selfshielding_toStar_CO, &
          c%par%f_selfshielding_toStar_H2O, &
          c%par%f_selfshielding_toStar_OH, &
          !
          c%par%SitesPerGrain, &
          c%par%GrainMaterialDensity_CGS, &
          c%par%GrainRadius_CGS, &
          !
          c%par%ratioDust2GasMass, &
          c%par%ratioDust2HnucNum, &
          c%par%dust_depletion, &
          c%par%MeanMolWeight, &
          !
          c%par%omega_Kepler, &
          c%par%velo_Kepler, &
          c%par%velo_gradient, &
          c%par%velo_width_turb, &
          c%par%coherent_length, &
          c%par%sound_speed, &
          !
          c%par%alpha_viscosity, &
          c%par%ambipolar_f, &
          c%par%ion_charge, &
          c%par%Neufeld_G, &
          c%par%Neufeld_dv_dz, &
          !
          c%par%t_final, &
          !
          c%par%X_H2, &
          c%par%X_HI, &
          c%par%X_CI, &
          c%par%X_Cplus, &
          c%par%X_OI, &
          c%par%X_CO, &
          c%par%X_H2O, &
          c%par%X_OH, &
          c%par%X_E, &
          c%par%X_Hplus, &
          c%par%X_gH, &
          c%par%flux_tot, &
          c%par%flux_Xray, &
          c%par%flux_UV, &
          c%par%flux_Lya, &
          c%par%flux_Vis, &
          c%par%flux_NIR, &
          c%par%flux_MIR, &
          c%par%flux_FIR, &
          c%par%dir_tot_r, &
          c%par%dir_tot_z, &
          c%par%dir_Xray_r, &
          c%par%dir_Xray_z, &
          c%par%dir_UV_r,  &
          c%par%dir_UV_z,  &
          c%par%dir_Lya_r, &
          c%par%dir_Lya_z, &
          c%par%dir_Vis_r, &
          c%par%dir_Vis_z, &
          c%par%dir_NIR_r, &
          c%par%dir_NIR_z, &
          c%par%dir_MIR_r, &
          c%par%dir_MIR_z, &
          c%par%dir_FIR_r, &
          c%par%dir_FIR_z, &
          c%par%aniso_tot, &
          c%par%aniso_Xray, &
          c%par%aniso_UV, &
          c%par%aniso_Lya, &
          c%par%aniso_Vis, &
          c%par%aniso_NIR, &
          c%par%aniso_MIR, &
          c%par%aniso_FIR, &
          c%par%pressure_thermal, &
          c%par%gravity_z, &
          c%par%gravity_acc_z
      end if
    end associate
  end do
  !
  flush(fU)
  close(fU, iostat=ios, status='KEEP')
end subroutine back_cells_physical_data



subroutine back_cells_physical_data_aux(dir_name, dump)
  character(len=*), intent(in) :: dir_name
  logical, intent(in), optional :: dump
  character(len=256) filename
  logical isdump
  integer i, fU, ios
  integer record_len
  !
  if (present(dump)) then
    isdump = dump
  else
    isdump = .true.
  end if
  !
  associate(c => leaves%list(1)%p, &
            lenDouble => kind(0D0))
    ! Get the size info
    record_len = &
      lenDouble * (size(c%val) + 25)
  end associate
  !
  filename = combine_dir_filename(dir_name, 'physical_data_aux.bin')
  !
  if (isdump) then
    call my_mkdir(dir_name)
    call openFileBinary(fU, filename, rw='w', &
           record_len=record_len, getu=1)
  else
    call openFileBinary(fU, filename, rw='r', &
           record_len=record_len, getu=1)
  end if
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (isdump) then
        write(fU, rec=i) &
          c%val, &
          c%h_c_rates%hc_net_rate, &
          c%h_c_rates%heating_photoelectric_small_grain_rate, &
          c%h_c_rates%heating_formation_H2_rate, &
          c%h_c_rates%heating_cosmic_ray_rate, &
          c%h_c_rates%heating_vibrational_H2_rate, &
          c%h_c_rates%heating_ionization_CI_rate, &
          c%h_c_rates%heating_photodissociation_H2_rate, &
          c%h_c_rates%heating_photodissociation_H2O_rate, &
          c%h_c_rates%heating_photodissociation_OH_rate, &
          c%h_c_rates%heating_Xray_Bethell_rate, &
          c%h_c_rates%heating_viscosity_rate, &
          c%h_c_rates%heating_chem, &
          c%h_c_rates%cooling_photoelectric_small_grain_rate, &
          c%h_c_rates%cooling_vibrational_H2_rate, &
          c%h_c_rates%cooling_gas_grain_collision_rate, &
          c%h_c_rates%cooling_OI_rate, &
          c%h_c_rates%cooling_CII_rate, &
          c%h_c_rates%cooling_Neufeld_H2O_rate_rot, &
          c%h_c_rates%cooling_Neufeld_H2O_rate_vib, &
          c%h_c_rates%cooling_Neufeld_CO_rate_rot, &
          c%h_c_rates%cooling_Neufeld_CO_rate_vib, &
          c%h_c_rates%cooling_Neufeld_H2_rot_rate, &
          c%h_c_rates%cooling_LymanAlpha_rate, &
          c%h_c_rates%cooling_free_bound_rate, &
          c%h_c_rates%cooling_free_free_rate
      else
        read(fU, rec=i) &
          c%val, &
          c%h_c_rates%hc_net_rate, &
          c%h_c_rates%heating_photoelectric_small_grain_rate, &
          c%h_c_rates%heating_formation_H2_rate, &
          c%h_c_rates%heating_cosmic_ray_rate, &
          c%h_c_rates%heating_vibrational_H2_rate, &
          c%h_c_rates%heating_ionization_CI_rate, &
          c%h_c_rates%heating_photodissociation_H2_rate, &
          c%h_c_rates%heating_photodissociation_H2O_rate, &
          c%h_c_rates%heating_photodissociation_OH_rate, &
          c%h_c_rates%heating_Xray_Bethell_rate, &
          c%h_c_rates%heating_viscosity_rate, &
          c%h_c_rates%heating_chem, &
          c%h_c_rates%cooling_photoelectric_small_grain_rate, &
          c%h_c_rates%cooling_vibrational_H2_rate, &
          c%h_c_rates%cooling_gas_grain_collision_rate, &
          c%h_c_rates%cooling_OI_rate, &
          c%h_c_rates%cooling_CII_rate, &
          c%h_c_rates%cooling_Neufeld_H2O_rate_rot, &
          c%h_c_rates%cooling_Neufeld_H2O_rate_vib, &
          c%h_c_rates%cooling_Neufeld_CO_rate_rot, &
          c%h_c_rates%cooling_Neufeld_CO_rate_vib, &
          c%h_c_rates%cooling_Neufeld_H2_rot_rate, &
          c%h_c_rates%cooling_LymanAlpha_rate, &
          c%h_c_rates%cooling_free_bound_rate, &
          c%h_c_rates%cooling_free_free_rate
      end if
    end associate
  end do
  !
  flush(fU)
  close(fU, iostat=ios, status='KEEP')
end subroutine back_cells_physical_data_aux



end module data_dump
