module data_dump

use trivials
use data_struct
use grid

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
          ab_count_dust, ab_count_water, &
          sc_count_dust, sc_count_HI, &
          !
          ndustcompo, &
          !
          rmin, rmax, rcen, dr, zmin, zmax, zcen, dz, &
          volume, surf_area, area_T, area_B, area_I, area_O, &
          Tgas, &
          Tdust, &
          !
          n_gas, &
          !
          mgas_cell, &
          !
          Tdusts, &
          rho_dusts, &
          n_dusts, &
          mp_dusts, &
          mdusts_cell, &
          !
          en_exchange_per_vol, &
          en_exchange, &
          en_exchange_tot, &
          !
          abso_wei, &
          !
          sig_dusts, &
          en_gains, &
          en_gains_abso, &
          en_prevs, &
          kphs, &
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
          Ncol, &
          dNcol, &
          !
          ab_en_water, &
          !
          phflux_Lya, &
          !
          flux_UV_star_unatten, &
          flux_Lya_star_unatten, &
          flux_Vis_star_unatten, &
          !
          G0_UV_toISM, &
          G0_UV_toStar, &
          G0_Lya_atten, &
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
          dust_depletion, &
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
          !
          X_H2, X_HI, X_CI, X_Cplus, X_OI, X_CO, &
          X_H2O, X_OH, X_E, X_Hplus, X_gH, &
          flux_tot, flux_Xray, flux_UV, flux_Lya, &
          flux_Vis, flux_NIR, flux_MIR, flux_FIR, &
          dir_tot_r, dir_tot_z, dir_Xray_r, dir_Xray_z, &
          dir_UV_r,  dir_UV_z,  dir_Lya_r, dir_Lya_z, &
          dir_Vis_r, dir_Vis_z, dir_NIR_r, dir_NIR_z, &
          dir_MIR_r, dir_MIR_z, dir_FIR_r, dir_FIR_z, &
          aniso_tot, aniso_Xray, aniso_UV, aniso_Lya, &
          aniso_Vis, aniso_NIR, aniso_MIR, aniso_FIR, &
          pressure_thermal, gravity_z, gravity_acc_z
      else
        read(fU) &
          ab_count_dust, ab_count_water, &
          sc_count_dust, sc_count_HI, &
          !
          ndustcompo, &
          !
          rmin, rmax, rcen, dr, zmin, zmax, zcen, dz, &
          volume, surf_area, area_T, area_B, area_I, area_O, &
          Tgas, &
          Tdust, &
          !
          n_gas, &
          !
          mgas_cell, &
          !
          Tdusts, &
          rho_dusts, &
          n_dusts, &
          mp_dusts, &
          mdusts_cell, &
          !
          en_exchange_per_vol, &
          en_exchange, &
          en_exchange_tot, &
          !
          abso_wei, &
          !
          sig_dusts, &
          en_gains, &
          en_gains_abso, &
          en_prevs, &
          kphs, &
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
          Ncol, &
          dNcol, &
          !
          ab_en_water, &
          !
          phflux_Lya, &
          !
          flux_UV_star_unatten, &
          flux_Lya_star_unatten, &
          flux_Vis_star_unatten, &
          !
          G0_UV_toISM, &
          G0_UV_toStar, &
          G0_Lya_atten, &
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
          dust_depletion, &
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
          !
          X_H2, X_HI, X_CI, X_Cplus, X_OI, X_CO, &
          X_H2O, X_OH, X_E, X_Hplus, X_gH, &
          flux_tot, flux_Xray, flux_UV, flux_Lya, &
          flux_Vis, flux_NIR, flux_MIR, flux_FIR, &
          dir_tot_r, dir_tot_z, dir_Xray_r, dir_Xray_z, &
          dir_UV_r,  dir_UV_z,  dir_Lya_r, dir_Lya_z, &
          dir_Vis_r, dir_Vis_z, dir_NIR_r, dir_NIR_z, &
          dir_MIR_r, dir_MIR_z, dir_FIR_r, dir_FIR_z, &
          aniso_tot, aniso_Xray, aniso_UV, aniso_Lya, &
          aniso_Vis, aniso_NIR, aniso_MIR, aniso_FIR, &
          pressure_thermal, gravity_z, gravity_acc_z
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
