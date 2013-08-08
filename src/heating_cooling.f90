module heating_cooling

use phy_const
use data_struct

implicit none

type, extends(type_cell_rz_phy_basic) :: type_heating_cooling_parameters
  double precision :: X_H2, X_HI, X_CI, X_Cplus, X_OI, X_CO, X_H2O, X_OH, X_E, X_Hplus, X_gH
  double precision :: Neufeld_G = 1D0, Neufeld_dv_dz
end type type_heating_cooling_parameters


type(type_heating_cooling_parameters) :: heating_cooling_params

type(type_heating_cooling_rates_list) :: heating_cooling_rates


contains


function heating_photoelectric_small_grain()
  ! Wolfire 1995, which is actually taken from
  ! Bakes 1994, equation 42 and 43.
  ! It is for very small graphitic grains and PAHs.
  ! Applicable for T up to 1D4 K.
  ! Note that in Bakes 1994 they made important changes from their
  ! section 4 to section 5.
  !
  ! Output unit = erg s-1 cm-3
  !
  double precision heating_photoelectric_small_grain
  associate( &
    chi => heating_cooling_params%UV_G0_factor * exp(-phy_UVext2Av*heating_cooling_params%Av), &
    n_gas   => heating_cooling_params%n_gas, &
    n_e   => heating_cooling_params%X_E * heating_cooling_params%n_gas, &
    Tgas  => heating_cooling_params%Tgas)
    associate(tmp => chi*sqrt(Tgas)/(n_e))
      heating_photoelectric_small_grain = &
        1D-24 * chi * n_gas * &
        ( &
        4.87D-2 / (1D0 + 4D-3 * (tmp)**0.73D0) &
        + &
        3.65D-2 * (1D-4*Tgas)**0.7 / (1D0 + 2D-4*tmp))
    end associate
  end associate
end function heating_photoelectric_small_grain


function heating_formation_H2()
  ! Rollig 2006
  ! Sternberg 1989, equation D5
  ! 2.4D-12 erg = 1/3 * 4.5 eV
  double precision heating_formation_H2
    heating_formation_H2 = &
      2.4D-12 * heating_cooling_params%R_H2_form_rate
end function heating_formation_H2


function heating_cosmic_ray()
  ! Bruderer 2009
  !Return unit = erg s-1 cm-3
  ! 1.5D-11 = 9 * 1.6D-19 / 1D-7
  double precision heating_cosmic_ray
  associate( &
        n_gas             => heating_cooling_params%n_gas, &
        zeta_cosmicray_H2 => heating_cooling_params%zeta_cosmicray_H2)
    heating_cosmic_ray = 1.5D-11 * zeta_cosmicray_H2 * n_gas
  end associate
end function heating_cosmic_ray


function heating_vibrational_H2()
  ! Rollig 2006, equation C.2 - C.3
  ! It causes problems: T cannot be solved if it is included.
  ! 9.4D-22 =  2.9e-10*23500.0*1.38e-16
  ! It seems this function duplicates another one.
  !   Yes, the other one is redundant.  That one is a
  !   special case of this one; there the f factor in
  !   Rollig 2006 is taken to be 2.8D-3.
  double precision heating_vibrational_H2
  associate( &
        n_gas   => heating_cooling_params%n_gas, &
        X_H2    => heating_cooling_params%X_H2, &
        chi => heating_cooling_params%UV_G0_factor * &
               heating_cooling_params%f_selfshielding_H2 * &
               exp(-phy_UVext2Av*heating_cooling_params%Av), &
        gamma_10 => 5.4D-13 * sqrt(heating_cooling_params%Tgas))
    heating_vibrational_H2 = &
      (n_gas * X_H2) * chi * 9.4D-22 / &
      (1D0 + (1.9D-6 + chi*4.7D-10) / (n_gas * gamma_10))
  end associate
end function heating_vibrational_H2


function heating_photodissociation_H2()
  ! Tielens 2005, P72, equation 3.18, 3.19
  double precision heating_photodissociation_H2
  associate( &
        chi => heating_cooling_params%UV_G0_factor * &
               heating_cooling_params%f_selfshielding_H2 * &
               exp(-phy_UVext2Av*heating_cooling_params%Av), &
        n_gas => heating_cooling_params%n_gas, &
        X_H2  => heating_cooling_params%X_H2)
    heating_photodissociation_H2 = &
      4D-14 * (n_gas*X_H2) * 3.4D-10 * chi
end associate
end function heating_photodissociation_H2


function heating_photodissociation_H2O()
  ! LyAlpha = 1.63e-11 erg = 10.2 eV
  ! H(H2O) - H(OH) - H(H) = 498.826e3 J mol-1 = 8.282e-12 erg = 5.18 eV.
  double precision heating_photodissociation_H2O
  associate( &
        chi => heating_cooling_params%LymanAlpha_flux_0 * &
          heating_cooling_params%f_selfshielding_H2O * &
          exp(-phy_UVext2Av * heating_cooling_params%Av), &
        LyAlpha_cross_H2O => const_LyAlpha_cross_H2O, &
        ph_disso_en_H2O   => 8.07D-12, &
        n_H2O => heating_cooling_params%n_gas * heating_cooling_params%X_H2O)
    heating_photodissociation_H2O = &
      ph_disso_en_H2O * n_H2O * LyAlpha_cross_H2O * chi
  end associate
end function heating_photodissociation_H2O


function heating_photodissociation_OH()
  ! LyAlpha = 1.63e-11 erg = 10.2 eV
  ! H(OH) - H(O) - H(H) = 428.188 J mol-1 = 7.11e-12 erg = 4.44 eV.
  double precision heating_photodissociation_OH
  associate( &
        chi => heating_cooling_params%LymanAlpha_flux_0 * &
          heating_cooling_params%f_selfshielding_OH * &
          exp(-phy_UVext2Av * heating_cooling_params%Av), &
        LyAlpha_cross_OH => const_LyAlpha_cross_OH, &
        ph_disso_en_OH   => 9.19D-12, &
        n_OH => heating_cooling_params%n_gas * heating_cooling_params%X_OH)
    heating_photodissociation_OH = &
      ph_disso_en_OH * n_OH * LyAlpha_cross_OH * chi
  end associate
end function heating_photodissociation_OH


function heating_ionization_CI()
  ! Tielens 2005, P66, equation 3.8
  double precision heating_ionization_CI
  associate( &
        chi => heating_cooling_params%UV_G0_factor * exp(-phy_UVext2Av*heating_cooling_params%Av), &
        n_gas => heating_cooling_params%n_gas, &
        X_CI  => heating_cooling_params%X_CI)
    heating_ionization_CI = &
      2.2D-22 * X_CI * n_gas * chi
  end associate
end function heating_ionization_CI


function heating_Xray_bethell()
  use load_Bethell_Xray_cross
  double precision heating_Xray_bethell
  double precision :: en_X = 1D0! keV
  double precision :: en_deposit = 18D0 * phy_eV2erg ! 18 eV; AGN paper
  double precision sigma
  sigma = sigma_Xray_Bethell(en_X, &
    heating_cooling_params%dust_depletion, &
    heating_cooling_params%ratioDust2HnucNum, &
    heating_cooling_params%GrainRadius_CGS)
  heating_Xray_bethell = sigma * heating_cooling_params%n_gas * en_deposit * &
    heating_cooling_params%Xray_flux_0 * exp(-sigma*heating_cooling_params%Ncol)
end function heating_Xray_bethell


function heating_viscosity()
  ! From the AGN paper.
  double precision heating_viscosity
  associate( &
        rho => heating_cooling_params%n_gas * phy_mProton_CGS * &
               heating_cooling_params%MeanMolWeight, &
        c2  => phy_kBoltzmann_CGS * heating_cooling_params%Tgas &
               / (phy_mProton_CGS * heating_cooling_params%MeanMolWeight))
    heating_viscosity = 2.25D0 * heating_cooling_params%alpha_viscosity * rho * c2 &
      * heating_cooling_params%omega_Kepler
  end associate
end function heating_viscosity


function cooling_photoelectric_small_grain()
  ! Wolfire 1995, which is actually taken from
  ! Bakes 1994, equation 44.
  !
  ! Note that the coefficient in Wolfire 1995 is 4.65D-30.
  ! Note sure why.
  !
  ! Output unit = erg s-1 cm-3
  double precision cooling_photoelectric_small_grain
  associate( &
    chi => heating_cooling_params%UV_G0_factor * exp(-phy_UVext2Av*heating_cooling_params%Av), &
    n_gas => heating_cooling_params%n_gas, &
    n_e   => heating_cooling_params%X_E * heating_cooling_params%n_gas, &
    Tgas  => heating_cooling_params%Tgas)
    associate(tmp => chi*sqrt(Tgas)/(n_e))
      cooling_photoelectric_small_grain = &
        3.49D-30 * Tgas**0.944D0 * tmp**(0.735D0 * Tgas**(-0.068)) * n_e * n_gas
    end associate
  end associate
end function cooling_photoelectric_small_grain


function cooling_LymanAlpha()
  ! Tielens 2005, page 53, equation 2.62
  double precision cooling_LymanAlpha
  associate( &
        Tgas => heating_cooling_params%Tgas, &
        n_HI => heating_cooling_params%n_gas * heating_cooling_params%X_HI, &
        n_E  => heating_cooling_params%n_gas * heating_cooling_params%X_E)
    cooling_LymanAlpha = &
      7.3D-19 * n_e * n_HI * exp(-118400D0 / Tgas)
  end associate
end function cooling_LymanAlpha


function cooling_free_bound()
  ! Draine 2011, page 139, equation 14.5 (optically thin case); and page 320, equation 27.22, 27.23
  double precision cooling_free_bound
  associate( &
        T   => heating_cooling_params%Tgas, &
        n_p => heating_cooling_params%n_gas * heating_cooling_params%X_Hplus, &
        n_E => heating_cooling_params%n_gas * heating_cooling_params%X_E, &
        T4  => heating_cooling_params%Tgas/1D4, &
        ZZ  => 1D0)
    cooling_free_bound = n_E * n_p * &
      4.13D-13 * ZZ * (T4/ZZ)**(-0.7131D0-0.0115D0*log(T4/ZZ)) * & ! alpha_A
      (0.787D0 - 0.0230D0*log(T4/ZZ)) * phy_kBoltzmann_CGS * T ! Err_A
  end associate
end function cooling_free_bound


function cooling_free_free()
  ! Tielens 2005, page 53, equation 2.65
  ! Essentially a simpler version of Draine equation (10.12).
  double precision cooling_free_free
  associate( &
        T   => heating_cooling_params%Tgas, &
        n_p => heating_cooling_params%n_gas * heating_cooling_params%X_Hplus, &
        n_E => heating_cooling_params%n_gas * heating_cooling_params%X_E, &
        ZZ  => 1D0)
    cooling_free_free = &
      1.4D-27 * ZZ * sqrt(T) * n_E * n_p
  end associate
end function cooling_free_free


function cooling_vibrational_H2()
  ! Rollig 2006, equation C.1
  ! 8.26D-13 = 5988 * 1.38D-16
  double precision cooling_vibrational_H2
  associate( &
        n_gas   => heating_cooling_params%n_gas, &
        Tgas    => heating_cooling_params%Tgas, &
        X_H2    => heating_cooling_params%X_H2, &
        chi => heating_cooling_params%UV_G0_factor * &
               heating_cooling_params%f_selfshielding_H2 * &
               exp(-phy_UVext2Av*heating_cooling_params%Av), &
        gamma_10 => 5.4D-13 * sqrt(heating_cooling_params%Tgas), &
        A_10 => 8.6D-7, &
        D_1 => 2.6D-11)
    cooling_vibrational_H2 = &
      8.26D-13 * gamma_10 * exp(-5988D0/Tgas) * &
        (n_gas*n_gas*X_H2) * &
        (A_10 + chi * D_1) / (gamma_10*n_gas + A_10 + chi*D_1)
  end associate
end function cooling_vibrational_H2


function cooling_Neufeld_H2_rot()
  use load_Neufeld_cooling_H2
  double precision cooling_Neufeld_H2_rot
  associate( &
    T => a_Neufeld_cooling_H2_params%T)
    T = heating_cooling_params%Tgas
    call get_H2_rot_cool_params
    associate( &
        n_H2  => heating_cooling_params%n_gas * heating_cooling_params%X_H2, &
        L0    => a_Neufeld_cooling_H2_params%L0, &
        L_LTE => a_Neufeld_cooling_H2_params%L_LTE, &
        n_12  => a_Neufeld_cooling_H2_params%n_12, &
        alpha => a_Neufeld_cooling_H2_params%alpha)
      if (alpha .GT. 0D0) then
        cooling_Neufeld_H2_rot = n_H2 * n_H2 / &
          (1D0/L0 + n_H2/L_LTE + &
           1D0/L0 * (n_H2/n_12)**alpha * (1D0 - n_12*L0/L_LTE))
      else
        cooling_Neufeld_H2_rot = n_H2 * n_H2 / (1D0/L0 + n_H2/L_LTE)
      end if
    end associate
  end associate
end function cooling_Neufeld_H2_rot


function cooling_gas_grain_collision()
  ! Rollig 2006, equation A.7
  ! Z should be the abundance of grains relative to a certain value
  double precision cooling_gas_grain_collision
  associate( &
        n_gas => heating_cooling_params%n_gas, &
        Tgas  => heating_cooling_params%Tgas, &
        Tdust => heating_cooling_params%Tdust, &
        r_g   => heating_cooling_params%GrainRadius_CGS, &
        Z     => heating_cooling_params%dust_depletion) ! Todo
    !cooling_gas_grain_collision = &
    !  3.5D-34 * sqrt(Tgas) * (Tgas - Tdust) * n_gas * n_gas * Z
    ! Hollenbach 1989, eq 2.15
    cooling_gas_grain_collision = &
      1.2D-31 * n_gas * n_gas * sqrt(Tgas/1D3 * 1D-6/r_g) * &
      (1D0 - 0.8D0*exp(-75D0/Tgas)) * (Tgas - Tdust) * Z
  end associate
end function cooling_gas_grain_collision


function cooling_OI()
  ! Rollig 2006, equation A.5 and A.6
  ! Z: abundance of O in 3D-4
  ! Need to use the real abundance of OI.
  double precision cooling_OI
  associate( &
    n_gas => heating_cooling_params%n_gas, &
    Tgas  => heating_cooling_params%Tgas, &
    Z     => heating_cooling_params%X_OI / 1.76D-4)
    associate( &
      tmp1 => n_gas + 0.5D0 * 1.66D-5 / (1.35D-11 * Tgas**0.45D0), &
      tmp2 => n_gas + 0.5D0 * 8.46D-5 / (4.37D-12 * Tgas**0.66D0), &
      tmp3 => exp(98D0/Tgas), &
      tmp4 => exp(228D0/Tgas))
      associate( &
        tmp5 => n_gas*n_gas + tmp3 * tmp1 * (3D0*n_gas + tmp4*5D0*tmp2))
        associate( &
            cooling_OI_63  => 3.15D-14 * 8.46D-5 * 0.5D0 * Z * &
              3D-4 * n_gas * tmp3 * 3D0 * n_gas * tmp1 &
              / tmp5, &
            cooling_OI_146 => 1.35D-14 * 1.66D-5 * 0.5D0 * Z * &
              3D-4 * n_gas * n_gas * n_gas &
              / tmp5)
          cooling_OI = cooling_OI_63 + cooling_OI_146
        end associate
      end associate
    end associate
  end associate
end function cooling_OI


function cooling_CII()
  ! Rollig 2006, equation A.2
  ! Z: Abudance of C in 1.4D-4
  ! Need to use the real abundance of CII
  ! Note that the value used in Rollig 2006 is 1.4D-4.
  ! I guess this is because they normalize relative to n(H2).
  double precision cooling_CII
  associate( &
    Tgas    => heating_cooling_params%Tgas, &
    n_gas   => heating_cooling_params%n_gas, &
    Z       => heating_cooling_params%X_Cplus/7.3D-5)
    associate( &
        cooling_CII_158 => 2.02D-24 * n_gas * Z / &
          (1D0 + 0.5D0 * exp(92D0/Tgas) * (1D0 + 1300D0/n_gas)))
      cooling_CII = cooling_CII_158
    end associate
  end associate
end function cooling_CII


function cooling_CI()
  double precision cooling_CI
  cooling_CI = 0D0
end function cooling_CI


function cooling_Neufeld_H2O_rot()
  use load_Neufeld_cooling_H2O
  double precision cooling_Neufeld_H2O_rot
  associate( &
        L0    => a_Neufeld_cooling_H2O_params%L0, &
        L_LTE => a_Neufeld_cooling_H2O_params%L_LTE, &
        n_12  => a_Neufeld_cooling_H2O_params%n_12, &
        alpha => a_Neufeld_cooling_H2O_params%alpha, &
        T     => a_Neufeld_cooling_H2O_params%T, &
        log10N=> a_Neufeld_cooling_H2O_params%log10N, &
        G     => heating_cooling_params%Neufeld_G, &
        n_M   => heating_cooling_params%n_gas * heating_cooling_params%X_H2O, &
        n_H2  => heating_cooling_params%n_gas * heating_cooling_params%X_H2, &
        dv_dz => heating_cooling_params%Neufeld_dv_dz)
    T = heating_cooling_params%Tgas
    log10N = log10(G * n_M / dv_dz)
    L0    = get_L0()
    L_LTE = get_L_LTE()
    n_12  = get_n_12()
    alpha = get_alpha()
    !
    cooling_Neufeld_H2O_rot = n_H2 * n_M / &
      (1D0/L0 + n_H2/L_LTE + &
       1D0/L0 * (n_H2/n_12)**alpha * (1D0 - n_12*L0/L_LTE))
  end associate
end function cooling_Neufeld_H2O_rot


function cooling_Neufeld_H2O_vib()
  use load_Neufeld_cooling_H2O
  double precision cooling_Neufeld_H2O_vib
  associate( &
    L0        => a_Neufeld_cooling_H2O_params%L0_vib, &
    L_LTE     => a_Neufeld_cooling_H2O_params%L_LTE_vib, &
    T         => a_Neufeld_cooling_H2O_params%T, &
    log10N    => a_Neufeld_cooling_H2O_params%log10N, &
    G         => heating_cooling_params%Neufeld_G, &
    n_M       => heating_cooling_params%n_gas * heating_cooling_params%X_H2O, &
    n_H2      => heating_cooling_params%n_gas * heating_cooling_params%X_H2, &
    dv_dz     => heating_cooling_params%Neufeld_dv_dz)
    !
    T      = heating_cooling_params%Tgas
    log10N = log10(G * n_M / dv_dz)
    L0    = get_L0_vib()
    L_LTE = get_L_LTE_vib()
    !
    cooling_Neufeld_H2O_vib = n_H2 * n_M / (1D0/L0 + n_H2/L_LTE)
  end associate
end function cooling_Neufeld_H2O_vib


function cooling_Neufeld_CO_rot()
  use load_Neufeld_cooling_CO
  double precision cooling_Neufeld_CO_rot
  associate( &
    ! L is the cooling coefficient.
    ! L     => a_Neufeld_cooling_CO_params%L, &
    ! The following four are to be interpolated or extrapolated from the Neufeld tables.
    L0    => a_Neufeld_cooling_CO_params%L0, &
    L_LTE => a_Neufeld_cooling_CO_params%L_LTE, &
    n_12  => a_Neufeld_cooling_CO_params%n_12, &
    alpha => a_Neufeld_cooling_CO_params%alpha, &
    T     => a_Neufeld_cooling_CO_params%T, &
    log10N=> a_Neufeld_cooling_CO_params%log10N, &
    G     => heating_cooling_params%Neufeld_G, &
    n_M   => heating_cooling_params%n_gas * heating_cooling_params%X_CO, &
    n_H2  => heating_cooling_params%n_gas * heating_cooling_params%X_H2, &
    dv_dz => heating_cooling_params%Neufeld_dv_dz)
    !
    T = heating_cooling_params%Tgas
    log10N = log10(G * n_M / dv_dz)
    L0    = get_L0()
    L_LTE = get_L_LTE()
    n_12  = get_n_12()
    alpha = get_alpha()
    !
    cooling_Neufeld_CO_rot = n_H2 * n_M / &
      (1D0/L0 + n_H2/L_LTE + &
       1D0/L0 * (n_H2/n_12)**alpha * (1D0 - n_12*L0/L_LTE))
  end associate
end function cooling_Neufeld_CO_rot


function cooling_Neufeld_CO_vib()
  use load_Neufeld_cooling_CO
  double precision cooling_Neufeld_CO_vib
  associate( &
    L0     => a_Neufeld_cooling_CO_params%L0_vib, &
    L_LTE  => a_Neufeld_cooling_CO_params%L_LTE_vib, &
    T      => a_Neufeld_cooling_CO_params%T, &
    log10N => a_Neufeld_cooling_CO_params%log10N, &
    G      => heating_cooling_params%Neufeld_G, &
    n_M    => heating_cooling_params%n_gas * heating_cooling_params%X_CO, &
    n_H2   => heating_cooling_params%n_gas * heating_cooling_params%X_H2, &
    dv_dz  => heating_cooling_params%Neufeld_dv_dz)
    !
    T      = heating_cooling_params%Tgas
    log10N = log10(G * n_M / dv_dz)
    L0    = get_L0_vib()
    L_LTE = get_L_LTE_vib()
    !
    cooling_Neufeld_CO_vib = n_H2 * n_M / (1D0/L0 + n_H2/L_LTE)
  end associate
end function cooling_Neufeld_CO_vib


function heating_minus_cooling()
  double precision heating_minus_cooling
  associate(r => heating_cooling_rates)
    r%heating_photoelectric_small_grain_rate = heating_photoelectric_small_grain()
    r%heating_formation_H2_rate              = heating_formation_H2()
    r%heating_cosmic_ray_rate                = heating_cosmic_ray()
    r%heating_vibrational_H2_rate            = heating_vibrational_H2()
    r%heating_ionization_CI_rate             = heating_ionization_CI()
    r%heating_photodissociation_H2_rate      = heating_photodissociation_H2()
    r%heating_photodissociation_H2O_rate     = heating_photodissociation_H2O()
    r%heating_photodissociation_OH_rate      = heating_photodissociation_OH()
    r%heating_Xray_Bethell_rate              = heating_Xray_Bethell()
    r%heating_viscosity_rate                 = heating_viscosity()
    r%cooling_photoelectric_small_grain_rate = cooling_photoelectric_small_grain()
    r%cooling_vibrational_H2_rate            = cooling_vibrational_H2()
    r%cooling_gas_grain_collision_rate       = cooling_gas_grain_collision()
    r%cooling_OI_rate                        = cooling_OI()
    r%cooling_CII_rate                       = cooling_CII()
    r%cooling_Neufeld_H2O_rate_rot           = cooling_Neufeld_H2O_rot()
    r%cooling_Neufeld_H2O_rate_vib           = cooling_Neufeld_H2O_vib()
    r%cooling_Neufeld_CO_rate_rot            = cooling_Neufeld_CO_rot()
    r%cooling_Neufeld_CO_rate_vib            = cooling_Neufeld_CO_vib()
    r%cooling_Neufeld_H2_rot_rate            = cooling_Neufeld_H2_rot()
    r%cooling_LymanAlpha_rate                = cooling_LymanAlpha()
    r%cooling_free_bound_rate                = cooling_free_bound()
    r%cooling_free_free_rate                 = cooling_free_free()
    heating_minus_cooling = &
        r%heating_photoelectric_small_grain_rate &  ! 1
      + r%heating_formation_H2_rate &               ! 2
      + r%heating_cosmic_ray_rate &                 ! 3
      + r%heating_vibrational_H2_rate &             ! 4
      + r%heating_ionization_CI_rate &              ! 5
      + r%heating_photodissociation_H2_rate &       ! 6
      + r%heating_photodissociation_H2O_rate &      ! 7
      + r%heating_photodissociation_OH_rate &       ! 8
      + r%heating_Xray_Bethell_rate         &       ! 9
      + r%heating_viscosity_rate            &       ! 10
      - r%cooling_photoelectric_small_grain_rate &  ! 1
      - r%cooling_vibrational_H2_rate &             ! 2
      - r%cooling_gas_grain_collision_rate &        ! 3
      - r%cooling_OI_rate &                         ! 4
      - r%cooling_CII_rate &                        ! 5
      - r%cooling_Neufeld_H2O_rate_rot &            ! 6
      - r%cooling_Neufeld_H2O_rate_vib &            ! 7
      - r%cooling_Neufeld_CO_rate_rot &             ! 8
      - r%cooling_Neufeld_CO_rate_vib &             ! 9
      - r%cooling_Neufeld_H2_rot_rate &             ! 10
      - r%cooling_LymanAlpha_rate &                 ! 11
      - r%cooling_free_bound_rate &                 ! 12
      - r%cooling_free_free_rate                    ! 13
  end associate
end function heating_minus_cooling


function heating_cooling_dust_solve_T(T0, n_iter, converged)
  double precision heating_cooling_dust_solve_T
  ! The problem is to choose a wise initial point.
  integer :: i, ii, n_max_iter=64
  double precision, intent(in) :: T0
  integer, intent(out), optional :: n_iter
  logical, intent(out), optional :: converged
  double precision Tn_m1, Tn, Tn_p1, fn_m1, fn, fn_p1
  double precision T_save
  double precision :: rtol = 1D-4
  double precision :: atol = 1D0
  integer, parameter :: n_max_iter_outer = 19
  double precision, dimension(n_max_iter_outer), parameter :: &
    ratios_try = (/1D0, 1.2D0, 0.8D0, 1.5D0, 0.5D0, 3D0, 0.3D0, &
        1D1, 1D-1, 1D2, 1D-2, 1D3, 1D-3, 1D4, 1D-4, 1D5, 1D-5, 1D6, 1D-6/)
  T_save = heating_cooling_params%Tgas
  if (present(converged)) then
    converged = .false.
  end if
  do ii=1, n_max_iter_outer
    Tn_m1 = T0 * ratios_try(ii)
    Tn    = Tn_m1 * 1.2D0 + 10D0
    heating_cooling_params%Tgas = Tn_m1
    fn_m1 = heating_minus_cooling()
    heating_cooling_params%Tgas = Tn
    fn = heating_minus_cooling()
    do i=1, n_max_iter
      Tn_p1 = Tn - fn * ((Tn - Tn_m1) / (fn - fn_m1))
      write(*,*) i, Tn_m1, Tn, Tn_p1, fn
      if (isnan(Tn_p1) .or. (Tn_p1 .le. 0D0)) then
        exit
      end if
      if (abs(Tn_p1 - Tn) .LE. (atol + rtol * (Tn_p1 + Tn))) then
        exit
      end if
      Tn_m1 = Tn
      Tn = Tn_p1
      fn_m1 = fn
      heating_cooling_params%Tgas = Tn_p1
      fn = heating_minus_cooling()
    end do
    ! If the result is illegitimate, redo the iteration.
    !if (.NOT. (isnan(Tn_p1) .OR. (Tn_p1 .LE. 0D0))) then
    if (.NOT. (isnan(Tn_p1) .or. (Tn_p1 .le. 0D0))) then
      if (present(converged)) then
        converged = .true.
      end if
      exit
    end if
  end do
  heating_cooling_dust_solve_T = Tn_p1
  heating_cooling_params%Tgas = T_save
  if (present(n_iter)) then
    n_iter = i
  end if
end function heating_cooling_dust_solve_T


function solve_bisect_T(T0, n_iter, converged)
  implicit none
  double precision solve_bisect_T, T0
  integer n_iter
  logical converged
  double precision T_save
  double precision x1, x2, f1, f2, dx, xmid, fmid
  integer i, j
  integer :: nmax_expand = 1024
  integer :: nmax_shrink = 64
  double precision :: expand_factor = 0.5D0
  double precision :: expand_factor_tmp, x_tmp
  logical found_diff_sign
  double precision :: rtol = 1D-2
  double precision :: atol = 1D0
  T_save = heating_cooling_params%Tgas
  x1 = T0 / 1.1D0
  x2 = T0 * 1.1D0
  heating_cooling_params%Tgas = x1
  f1 = heating_minus_cooling()
  heating_cooling_params%Tgas = x2
  f2 = heating_minus_cooling()
  found_diff_sign = .false.
  converged = .false.
  solve_bisect_T = -1D0
  do i=1, nmax_expand
    if (f1*f2 .le. 0D0) then
      found_diff_sign = .true.
      exit
    end if
    if (abs(f1) .lt. abs(f2)) then
      x1 = max(1D0, x1 + expand_factor * (x1 - x2))
      heating_cooling_params%Tgas = x1
      f1 = heating_minus_cooling()
    else if (abs(f1) .ge. abs(f2)) then
      x2 = max(1D0, x2 + expand_factor * (x2 - x1))
      heating_cooling_params%Tgas = x2
      f2 = heating_minus_cooling()
    else
      if (isnan(f1)) then
        x1 = max(1D0, (x1 + expand_factor*x2)/(1D0+expand_factor))
        expand_factor_tmp = expand_factor / 2D0
        do j=1, 30
          x_tmp = max(1D0, x1 + expand_factor_tmp * (x1 - x2))
          heating_cooling_params%Tgas = x_tmp
          f1 = heating_minus_cooling()
          if (.not. isnan(f1)) then
            x1 = x_tmp
            exit
          else
            x_tmp = x1
            expand_factor_tmp  = expand_factor_tmp / 2D0
          end if
        end do
      else
        x2 = max(1D0, (x2 + expand_factor*x1)/(1D0+expand_factor))
        expand_factor_tmp = expand_factor / 2D0
        do j=1, 30
          x_tmp = max(1D0, x2 + expand_factor_tmp * (x2 - x1))
          heating_cooling_params%Tgas = x_tmp
          f2 = heating_minus_cooling()
          if (.not. isnan(f2)) then
            x2 = x_tmp
            exit
          else
            x_tmp = x2
            expand_factor_tmp  = expand_factor_tmp / 2D0
          end if
        end do
      end if
    end if
  end do
  if (.not. found_diff_sign) then
    write(*,*) 'Fail to find an appropriate bracketing interval.'
    converged = .false.
    heating_cooling_params%Tgas = T_save
    n_iter = 0
    return
  else
    if (f1 .eq. 0D0) then
      solve_bisect_T = x1
      converged = .true.
      heating_cooling_params%Tgas = T_save
      n_iter = 0
      return
    else if (f2 .eq. 0D0) then
      solve_bisect_T = x2
      converged = .true.
      heating_cooling_params%Tgas = T_save
      n_iter = 0
      return
    else
      do i=1, nmax_shrink
        dx = x2 - x1
        xmid = 0.5D0 * (x1+x2)
        if (abs(dx) .lt. (rtol*xmid + atol)) then
          solve_bisect_T = xmid
          converged = .true.
          heating_cooling_params%Tgas = T_save
          n_iter = i
          return
        end if
        heating_cooling_params%Tgas = xmid
        fmid = heating_minus_cooling()
        if (isnan(fmid)) then
          xmid = x2 - dx * 0.1D0
          heating_cooling_params%Tgas = xmid
          fmid = heating_minus_cooling()
        end if
        if (fmid*f1 .lt. 0D0) then
          x2 = xmid
          f2 = fmid
        else if (fmid .eq. 0D0) then
          solve_bisect_T = xmid
          converged = .true.
          heating_cooling_params%Tgas = T_save
          n_iter = i
          return
        else
          x1 = xmid
          f1 = fmid
        end if
      end do
    end if
  end if
end function solve_bisect_T



end module heating_cooling
