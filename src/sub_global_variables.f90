module phy_const
  implicit none
  double precision, parameter :: phy_Pi = 3.1415926535897932384626433D0
  double precision, parameter :: phy_Pi_2 = 1.57079632679489661923132D0
  double precision, parameter :: phy_2Pi = 6.283185307179586476925D0
  double precision, parameter :: phy_sqrt2Pi = 2.5066282746310005024D0
  double precision, parameter :: phy_NaN = transfer(X'FFFFFFFFFFFFFFFF', 0D0)
  !
  double precision, parameter :: phy_elementaryCharge_SI = 1.602176487D-19
  double precision, parameter :: phy_electronClassicalRadius_SI = 2.8179403267D-15
  double precision, parameter :: phy_electronClassicalRadius_CGS = 2.8179403267D-13
  double precision, parameter :: phy_CoulombConst_SI = 8.9875517873681764D9
  double precision, parameter :: phy_mProton_SI  = 1.67262158D-27 ! kg
  double precision, parameter :: phy_mProton_CGS = 1.67262158D-24 ! g
  double precision, parameter :: phy_mElectron_SI  = 9.10938188E-31 ! kg
  double precision, parameter :: phy_mElectron_CGS = 9.10938188D-28 ! g
  double precision, parameter :: phy_kBoltzmann_SI  = 1.3806503D-23
  double precision, parameter :: phy_kBoltzmann_CGS = 1.3806503D-16
  double precision, parameter :: phy_hPlanck_SI  = 6.62606896D-34
  double precision, parameter :: phy_hPlanck_CGS = 6.62606896D-27
  double precision, parameter :: phy_hbarPlanck_SI  = 1.054571628D-34
  double precision, parameter :: phy_hbarPlanck_CGS = 1.054571628D-27
  double precision, parameter :: phy_GravitationConst_SI  = 6.67428D-11
  double precision, parameter :: phy_GravitationConst_CGS = 6.67428D-8
  double precision, parameter :: phy_SpeedOfLight_SI  = 299792458D0
  double precision, parameter :: phy_SpeedOfLight_CGS = 299792458D2
  double precision, parameter :: phy_StefanBoltzmann_SI = 5.6704D-8
  double precision, parameter :: phy_StefanBoltzmann_CGS = 5.670373D-5
  double precision, parameter :: phy_IdealGasConst_SI = 8.314472D0
  double precision, parameter :: phy_ThomsonScatterCross_CGS = 6.6524574D-25
  !
  double precision, parameter :: phy_Lsun_SI = 3.839D26 ! J s-1
  double precision, parameter :: phy_Lsun_CGS = 3.839D33 ! erg s-1
  double precision, parameter :: phy_Msun_SI = 1.9891D30 ! kg
  double precision, parameter :: phy_Msun_CGS = 1.9891D33 ! g
  double precision, parameter :: phy_Rsun_SI = 6.955D8 ! m
  double precision, parameter :: phy_Rsun_CGS = 6.955D10 ! cm
  !
  double precision, parameter :: phy_Rearth_CGS = 6371D5 ! cm
  double precision, parameter :: phy_Mearth_CGS = 5.97219D27 ! g
  double precision, parameter :: phy_Mmoon_CGS = 7.34767309D25 ! g
  !
  double precision, parameter :: phy_SecondsPerYear = 3600D0*24D0*365D0
  double precision, parameter :: phy_Deg2Rad = phy_Pi/180D0
  double precision, parameter :: phy_erg2joule = 1D-7
  double precision, parameter :: phy_m2cm = 1D2
  double precision, parameter :: phy_kg2g = 1D3
  double precision, parameter :: phy_eV2erg = 1.60217657D-12
  double precision, parameter :: phy_cm_1_2erg = phy_hPlanck_CGS * phy_SpeedOfLight_CGS
  double precision, parameter :: phy_cm_1_2K = phy_cm_1_2erg/phy_kBoltzmann_CGS
  double precision, parameter :: phy_AvogadroConst = 6.02214179D23
  double precision, parameter :: phy_AU2cm = 1.49597871D13
  double precision, parameter :: phy_AU2m  = 1.49597871D11
  double precision, parameter :: phy_pc2m  = 3.08567758D16
  double precision, parameter :: phy_pc2cm = 3.08567758D18
  double precision, parameter :: phy_Angstrom2micron  = 1D-4
  double precision, parameter :: phy_Angstrom2cm  = 1D-8
  double precision, parameter :: phy_micron2cm  = 1D-4
  !
  double precision, parameter :: phy_jansky2CGS = 1D-23
  double precision, parameter :: phy_jansky2SI  = 1D-26
  !
  double precision, parameter :: phy_CMB_T = 2.72548D0
  !
  double precision, parameter :: phy_ratioDust2GasMass_ISM = 0.01D0
  double precision, parameter :: phy_Habing_photon_energy_CGS = 1.99D-11
  double precision, parameter :: phy_LyAlpha_energy_CGS = 1.64D-11
  double precision, parameter :: phy_UV_cont_energy_CGS = phy_Habing_photon_energy_CGS
  double precision, parameter :: phy_Habing_energy_density_CGS = 5.29D-14 ! Draine 2011 book, equation 12.6
  double precision, parameter :: phy_Habing_photon_flux_CGS = 6D7 ! cm-2 s-1
  double precision, parameter :: phy_Habing_energy_flux_CGS = 1.194D-3 ! erg cm-2 s-1
  double precision, parameter :: phy_UVext2Av = 2.6D0 ! Tielens 2005, eq 3.19
  !
  double precision, parameter :: phy_LyAlpha_nu0 = 2.4660718D15
  double precision, parameter :: phy_LyAlpha_l0 = 1215.668D0
  double precision, parameter :: phy_LyAlpha_dnul = 9.938D7
  double precision, parameter :: phy_LyAlpha_f12 = 0.4162D0
  !
  double precision, parameter :: const_LyAlpha_cross_H2O = 1.2D-17 ! Van Dishoeck 2006, Table 1
  double precision, parameter :: const_LyAlpha_cross_OH  = 1.8D-18 ! Van Dishoeck 2006, Table 1
  !
  double precision, parameter :: const_cosmicray_attenuate_N = 5.75D25 ! 96 g cm-2, Nomura 2007
  !
  !double precision :: colDen2Av_coeff = 1D-21 ! Sun Kwok, eq 10.21
  !double precision :: colDen2Av_coeff = 5.3D-22 ! Draine 2011, eq 21.7
end module phy_const
