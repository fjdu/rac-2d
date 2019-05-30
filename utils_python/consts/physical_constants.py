phy_Pi = 3.1415926535897932384626433e0
phy_Pi_2 = 1.57079632679489661923132e0
phy_2Pi = 6.283185307179586476925e0
phy_sqrt2Pi = 2.5066282746310005024e0

phy_max_exp = 222e0

phy_elementaryCharge_SI = 1.602176487e-19
phy_electronClassicalRadius_SI = 2.8179403267e-15
phy_electronClassicalRadius_CGS = 2.8179403267e-13
phy_CoulombConst_SI = 8.9875517873681764e9
phy_mProton_SI  = 1.67262158e-27 # kg
phy_mProton_CGS = 1.67262158e-24 # g
phy_mElectron_SI  = 9.10938188E-31 # kg
phy_mElectron_CGS = 9.10938188e-28 # g
phy_kBoltzmann_SI  = 1.3806503e-23
phy_kBoltzmann_CGS = 1.3806503e-16
phy_hPlanck_SI  = 6.62606896e-34
phy_hPlanck_CGS = 6.62606896e-27
phy_hbarPlanck_SI  = 1.054571628e-34
phy_hbarPlanck_CGS = 1.054571628e-27
phy_GravitationConst_SI  = 6.67428e-11
phy_GravitationConst_CGS = 6.67428e-8
phy_SpeedOfLight_SI  = 299792458e0
phy_SpeedOfLight_CGS = 299792458e2
phy_StefanBoltzmann_SI = 5.6704e-8
phy_StefanBoltzmann_CGS = 5.670373e-5
phy_IdealGasConst_SI = 8.314472e0
phy_ThomsonScatterCross_CGS = 6.6524574e-25

phy_Lsun_SI = 3.839e26 # J s-1
phy_Lsun_CGS = 3.839e33 # erg s-1
phy_Msun_SI = 1.9891e30 # kg
phy_Msun_CGS = 1.9891e33 # g
phy_Rsun_SI = 6.955e8 # m
phy_Rsun_CGS = 6.955e10 # cm

phy_Rearth_CGS = 6371e5 # cm
phy_Mearth_CGS = 5.97219e27 # g
phy_Mmoon_CGS = 7.34767309e25 # g

phy_SecondsPerYear = 3600e0*24e0*365e0
phy_Deg2Rad = phy_Pi/180e0
phy_erg2joule = 1e-7
phy_m2cm = 1e2
phy_kg2g = 1e3
phy_eV2erg = 1.60217657e-12
phy_cm_1_2erg = phy_hPlanck_CGS * phy_SpeedOfLight_CGS
phy_cm_1_2K = phy_cm_1_2erg/phy_kBoltzmann_CGS
phy_AvogadroConst = 6.02214179e23
phy_AU2cm = 1.49597871e13
phy_AU2m  = 1.49597871e11
phy_pc2m  = 3.08567758e16
phy_pc2cm = 3.08567758e18
phy_Angstrom2micron  = 1e-4
phy_Angstrom2cm  = 1e-8
phy_micron2cm  = 1e-4

phy_jansky2CGS = 1e-23
phy_jansky2SI  = 1e-26

phy_CMB_T = 2.72548e0

phy_ratioDust2GasMass_ISM = 0.01e0
phy_Habing_photon_energy_CGS = 1.99e-11
phy_LyAlpha_energy_CGS = 1.64e-11
phy_UV_cont_energy_CGS = phy_Habing_photon_energy_CGS
phy_Habing_energy_density_CGS = 5.29e-14 # Draine 2011 book, equation 12.6
phy_Habing_photon_flux_CGS = 6e7 # cm-2 s-1
phy_Habing_energy_flux_CGS = 1.194e-3 # erg cm-2 s-1
phy_UVext2Av = 2.6e0 # Tielens 2005, eq 3.19

phy_LyAlpha_nu0 = 2.4660718e15
phy_LyAlpha_l0 = 1215.668e0
phy_LyAlpha_dnul = 9.938e7
phy_LyAlpha_f12 = 0.4162e0

const_LyAlpha_cross_H2O = 1.2e-17 # Van Dishoeck 2006, Table 1
const_LyAlpha_cross_OH  = 1.8e-18 # Van Dishoeck 2006, Table 1

const_cosmicray_attenuate_N = 5.75e25 # 96 g cm-2, Nomura 2007
const_PAH_abundance_0 = 1.6e-7
const_SitesDensity_CGS = 1e15
#
#double precision :: colDen2Av_coeff = 1e-21 # Sun Kwok, eq 10.21
phy_colDen2Av_coeff = 5.3e-22 # Draine 2011, eq 21.7
