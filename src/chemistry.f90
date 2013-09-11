module chemistry

use trivials
use phy_const
use data_struct


implicit none


integer, parameter, public :: const_len_species_name       = 12
integer, parameter, private :: const_len_reactionfile_row   = 150
integer, parameter, public :: const_len_init_abun_file_row = 64
integer, parameter, private :: const_nSpecies_guess         = 512
integer, parameter, private :: const_n_dupli_max_guess      = 8
integer, parameter, private :: const_n_reac_max             = 3
integer, parameter, private :: const_n_prod_max             = 4
!
integer, parameter, public :: const_nElement               = 17
character(LEN=8), dimension(const_nElement), parameter :: &
  const_nameElements = &
    (/'+-      ', 'E       ', 'Grain   ', 'H       ', &
      'D       ', 'He      ', 'C       ', 'N       ', &
      'O       ', 'Si      ', 'S       ', 'Fe      ', &
      'Na      ', 'Mg      ', 'Cl      ', 'P       ', &
      'F       '/)
double precision, dimension(const_nElement), parameter :: &
  const_ElementMassNumber = &
    (/0D0,        5.45D-4,    0D0,        1D0,        &
      2D0,        4D0,        12D0,       14D0,       &
      16D0,       28D0,       32D0,       56D0,       &
      23D0,       24D0,       35.5D0,     31D0,       &
      19D0/)


type :: type_chemical_evol_idx_species
  integer i_H2, i_HI, i_E, i_CI, i_Cplus, i_OI, i_O2, i_CO, i_H2O, i_OH, i_Hplus, i_gH
  integer iiH2, iiHI, iiE, iiCI, iiCplus, iiOI, iiO2, iiCO, iiH2O, iiOH, iiHplus, iigH
  integer :: nItem = 12
  integer, dimension(:), allocatable :: idx
  character(len=8), dimension(12) :: names = &
    (/'H2      ', 'H       ', 'E-      ', 'C       ', 'C+      ', &
      'O       ', 'O2      ', 'CO      ', 'H2O     ', 'OH      ', &
      'H+      ', 'gH      '/)
end type type_chemical_evol_idx_species


type :: type_chemical_evol_reactions_str
  integer nReactions
  character :: commentChar = '!'
  character(len=const_len_reactionfile_row), dimension(:), allocatable :: list
end type type_chemical_evol_reactions_str


type :: type_chemical_evol_a_list
  integer nItem
  integer, dimension(:), allocatable :: list
end type type_chemical_evol_a_list


type :: type_chemical_evol_a_list_big
  integer nItem
  integer, dimension(:), allocatable :: list
  integer, dimension(:), allocatable :: n_repeat
  double precision, dimension(:), allocatable :: contri
end type type_chemical_evol_a_list_big


type :: type_chemical_evol_reactions
  integer nReactions
  character(len=const_len_species_name), dimension(:,:), allocatable :: &
    reac_names, prod_names
  integer, dimension(:,:), allocatable :: reac, prod
  integer, dimension(:), allocatable :: n_reac, n_prod
  double precision, dimension(:,:), allocatable :: ABC
  double precision, dimension(:,:), allocatable :: T_range
  integer, dimension(:), allocatable :: itype
  character(len=2), dimension(:), allocatable :: ctype
  character, dimension(:), allocatable :: quality
  double precision, dimension(:), allocatable :: rates
  type(type_chemical_evol_a_list), dimension(:), allocatable :: dupli
end type type_chemical_evol_reactions


type :: type_chemical_evol_species
  integer nSpecies
  character(len=const_len_species_name), dimension(:), allocatable :: names
  double precision, dimension(:), allocatable :: mass_num
  double precision, dimension(:), allocatable :: vib_freq
  double precision, dimension(:), allocatable :: Edesorb
  double precision, dimension(:), allocatable :: adsorb_coeff
  double precision, dimension(:), allocatable :: desorb_coeff
  integer, dimension(:), allocatable :: idx_gasgrain_counterpart
  integer, dimension(:,:), allocatable :: elements
  type(type_chemical_evol_a_list_big), dimension(:), allocatable :: produ, destr
end type type_chemical_evol_species


type :: type_chemical_evol_solver_params
  character(len=128) chem_files_dir, filename_chemical_network, filename_initial_abundances
  double precision RTOL, ATOL
  double precision t_max, dt_first_step, ratio_tstep
  logical allow_stop_before_t_max
  logical H2_form_use_moeq
  logical neutralize
  real :: max_runtime_allowed = 3600.0 ! seconds
  integer :: mxstep_per_interval = 2000
  integer n_record, n_record_real
  integer NEQ, ITOL, ITASK, ISTATE, IOPT, LIW, LRW, MF, NNZ
  integer NERR
end type type_chemical_evol_solver_params


type :: type_chemical_evol_solver_storage
  double precision, dimension(:), allocatable :: RWORK
  integer, dimension(:), allocatable :: IWORK
  logical, dimension(:, :), allocatable :: sparseMaskJac
  double precision, dimension(:), allocatable :: y, y0
  double precision, dimension(:), allocatable :: ydot
  double precision, dimension(:), allocatable :: touts
  double precision, dimension(:, :), allocatable :: record
end type type_chemical_evol_solver_storage


type :: type_elemental_residence
  integer ielement
  double precision, dimension(:), allocatable :: ele_frac, ele_accu
  integer, dimension(:), allocatable :: iSpecies
  integer n_nonzero
end type type_elemental_residence


type(type_chemical_evol_reactions_str)  :: chem_reac_str

type(type_chemical_evol_idx_species)    :: chem_idx_some_spe

type(type_chemical_evol_reactions)      :: chem_net

type(type_chemical_evol_species)        :: chem_species

type(type_chemical_evol_solver_params)  :: chem_solver_params
type(type_chemical_evol_solver_storage) :: chem_solver_storage

! This thing is specific to each cell.
type(type_cell_rz_phy_basic), pointer   :: chem_params => null()


type(type_elemental_residence), dimension(:), allocatable :: chem_ele_resi

integer, parameter :: chem_ele_resi_nmax = 10
double precision, parameter :: chem_ele_frac_threshold_to_sum = 0.999D0
double precision, parameter :: chem_ele_frac_threshold_to_max = 1D-5

double precision, parameter :: const_cosmicray_intensity_0 = 1.36D-17 ! UMIST paper
double precision, parameter :: CosmicDesorpPreFactor = 3.16D-19
double precision, parameter :: CosmicDesorpGrainT = 70D0
double precision, parameter :: SitesDensity_CGS = 1D15

namelist /chemistry_configure/ &
  chem_solver_params


contains


subroutine chem_load_initial_abundances
  integer fU, i, ios
  character(len=const_len_init_abun_file_row) str
  if (.NOT. getFileUnit (fU)) then
    write(*,*) 'Cannot get a file unit!  In chem_load_initial_abundances.'
    stop
  end if
  call openFileSequentialRead(fU, combine_dir_filename(chem_solver_params%chem_files_dir, &
                                  chem_solver_params%filename_initial_abundances), 999)
  chem_solver_storage%y = 0D0
  do
    read(fU, FMT='(A)', IOSTAT=ios) str
    if (ios .NE. 0) then
      exit
    end if
    do i=1, chem_species%nSpecies
      if (trim(str(1:const_len_species_name)) .EQ. chem_species%names(i)) then
        read(str(const_len_species_name+1:const_len_init_abun_file_row), &
          '(ES16.6)') chem_solver_storage%y(i)
        exit
      end if
    end do
  end do
  close(fU)
  chem_solver_storage%y0 = chem_solver_storage%y
end subroutine chem_load_initial_abundances


subroutine chem_neutralize
  double precision totalcharge
  totalcharge = sum(chem_solver_storage%y(:) * dble(chem_species%elements(1,:)))
  chem_solver_storage%y(chem_idx_some_spe%i_E) = &
    chem_solver_storage%y(chem_idx_some_spe%i_E) + totalcharge
end subroutine chem_neutralize


subroutine chem_evol_solve_prepare
  !chem_solver_params%allow_stop_before_t_max = .FALSE.
  !chem_solver_params%dt_first_step = 1D-3
  !chem_solver_params%ratio_tstep = 1.4D0
  !chem_solver_params%max_runtime_allowed = 60.0
  chem_solver_params%n_record = ceiling( &
    log(chem_solver_params%t_max / chem_solver_params%dt_first_step * &
        (chem_solver_params%ratio_tstep - 1D0) + 1D0) &
    / &
    log(chem_solver_params%ratio_tstep))
  if (.NOT. allocated(chem_solver_storage%y)) then
    allocate(&
      chem_solver_storage%y(chem_species%nSpecies), &
      chem_solver_storage%y0(chem_species%nSpecies), &
      chem_solver_storage%ydot(chem_species%nSpecies), &
      chem_solver_storage%touts(chem_solver_params%n_record), &
      chem_solver_storage%record(chem_species%nSpecies, chem_solver_params%n_record))
  end if
end subroutine chem_evol_solve_prepare


subroutine chem_set_solver_flags
  chem_solver_params%IOPT = 1 ! 1: allow optional input; 0: disallow
  chem_solver_params%MF = 021 ! Line 2557 of opkdmain.f
  chem_solver_params%NERR = 0 ! for counting number of errors in iteration
  chem_solver_params%ITOL = 1 ! scalar control of tolerance
  chem_solver_params%ITASK = 1 ! normal, allow overshoot
  chem_solver_params%ISTATE = 1 ! first call
  !chem_solver_params%RTOL = 1D-6
  !chem_solver_params%ATOL = 1D-30
end subroutine chem_set_solver_flags


subroutine chem_evol_solve
  use my_timer
  external chem_ode_f, chem_ode_jac
  integer i, itmp
  double precision t, tout, t_step, t_scale_min
  double precision :: const_factor = 1D3
  type(atimer) timer
  real time_thisstep, runtime_thisstep, time_laststep, runtime_laststep
  !--
  character(len=128) :: chem_evol_save_filename = 'chem_evol_tmp.dat'
  character(len=32) fmtstr
  logical flag_chem_evol_save
  integer fU_chem_evol_save
  !
  flag_chem_evol_save = .false.
  if (flag_chem_evol_save) then
    if (.not. getFileUnit(fU_chem_evol_save)) then
      write(*,*) 'Cannot get a unit for output!  In chem_evol_solve.'
      stop
    end if
    call openFileSequentialWrite(fU_chem_evol_save, chem_evol_save_filename, 99999)
    write(fmtstr, '("(", I4, "A14)")') chem_species%nSpecies+1
    write(fU_chem_evol_save, fmtstr) '!Time', chem_species%names(1:chem_species%nSpecies)
    write(fmtstr, '("(", I4, "ES14.4E4)")') chem_species%nSpecies+1
  end if
  !--
  t = 0D0
  tout = chem_solver_params%dt_first_step
  t_step = chem_solver_params%dt_first_step
  chem_solver_storage%touts(1) = tout
  chem_solver_storage%record(:,1) = chem_solver_storage%y
  !
  call timer%init('Chem')
  time_laststep = timer%elapsed_time()
  runtime_laststep = huge(0.0)
  !
  do i=2, chem_solver_params%n_record
    write (*, '(A, 25X, "Solving chemistry... ", I5, " (", F5.1, "%)", &
              &"  t = ", ES9.2, "  tStep = ", ES9.2)') &
      CHAR(27)//'[A', i, real(i*100)/real(chem_solver_params%n_record), t, t_step
    call DLSODES( &
         chem_ode_f, &
         !
         chem_solver_params%NEQ, &
         chem_solver_storage%y, &
         !
         t, &
         tout, &
         !
         chem_solver_params%ITOL, &
         chem_solver_params%RTOL, &
         chem_solver_params%ATOL, &
         chem_solver_params%ITASK, &
         chem_solver_params%ISTATE, &
         chem_solver_params%IOPT, &
         chem_solver_storage%RWORK, &
         chem_solver_params%LRW, &
         chem_solver_storage%IWORK, &
         chem_solver_params%LIW, &
         !
         chem_ode_jac, &
         !
         chem_solver_params%MF)
    !--
    if (flag_chem_evol_save) then
      write(fU_chem_evol_save, fmtstr) tout, chem_solver_storage%y
    end if
    !--
    time_thisstep = timer%elapsed_time()
    runtime_thisstep = time_thisstep - time_laststep
    if ((runtime_thisstep .gt. max(5.0*runtime_laststep, 0.1*chem_solver_params%max_runtime_allowed)) &
        .or. &
        (time_thisstep .gt. chem_solver_params%max_runtime_allowed)) then
      write(*, '(A, ES9.2/)') 'Premature finish: t = ', t
      exit
    end if
    time_laststep = time_thisstep
    runtime_laststep = runtime_thisstep
    !
    if (chem_solver_params%ISTATE .LT. 0) then
      chem_solver_params%NERR = chem_solver_params%NERR + 1
      write(*, '(A, I3/)') 'Error: ', chem_solver_params%ISTATE
      chem_solver_params%ISTATE = 3
    end if
    if (chem_solver_params%allow_stop_before_t_max) then
      call dintdy(t, 1, &
        chem_solver_storage%RWORK(chem_solver_storage%IWORK(22)), &
        chem_solver_params%NEQ, chem_solver_storage%ydot, itmp)
      !where (chem_solver_storage%y .LE. 1D-15)
      !  chem_solver_storage%ydot = 1D-200
      !end where
      t_scale_min = minval(abs(chem_solver_storage%y(chem_idx_some_spe%idx) / &
                     chem_solver_storage%ydot(chem_idx_some_spe%idx)))
      if (t_scale_min .GT. const_factor * chem_solver_params%t_max) then
        chem_solver_params%n_record_real = i
        write(*, '(A, ES9.2, 2X, A, ES9.2/)') 'Early finishing: t_scale_min = ', &
          t_scale_min, ' t = ', t
        exit
      end if
    end if
    chem_solver_storage%touts(i) = tout
    chem_solver_storage%record(:,i) = chem_solver_storage%y
    t_step = t_step * chem_solver_params%ratio_tstep
    tout = t + t_step
  end do
  !--
  if (flag_chem_evol_save) then
    close(fU_chem_evol_save)
  end if
  !--
end subroutine chem_evol_solve


subroutine chem_cal_rates
  integer i, j, k, i1, i2
  double precision T300!, TemperatureReduced, JNegaPosi, JChargeNeut
  double precision, dimension(4) :: tmpVecReal
  integer, dimension(1) :: tmpVecInt
  double precision :: SitesPerGrain, tmp
  T300 = chem_params%Tgas / 300D0
  ! TemperatureReduced = phy_kBoltzmann_SI * chem_params%Tgas / &
  !   (phy_elementaryCharge_SI**2 * phy_CoulombConst_SI / &
  !   (chem_params%GrainRadius_CGS*1D-4))
  ! Pagani 2009, equation 11, 12, 13; not activated here.
  ! JNegaPosi = (1D0 + 1D0/TemperatureReduced) * &
  !             (1D0 + sqrt(2D0/(2D0+TemperatureReduced)))
  ! JChargeNeut = (1D0 + sqrt(phy_Pi/2D0/TemperatureReduced))
  !
  SitesPerGrain = SitesDensity_CGS * (4D0 * phy_Pi * chem_params%GrainRadius_CGS**2)
  !
  ! if (chem_params%ratioDust2HnucNum .LT. 1D-80) then ! If not set, calculate it.
  !   chem_params%ratioDust2HnucNum = & ! n_Grain/n_H
  !     chem_params%ratioDust2GasMass * (phy_mProton_CGS * chem_params%MeanMolWeight) &
  !     / (4.0D0*phy_Pi/3.0D0 * (chem_params%GrainRadius_CGS)**3 * &
  !        chem_params%GrainMaterialDensity_CGS)
  ! end if
  ! Le Petit 2009, equation 46 (not quite clear)
  ! Le Bourlot 1995, Appendix A
  ! Formation rate of H2:
  !   d/dt n(H2) = chem_params%R_H2_form_rate * n(H) * n_gas
  !if (chem_params%R_H2_form_rate .LT. 1D-80) then ! If not set, calculate it.
  !  ! Le Petit 2006
  !  chem_params%stickCoeffH = sqrt(10D0/max(10D0, chem_params%Tgas))
  !  chem_params%R_H2_form_rate = 0.5D0 * chem_params%stickCoeffH &
  !    * 3D0 / 4D0 * chem_params%MeanMolWeight * phy_mProton_CGS &
  !    * chem_params%ratioDust2GasMass / chem_params%GrainMaterialDensity_CGS &
  !    / sqrt(chem_params%aGrainMin_CGS * chem_params%aGrainMax_CGS) &
  !    * sqrt(8D0*phy_kBoltzmann_CGS*chem_params%Tgas/(phy_Pi * phy_mProton_CGS))
  !end if
  !
  do i=1, chem_net%nReactions
    ! Set the default value.
    chem_net%rates(i) = 0D0
    ! Reactions with very negative barriers AND with a temperature range not
    ! quite applicable are discarded.
    ! To check the most up-to-date UMIST document to see if there is any
    ! more detailed prescription for this.
    !
    ! Reactions with negative barriers are evil!
    !if (chem_net%ABC(3, i) .LT. -100D0) then
    !  if ((minval(chem_net%T_range(:,i)/chem_params%Tgas) .GE. 2D0) .OR. &
    !      (minval(chem_params%Tgas/chem_net%T_range(:,i)) .GE. 2D0)) then
    !    cycle
    !  end if
    !endif
    if (chem_net%ABC(3, i) .LT. 0D0) then
      if ((chem_net%T_range(1, i) .gt. chem_params%Tgas) .or. &
          (chem_net%T_range(2, i) .lt. chem_params%Tgas)) then
        cycle
      end if
    end if
    select case (chem_net%itype(i))
      case (5, 53) !- Reactions with itype=53 need not be included.
        chem_net%rates(i) = chem_net%ABC(1, i) * (T300**chem_net%ABC(2, i)) &
            * exp(-chem_net%ABC(3, i)/chem_params%Tgas)
      case (1)
        chem_net%rates(i) = (chem_params%zeta_cosmicray_H2/const_cosmicray_intensity_0) &
            * exp(-chem_params%Ncol / const_cosmicray_attenuate_N) &
            * chem_net%ABC(1, i)
      case (2) ! Todo
        chem_net%rates(i) = (chem_params%zeta_cosmicray_H2/const_cosmicray_intensity_0) &
            * exp(-chem_params%Ncol / const_cosmicray_attenuate_N) &
            * chem_net%ABC(1, i) * (T300**chem_net%ABC(2, i)) &
            * chem_net%ABC(3, i) / (1D0 - chem_params%omega_albedo)
      case (3) ! Todo
        chem_net%rates(i) = chem_params%UV_G0_factor &
            * chem_net%ABC(1, i) &
            * exp(-chem_net%ABC(3, i) * chem_params%Av) &
            * f_selfshielding(i)
      case (13)
        chem_net%rates(i) = chem_params%LymanAlpha_number_flux_0 &
            * chem_net%ABC(1, i) &
            * exp(-chem_net%ABC(3, i) * chem_params%Av) &
            * f_selfshielding(i)
      !case (0)
      !  chem_net%rates(i) = chem_net%ABC(1, i) * chem_params%R_H2_form_rate
      !
      ! itype > 60 are for grain chemistry
      !
      case (61) ! Adsorption
        ! dt(n(A)) = -sigma * v * n(A) * n(dust)
        ! dt(X(A)) = -sigma * v * X(A) * n(dust)
        ! adsorb_coeff = sigma * v * n(dust)
        chem_net%rates(i) = &
          chem_net%ABC(1, i) * phy_Pi * chem_params%GrainRadius_CGS**2 &
          * sqrt(8D0/phy_Pi*phy_kBoltzmann_CGS*chem_params%Tgas &
                 / (chem_species%mass_num(chem_net%reac(1, i)) * phy_mProton_CGS)) &
          * chem_params%n_dust
        chem_species%adsorb_coeff(chem_net%reac(1, i)) = chem_net%rates(i)
      case (62) ! Desorption
        ! dt(N(gA)) = -desorb_coeff * N(gA)
        ! dt(X(gA)) = -desorb_coeff * X(gA)
        ! Cosmic ray desorption rate from Hasegawa1993.
        ! <timestamp>2013-08-28 Wed 12:47:30</timestamp>
        ! Error corrected: the factor vib_freq is missing for the cosmic-ray contribution
        chem_net%rates(i) = &
          chem_species%vib_freq(chem_net%reac(1, i)) &
          * (exp(-chem_net%ABC(3, i)/chem_params%Tdust) &
             + CosmicDesorpPreFactor * &
               (chem_params%zeta_cosmicray_H2/const_cosmicray_intensity_0) &
                * exp(-chem_params%Ncol / const_cosmicray_attenuate_N) &
                * exp(-chem_net%ABC(3, i)/CosmicDesorpGrainT))
        chem_species%desorb_coeff(chem_net%reac(1, i)) = chem_net%rates(i)
      case (63) ! A + A -> B
        ! Moment equation:
        ! dt(A) = -k_AA * <A(A-1)>
        ! k_AA = k_mig(A) / N_site
        !
        ! dt(X(B)) = -k_AA / D2G * X(A)**2
        !
        ! dt(N(B)) = -k_AA / (k_AA + k_des) * sigma * v * n(a) * N(A)
        ! dt(X(B)) = -k_AA / (k_AA + k_des) * sigma * v * n(a) * X(A)
        !          = -k_AA / (k_AA + k_des) * sigma * v * X(a) * X(A) * n_gas
        !          = -k_AA / (k_AA + k_des) * sigma * v * X(a) * X(A) * n_dust / D2G
        !          = -k_AA * X(a) * X(A) * [sigma * v * n_dust] / (k_AA + k_des) / D2G
        i1 = chem_net%reac(1, i)
        tmp = getMobility(chem_species%vib_freq(i1), &
                          chem_species%mass_num(i1), &
                          chem_species%Edesorb(i1), &
                          chem_params%Tdust) / SitesPerGrain
        if (chem_solver_params%H2_form_use_moeq) then
          i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
          chem_net%rates(i) = &
            tmp / (tmp + chem_species%desorb_coeff(chem_net%reac(1, i))) * &
            chem_species%adsorb_coeff(i1) / chem_params%ratioDust2HnucNum
        else
          chem_net%rates(i) = tmp / chem_params%ratioDust2HnucNum
        end if
        ! Todo
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          chem_params%R_H2_form_rate_coeff = chem_net%rates(i)
        end if
      case (64) ! A + B -> xxx
        ! dt(A) = k_AB * <A.B>
        ! k_AB = (k_mig(A) + k_mig(B)) / N_site
        ! X(A) = N(A) * D2G; D2G: dust-to-grain-number-ratio
        ! dt(X(A)) = k_AB / D2G * X(A) * X(B)
        i1 = chem_net%reac(1, i)
        i2 = chem_net%reac(2, i)
        chem_net%rates(i) = ( &
          getMobility(chem_species%vib_freq(i1), &
                      chem_species%mass_num(i1), &
                      chem_species%Edesorb(i1), &
                      chem_params%Tdust) &
          + &
          getMobility(chem_species%vib_freq(i2), &
                      chem_species%mass_num(i2), &
                      chem_species%Edesorb(i2), &
                      chem_params%Tdust)) &
          / (SitesPerGrain * chem_params%ratioDust2HnucNum)
      case default
        chem_net%rates(i) = 0D0
    end select
    ! Change the time unit from seconds into years.
    chem_net%rates(i) = chem_net%rates(i) * phy_SecondsPerYear
    ! dn/dt = k n1 n2 => dx/dt := d(n/n_H)/dt = k*n_H x1 x2
    if ((chem_net%n_reac(i) .EQ. 2) .and. (chem_net%itype(i) .lt. 60)) then
      chem_net%rates(i) = chem_net%rates(i) * chem_params%n_gas
    end if
    ! Choose the reaction with temperature range that
    ! best matches the current temperature.
    ! Since different locations can have very different temperatures,
    ! it is important to do this for each location.
    if (chem_net%dupli(i)%nItem .GT. 0) then
      do j=1, chem_net%dupli(i)%nItem
        k = chem_net%dupli(i)%list(j)
        tmpVecReal(1:2) = abs(chem_net%T_range(:,k) - chem_params%Tgas)
        tmpVecReal(3:4) = abs(chem_net%T_range(:,i) - chem_params%Tgas)
        tmpVecInt = minloc(tmpVecReal)
        i1 = tmpVecInt(1)
        if ((i1 .EQ. 1) .OR. (i1 .EQ. 2)) then
          chem_net%rates(i) = 0D0
          exit
        end if
        if ((i1 .EQ. 3) .OR. (i1 .EQ. 4)) then
          chem_net%rates(k) = 0D0
          cycle
        end if
      end do
    end if
  end do
end subroutine chem_cal_rates


function f_selfshielding(iReac)
  double precision f_selfshielding
  integer iReac
  if ((chem_net%ctype(iReac) .NE. 'PH') .OR. &
      (chem_net%ctype(iReac) .NE. 'LA')) then
    f_selfshielding = 1D0
    return
  end if
  select case (chem_species%names(chem_net%reac(1, iReac)))
    case ('H2')
      f_selfshielding = chem_params%f_selfshielding_H2
    case ('H2O')
      f_selfshielding = chem_params%f_selfshielding_H2O
    case ('OH')
      f_selfshielding = chem_params%f_selfshielding_OH
    case ('CO')
      f_selfshielding = chem_params%f_selfshielding_CO
    case default
      f_selfshielding = 1D0
  end select
end function f_selfshielding


subroutine chem_prepare_solver_storage
  integer i, j, k
  chem_solver_params%LRW = &
    20 + 4 * chem_solver_params%NNZ + 28 * chem_solver_params%NEQ
    !20 + chem_solver_params%NEQ * (12 + 1) &
    !+ 3 * chem_solver_params%NEQ + 4 * chem_solver_params%NNZ &
    !+ 2 * chem_solver_params%NEQ + chem_solver_params%NNZ &
    !+ 10 * chem_solver_params%NEQ
  chem_solver_params%LIW = 31 + chem_solver_params%NEQ + chem_solver_params%NNZ
  allocate( &
    chem_solver_storage%RWORK(chem_solver_params%LRW), &
    chem_solver_storage%IWORK(chem_solver_params%LIW))
  chem_solver_storage%RWORK(5:10) = 0D0
  chem_solver_storage%IWORK(5:10) = 0
  chem_solver_storage%IWORK(6) = chem_solver_params%mxstep_per_interval
  chem_solver_storage%IWORK(31) = 1
  k = 1
  do i=1, chem_solver_params%NEQ
    do j=1, chem_solver_params%NEQ
      if (chem_solver_storage%sparseMaskJac(j, i)) then
        chem_solver_storage%IWORK(31 + chem_solver_params%NEQ + k) = j
        k = k + 1
      end if
    end do
    chem_solver_storage%IWORK(31+i) = k
  end do
  deallocate(chem_solver_storage%sparseMaskJac)
end subroutine chem_prepare_solver_storage


subroutine chem_make_sparse_structure
  integer i, j, k
  chem_solver_params%NEQ = chem_species%nSpecies
  allocate(chem_solver_storage%sparseMaskJac(chem_species%nSpecies, &
    chem_species%nSpecies))
  chem_solver_storage%sparseMaskJac = .FALSE.
  do i=1, chem_net%nReactions
    do j=1, chem_net%n_reac(i)
      do k=1, chem_net%n_reac(i)
        chem_solver_storage%sparseMaskJac &
          (chem_net%reac(k, i), chem_net%reac(j, i)) = .TRUE.
      end do
      do k=1, chem_net%n_prod(i)
        chem_solver_storage%sparseMaskJac &
          (chem_net%prod(k, i), chem_net%reac(j, i)) = .TRUE.
      end do
    end do
  end do
  chem_solver_params%NNZ = count(chem_solver_storage%sparseMaskJac)
end subroutine chem_make_sparse_structure


subroutine chem_get_idx_for_special_species
  integer i
  allocate(chem_idx_some_spe%idx(chem_idx_some_spe%nItem))
  do i=1, chem_species%nSpecies
    select case (trim(chem_species%names(i)))
      case ('H2')
        chem_idx_some_spe%i_H2 = i
        chem_idx_some_spe%iiH2 = 1
        chem_idx_some_spe%idx(1) = i
      case ('H')
        chem_idx_some_spe%i_HI = i
        chem_idx_some_spe%iiHI = 2
        chem_idx_some_spe%idx(2) = i
      case ('E-')
        chem_idx_some_spe%i_E = i
        chem_idx_some_spe%iiE = 3
        chem_idx_some_spe%idx(3) = i
      case ('C')
        chem_idx_some_spe%i_CI = i
        chem_idx_some_spe%iiCI = 4
        chem_idx_some_spe%idx(4) = i
      case ('C+')
        chem_idx_some_spe%i_Cplus = i
        chem_idx_some_spe%iiCplus = 5
        chem_idx_some_spe%idx(5) = i
      case ('O')
        chem_idx_some_spe%i_OI = i
        chem_idx_some_spe%iiOI = 6
        chem_idx_some_spe%idx(6) = i
      case ('O2')
        chem_idx_some_spe%i_O2 = i
        chem_idx_some_spe%iiO2 = 7
        chem_idx_some_spe%idx(7) = i
      case ('CO')
        chem_idx_some_spe%i_CO = i
        chem_idx_some_spe%iiCO = 8
        chem_idx_some_spe%idx(8) = i
      case ('H2O')
        chem_idx_some_spe%i_H2O = i
        chem_idx_some_spe%iiH2O = 9
        chem_idx_some_spe%idx(9) = i
      case ('OH')
        chem_idx_some_spe%i_OH = i
        chem_idx_some_spe%iiOH = 10
        chem_idx_some_spe%idx(10) = i
      case ('H+')
        chem_idx_some_spe%i_Hplus = i
        chem_idx_some_spe%iiHplus = 11
        chem_idx_some_spe%idx(11) = i
      case ('gH')
        chem_idx_some_spe%i_gH = i
        chem_idx_some_spe%iigH = 12
        chem_idx_some_spe%idx(12) = i
    end select
  end do
end subroutine chem_get_idx_for_special_species


subroutine chem_get_dupli_reactions
  ! Find out the duplicated reactions.
  ! Among all the reaction belonging to a single "duplicated set",
  ! one and only one will be used.
  ! For each reaction, this subroutine finds out all the reactions with smaller
  ! indices and the same reactants, products, and reaction types.
  integer i, j, i1
  integer, dimension(const_n_dupli_max_guess) :: indices
  do i=1, chem_net%nReactions
    chem_net%dupli(i)%nItem = 0
    i1 = 0
    do j=1, i-1
      if ((chem_net%ctype(i) .EQ. chem_net%ctype(j)) .AND. &
          (chem_net%itype(i) .EQ. chem_net%itype(j)) .AND. &
          (sum(abs(chem_net%reac(:, j)-chem_net%reac(:, i))) .EQ. 0) .AND. &
          (sum(abs(chem_net%prod(:, j)-chem_net%prod(:, i))) .EQ. 0)) then
        i1 = i1 + 1
        indices(i1) = j
      end if
    end do
    if (i1 .GT. 0) then
      chem_net%dupli(i)%nItem = i1
      allocate(chem_net%dupli(i)%list(i1))
      chem_net%dupli(i)%list = indices(1:i1) 
    end if
  end do
end subroutine chem_get_dupli_reactions


subroutine chem_parse_reactions
  integer i, j, k, n_tmp
  logical flag
  character(len=const_len_species_name), &
    dimension(const_nSpecies_guess) :: names_tmp
  chem_net%reac = 0
  chem_net%prod = 0
  names_tmp(1) = chem_net%reac_names(1, 1)
  n_tmp = 1
  do i=1, chem_net%nReactions
    do k=1, chem_net%n_reac(i)
      flag = .TRUE.
      do j=1, n_tmp
        if (trim(names_tmp(j)) .EQ. trim(chem_net%reac_names(k, i))) then
          flag = .FALSE.
          chem_net%reac(k, i) = j
          exit
        end if
      end do
      if (flag) then
        n_tmp = n_tmp + 1
        names_tmp(n_tmp) = chem_net%reac_names(k, i)
        chem_net%reac(k, i) = n_tmp
      end if
    end do
    do k=1, chem_net%n_prod(i)
      flag = .TRUE.
      do j=1, n_tmp
        if (trim(names_tmp(j)) .EQ. trim(chem_net%prod_names(k, i))) then
          flag = .FALSE.
          chem_net%prod(k, i) = j
          exit
        end if
      end do
      if (flag) then
        n_tmp = n_tmp + 1
        names_tmp(n_tmp) = chem_net%prod_names(k, i)
        chem_net%prod(k, i) = n_tmp
      end if
    end do
  end do
  !
  chem_species%nSpecies = n_tmp
  allocate( &
        chem_species%names(n_tmp), &
        chem_species%mass_num(n_tmp), &
        chem_species%vib_freq(n_tmp), &
        chem_species%Edesorb(n_tmp), &
        chem_species%adsorb_coeff(n_tmp), &
        chem_species%desorb_coeff(n_tmp), &
        chem_species%idx_gasgrain_counterpart(n_tmp), &
        chem_species%elements(const_nElement, n_tmp))
  !
  chem_species%names = names_tmp(1:n_tmp)
  ! Initialize to invalid values so that error will be explicitly shown (just in case).
  chem_species%vib_freq = dblNaN()
  chem_species%Edesorb = dblNaN()
  chem_species%adsorb_coeff = dblNaN()
  chem_species%desorb_coeff = dblNaN()
  chem_species%idx_gasgrain_counterpart = -1
  do i=1, chem_species%nSpecies
    call getElements(chem_species%names(i), const_nameElements, &
                     const_nElement, chem_species%elements(:, i))
    chem_species%mass_num(i) = sum(dble(chem_species%elements(:, i)) * &
                                   const_ElementMassNumber)
  end do
  !
  ! deallocate(chem_species%elements) ! Not used later.
  !
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (62)
        chem_species%vib_freq(chem_net%reac(1, i)) = &
          getVibFreq(chem_species%mass_num(chem_net%reac(1, i)), &
                     chem_net%ABC(3, i))
        chem_species%Edesorb(chem_net%reac(1, i)) = chem_net%ABC(3, i)
        chem_species%idx_gasgrain_counterpart(chem_net%prod(1, i)) = chem_net%reac(1, i)
        chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i)) = chem_net%prod(1, i)
      case default
        cycle
    end select
  end do
end subroutine chem_parse_reactions


subroutine chem_load_reactions
  integer i, j, k, ios
  chem_net%nReactions = chem_reac_str%nReactions
  allocate(chem_net%reac_names(const_n_reac_max, chem_net%nReactions), &
           chem_net%prod_names(const_n_prod_max, chem_net%nReactions), &
           chem_net%reac(const_n_reac_max, chem_net%nReactions), &
           chem_net%prod(const_n_prod_max, chem_net%nReactions), &
           chem_net%n_reac(chem_net%nReactions), &
           chem_net%n_prod(chem_net%nReactions), &
           chem_net%ABC(3, chem_net%nReactions), &
           chem_net%T_range(2, chem_net%nReactions), &
           chem_net%itype(chem_net%nReactions), &
           chem_net%ctype(chem_net%nReactions), &
           chem_net%quality(chem_net%nReactions), &
           chem_net%rates(chem_net%nReactions), &
           chem_net%dupli(chem_net%nReactions))
  chem_net%reac_names = ' '
  chem_net%prod_names = ' '
  chem_net%n_reac = 0
  chem_net%n_prod = 0
  do i=1, chem_net%nReactions
    read(chem_reac_str%list(i), FMT = &
      '(7(A12), 3F9.0, 2F6.0, I3, X, A1, X, A2)', IOSTAT=ios) &
      chem_net%reac_names(:,i), &
      chem_net%prod_names(:,i), &
      chem_net%ABC(:,i), &
      chem_net%T_range(:,i), &
      chem_net%itype(i), &
      chem_net%quality(i), &
      chem_net%ctype(i)
    do j=1, const_n_reac_max
      do k=1, const_len_species_name
        if (chem_net%reac_names(j, i)(k:k) .NE. ' ') then
          chem_net%n_reac(i) = chem_net%n_reac(i) + 1
          exit
        end if
      end do
      if (trim(chem_net%reac_names(j, i)) .EQ. 'PHOTON') then
        chem_net%n_reac(i) = chem_net%n_reac(i) - 1
      end if
      if (trim(chem_net%reac_names(j, i)) .EQ. 'CRPHOT') then
        chem_net%n_reac(i) = chem_net%n_reac(i) - 1
      end if
      if (trim(chem_net%reac_names(j, i)) .EQ. 'CRP') then
        chem_net%n_reac(i) = chem_net%n_reac(i) - 1
      end if
    end do
    do j=1, const_n_prod_max
      do k=1, const_len_species_name
        if (chem_net%prod_names(j, i)(k:k) .NE. ' ') then
          chem_net%n_prod(i) = chem_net%n_prod(i) + 1
          exit
        end if
      end do
      if (trim(chem_net%prod_names(j, i)) .EQ. 'PHOTON') then
        chem_net%n_prod(i) = chem_net%n_prod(i) - 1
      end if
    end do
  end do
end subroutine chem_load_reactions


subroutine chem_read_reactions()
  integer fU, i, ios
  chem_reac_str%nReactions = &
    GetFileLen_comment_blank(combine_dir_filename(chem_solver_params%chem_files_dir, &
                             chem_solver_params%filename_chemical_network), &
    chem_reac_str%commentChar)
  if (.NOT. getFileUnit (fU)) then
    write(*,'(/A/)') 'Cannot get a file unit!  In chem_read_reactions.'
    stop
  end if
  call openFileSequentialRead(fU, combine_dir_filename(chem_solver_params%chem_files_dir, &
                                  chem_solver_params%filename_chemical_network), 999)
  allocate(chem_reac_str%list(chem_reac_str%nReactions))
  i = 1
  do
    read(fU, FMT='(A)', IOSTAT=ios) chem_reac_str%list(i)
    if (ios .NE. 0) then
      exit
    end if
    if (.NOT. &
        ( &
         (chem_reac_str%list(i)(1:1) .EQ. chem_reac_str%commentChar) .OR. &
         (chem_reac_str%list(i)(1:1) .EQ. ' ') &
        )) then
      i = i + 1
    end if
    if (i .GT. chem_reac_str%nReactions) then
      exit
    end if
  end do
  close(fU)
end subroutine chem_read_reactions


! Get the elemental composition of each molecule.
subroutine getElements &
  (nameSpec, listElements, nElements, arrNElements)
character(len=*) nameSpec, listElements(nElements)
integer, dimension(nElements) :: arrNElements
integer i, j, k, ntmp, lenName, lenEle, nElements
integer, dimension(32) :: belongto
logical, dimension(32) :: used
logical flagReplace
integer, parameter :: chargePos = 1
arrNElements = 0
lenName = len(trim(nameSpec))
belongto = 0
used = .FALSE.
do i=1, nElements
  lenEle = len(trim(listElements(i)))
  do j=1, lenName-lenEle+1
    if (nameSpec(j:(j+lenEle-1)) .EQ. &
        listElements(i)(1:(lenEle))) then
      flagReplace = .TRUE.
      do k=j, (j+lenEle-1)
        if (used(k)) then
          if (len(trim(listElements(belongto(k)))) .GE. &
              len(trim(listElements(i)))) then
            flagReplace = .FALSE.
            exit
          else
            arrNElements(belongto(k)) = &
              arrNElements(belongto(k)) - 1
          end if
        end if
      end do
      if (flagReplace) then
        belongto(j:(j+lenEle-1)) = i
        used(j:(j+lenEle-1)) = .TRUE.
        arrNElements(i) = arrNElements(i) + 1
      end if
    end if
  end do
end do
!
do i=2, lenName
  if (.NOT. used(i)) then
    do j=1, (i-1)
      if (used(i-j)) then
        belongto(i) = belongto(i-j)
        exit
      end if
    end do
    if (((nameSpec(i-1:i-1) .GT. '9') .OR. &
      (nameSpec(i-1:i-1) .LT. '0')) .AND. &
      (nameSpec(i:i) .LE. '9') .AND. &
      (nameSpec(i:i) .GE. '0')) then
      if ((nameSpec(i+1:i+1) .LE. '9') .AND. &
        (nameSpec(i+1:i+1) .GE. '0')) then
        read (nameSpec(i:i+1), '(I2)') ntmp
        if (ntmp .EQ. 0) cycle
        arrNElements(belongto(i)) = &
          arrNElements(belongto(i)) + ntmp - 1
      else
        read (nameSpec(i:i), '(I1)') ntmp
        if (ntmp .EQ. 0) cycle
        arrNElements(belongto(i)) = &
          arrNElements(belongto(i)) + ntmp - 1
      end if
    else if (nameSpec(i:i) .EQ. '+') then
      arrNElements(chargePos) = 1
    else if (nameSpec(i:i) .EQ. '-') then
      arrNElements(chargePos) = -1
    end if
  end if
end do
end subroutine getElements


function getVibFreq(massnum, Edesorb)
double precision getVibFreq
double precision, intent(in) :: massnum, Edesorb
getVibFreq = &
   sqrt(2D0 * SitesDensity_CGS * &
      phy_kBoltzmann_CGS * Edesorb / (phy_Pi**2) / &
      (phy_mProton_CGS * massnum))
end function getVibFreq


function getMobility(vibfreq, massnum, Edesorb, Tdust)
  double precision getMobility
  double precision, intent(in) :: vibfreq, massnum, Edesorb, Tdust
  double precision, parameter :: Diff2DesorRatio = 0.5D0
  double precision, parameter :: DiffBarrierWidth_CGS = 1D-8
  ! Edesorb is in Kelvin.
  ! 2a/hbar * sqrt(2mE)
  getMobility = vibfreq * exp(max( &
              -Edesorb * Diff2DesorRatio / Tdust, &
              !
              -2D0 * DiffBarrierWidth_CGS / phy_hbarPlanck_CGS * &
                sqrt(2D0 * massnum * phy_mProton_CGS &
                  * phy_kBoltzmann_CGS * Edesorb * Diff2DesorRatio)))
end function getMobility


subroutine chem_elemental_residence
  use quick_sort
  double precision, dimension(:,:), allocatable :: ele_spe, tmp
  integer i, j, i0
  double precision accum, accum_total
  allocate(ele_spe(const_nElement, chem_species%nSpecies), &
           tmp(2, chem_species%nSpecies))
  if (.not. allocated(chem_ele_resi)) then
    allocate(chem_ele_resi(const_nElement))
  end if
  do i=1, chem_species%nSpecies
    ele_spe(:, i) = chem_solver_storage%y(i) * dble(chem_species%elements(:,i))
  end do
  do i=1, const_nElement
    if (.not. allocated(chem_ele_resi(i)%ele_frac)) then
      allocate(chem_ele_resi(i)%ele_frac(chem_ele_resi_nmax), &
              chem_ele_resi(i)%ele_accu(chem_ele_resi_nmax), &
              chem_ele_resi(i)%iSpecies(chem_ele_resi_nmax))
    end if
    do j=1, chem_species%nSpecies
      tmp(1, j) = dble(j)
    end do
    tmp(2, :) = -abs(ele_spe(i, :))
    call quick_sort_array(tmp, 2, chem_species%nSpecies, 1, (/2/))
    accum = 0D0
    accum_total = sum(abs(ele_spe(i, :)))
    chem_ele_resi(i)%n_nonzero = chem_ele_resi_nmax
    do j=1, chem_ele_resi_nmax
      i0 = int(tmp(1, j))
      accum = accum + abs(ele_spe(i, i0))
      chem_ele_resi(i)%ele_frac(j) = ele_spe(i, i0) / accum_total
      chem_ele_resi(i)%ele_accu(j) = accum / accum_total
      chem_ele_resi(i)%iSpecies(j) = i0
      if ((accum .ge. chem_ele_frac_threshold_to_sum * accum_total) .or. &
          (abs(chem_ele_resi(i)%ele_frac(j)) .le. &
           chem_ele_frac_threshold_to_max * chem_ele_resi(i)%ele_frac(1))) then
        chem_ele_resi(i)%n_nonzero = j
        exit
      end if
    end do
  end do
  deallocate(ele_spe, tmp)
end subroutine chem_elemental_residence


subroutine get_species_produ_destr
  integer i, j, k, i0
  logical flag_repeat
  integer, dimension(:), allocatable :: counter1, counter2
  allocate(chem_species%produ(chem_species%nSpecies), &
           chem_species%destr(chem_species%nSpecies))
  chem_species%produ%nItem = 0
  chem_species%destr%nItem = 0
  do i=1, chem_net%nReactions
    do j=1, chem_net%n_reac(i)
      flag_repeat = .false.
      do k=1, j-1
        if (chem_net%reac(k, i) .eq. chem_net%reac(j, i)) then
          flag_repeat = .true.
          exit
        end if
      end do
      if (flag_repeat) then
        cycle
      end if
      i0 = chem_net%reac(j, i)
      chem_species%destr(i0)%nItem = chem_species%destr(i0)%nItem + 1
    end do
    do j=1, chem_net%n_prod(i)
      flag_repeat = .false.
      do k=1, j-1
        if (chem_net%prod(k, i) .eq. chem_net%prod(j, i)) then
          flag_repeat = .true.
          exit
        end if
      end do
      if (flag_repeat) then
        cycle
      end if
      i0 = chem_net%prod(j, i)
      chem_species%produ(i0)%nItem = chem_species%produ(i0)%nItem + 1
    end do
  end do
  do i=1, chem_species%nSpecies
    allocate(chem_species%produ(i)%list(chem_species%produ(i)%nItem), &
             chem_species%destr(i)%list(chem_species%destr(i)%nItem), &
             chem_species%produ(i)%n_repeat(chem_species%produ(i)%nItem), &
             chem_species%destr(i)%n_repeat(chem_species%destr(i)%nItem), &
             chem_species%produ(i)%contri(chem_species%produ(i)%nItem), &
             chem_species%destr(i)%contri(chem_species%destr(i)%nItem))
  end do
  allocate(counter1(chem_species%nSpecies), counter2(chem_species%nSpecies))
  counter1 = 0
  counter2 = 0
  do i=1, chem_net%nReactions
    do j=1, chem_net%n_reac(i)
      flag_repeat = .false.
      do k=1, j-1
        if (chem_net%reac(k, i) .eq. chem_net%reac(j, i)) then
          flag_repeat = .true.
          exit
        end if
      end do
      if (flag_repeat) then
        cycle
      end if
      i0 = chem_net%reac(j, i)
      counter1(i0) = counter1(i0) + 1
      chem_species%destr(i0)%list(counter1(i0)) = i
    end do
    do j=1, chem_net%n_prod(i)
      flag_repeat = .false.
      do k=1, j-1
        if (chem_net%prod(k, i) .eq. chem_net%prod(j, i)) then
          flag_repeat = .true.
          exit
        end if
      end do
      if (flag_repeat) then
        cycle
      end if
      i0 = chem_net%prod(j, i)
      counter2(i0) = counter2(i0) + 1
      chem_species%produ(i0)%list(counter2(i0)) = i
    end do
  end do
  do i=1, chem_species%nSpecies
    chem_species%produ(i)%n_repeat = 0
    chem_species%destr(i)%n_repeat = 0
    do j=1, chem_species%produ(i)%nItem
      i0 = chem_species%produ(i)%list(j)
      do k=1, chem_net%n_prod(i0)
        if (chem_net%prod(k, i0) .eq. i) then
          chem_species%produ(i)%n_repeat(j) = &
            chem_species%produ(i)%n_repeat(j) + 1
        end if
      end do
    end do
    do j=1, chem_species%destr(i)%nItem
      i0 = chem_species%destr(i)%list(j)
      do k=1, chem_net%n_reac(i0)
        if (chem_net%reac(k, i0) .eq. i) then
          chem_species%destr(i)%n_repeat(j) = &
            chem_species%destr(i)%n_repeat(j) + 1
        end if
      end do
    end do
  end do
  deallocate(counter1, counter2)
end subroutine get_species_produ_destr


subroutine get_contribution_each
  use quick_sort
  integer i, j, ireac
  double precision, dimension(:), allocatable :: rates_all
  double precision, dimension(:,:), allocatable :: tmp
  allocate(rates_all(chem_net%nReactions), &
           tmp(2, chem_net%nReactions))
  call chem_ode_f_alt(chem_net%nReactions, rates_all, &
                      chem_species%nSpecies, chem_solver_storage%y)
  do i=1, chem_species%nSpecies
    do j=1, chem_species%produ(i)%nItem
      ireac = chem_species%produ(i)%list(j)
      chem_species%produ(i)%contri(j) = &
        dble(chem_species%produ(i)%n_repeat(j)) * rates_all(ireac)
      tmp(1, j) = dble(ireac)
    end do
    tmp(2, 1:chem_species%produ(i)%nItem) = -chem_species%produ(i)%contri
    call quick_sort_array(tmp(:, 1:chem_species%produ(i)%nItem), &
        2, chem_species%produ(i)%nItem, 1, (/2/))
    chem_species%produ(i)%contri = -tmp(2, 1:chem_species%produ(i)%nItem)
    chem_species%produ(i)%list   = tmp(1, 1:chem_species%produ(i)%nItem)
    do j=1, chem_species%destr(i)%nItem
      ireac = chem_species%destr(i)%list(j)
      chem_species%destr(i)%contri(j) = &
        dble(chem_species%destr(i)%n_repeat(j)) * rates_all(ireac)
      tmp(1, j) = dble(ireac)
    end do
    tmp(2, 1:chem_species%destr(i)%nItem) = -chem_species%destr(i)%contri
    call quick_sort_array(tmp(:, 1:chem_species%destr(i)%nItem), &
        2, chem_species%destr(i)%nItem, 1, (/2/))
    chem_species%destr(i)%contri = -tmp(2, 1:chem_species%destr(i)%nItem)
    chem_species%destr(i)%list   = tmp(1, 1:chem_species%destr(i)%nItem)
  end do
  deallocate(rates_all, tmp)
end subroutine get_contribution_each


subroutine chem_ode_f_alt(nr, r, ny, y)
  double precision, dimension(nr), intent(out) :: r
  double precision, dimension(ny), intent(in) :: y
  integer, intent(in) :: nr, ny
  integer i, i1
  r = 0D0
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 64) ! A + B -> C ! 53
        r(i) = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(2, i))
      case (1, 2, 3, 13, 61, 62, 0) ! A -> B
        r(i) = chem_net%rates(i) * y(chem_net%reac(1, i))
      case (63) ! gA + gA -> gB
        ! dt(N(H2)) = k_HH * <H(H-1)>
        ! Moment equation:
        ! dt(N(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * N(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * X(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * X(H) * X(gH) * n_dust / D2G
        ! Rate equation:
        ! dt(X(H2)) = k_HH * X(H)**2 / D2G
        if (chem_solver_params%H2_form_use_moeq) then
          i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
          r(i) = chem_net%rates(i) * y(i1) * y(chem_net%reac(1, i))
          !ydot(i1) = ydot(i1) - rtmp ! It's like H + gH -> gH2. So dt(H) -= rtmp, dt(gH) += rtmp
          !ydot(chem_net%reac(1, i)) = ydot(chem_net%reac(1, i)) + rtmp
        else
          r(i) = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
        end if
      case default
        cycle
    end select
  end do
end subroutine chem_ode_f_alt



end module chemistry



subroutine chem_ode_f(NEQ, t, y, ydot)
  use chemistry
  implicit none
  integer NEQ, i, j, i1
  double precision t, y(NEQ), ydot(NEQ), rtmp
  ydot = 0D0
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 64) ! A + B -> C ! 53
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(2, i))
      case (1, 2, 3, 13, 61, 62, 0) ! A -> B
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case (63) ! gA + gA -> gB
        ! dt(N(H2)) = k_HH * <H(H-1)>
        ! Moment equation:
        ! dt(N(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * N(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * X(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * X(H) * X(gH) * n_dust / D2G
        ! Rate equation:
        ! dt(X(H2)) = k_HH * X(H)**2 / D2G
        if (chem_solver_params%H2_form_use_moeq) then
          i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
          rtmp = chem_net%rates(i) * y(i1) * y(chem_net%reac(1, i))
          ydot(i1) = ydot(i1) - rtmp ! It's like H + gH -> gH2. So dt(H) -= rtmp, dt(gH) += rtmp
          ydot(chem_net%reac(1, i)) = ydot(chem_net%reac(1, i)) + rtmp
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
end subroutine chem_ode_f


subroutine chem_ode_jac(NEQ, t, y, j, ian, jan, pdj)
  use chemistry
  implicit none
  double precision t, rtmp
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, i, j, k, i1
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 64) ! A + B -> C
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
      case (1, 2, 3, 13, 61, 62, 0) ! A -> B
        if (j .EQ. chem_net%reac(1, i)) then
          rtmp = chem_net%rates(i)
        else
          rtmp = 0D0
        end if
      case (63) ! gA + gA -> gB
        if (chem_solver_params%H2_form_use_moeq) then
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
    !if ((j .EQ. chem_net%reac(1, i)) .OR. &
    !    (j .EQ. chem_net%reac(2, i))) then
    !  rtmp = 0D0
    !  if (chem_net%n_reac(i) .EQ. 1) then
    !    rtmp = chem_net%rates(i)
    !  else if (chem_net%n_reac(i) .EQ. 2) then
    !    if (.not. ((chem_net%itype(i) .eq. 0) .or. (chem_net%itype(i) .eq. 63))) then
    !      if (j .EQ. chem_net%reac(1, i)) then
    !        rtmp = chem_net%rates(i) * y(chem_net%reac(2, i))
    !      else if (j .EQ. chem_net%reac(2, i)) then
    !        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
    !      end if
    !    else
    !      rtmp = chem_net%rates(i)
    !    end if
    !  end if
    !end if
  end do
end subroutine chem_ode_jac
