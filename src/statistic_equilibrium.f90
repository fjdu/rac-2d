module statistic_equilibrium

use trivials
use phy_const
implicit none

integer, parameter, private :: const_len_energy_level = 12
integer, parameter, private :: const_len_molecule = 12


type :: type_energy_level
  character(len=const_len_energy_level) :: name_energy
  integer id
  double precision :: energy
  double precision :: weight
end type type_energy_level


type :: type_rad_transition
  double precision Eup, Elow, dE, freq, lambda
  double precision Aul, Bul, Blu, beta, J_ave, cooling_rate
  integer iup, ilow
end type type_rad_transition


type :: type_collisional_transition
  character(len=const_len_molecule) :: name_partner
  double precision dens_partner
  integer n_transition, n_T
  integer, dimension(:), allocatable :: iup, ilow
  double precision, dimension(:), allocatable :: T_coll
  double precision, dimension(:,:), allocatable :: Cul
end type type_collisional_transition


type :: type_rad_set
  integer n_transition
  type(type_rad_transition), dimension(:), allocatable :: list
end type type_rad_set


type :: type_colli_set
  integer n_partner
  type(type_collisional_transition), dimension(:), allocatable :: list
end type type_colli_set


type :: type_molecule_energy_set
  character(len=const_len_molecule) name_molecule
  double precision Tkin, density_mol, dv, length_scale, cooling_rate_total
  integer n_level
  type(type_energy_level), dimension(:), allocatable :: level_list
  double precision, dimension(:), allocatable :: f_occupation
  type(type_rad_set), allocatable :: rad_data
  type(type_colli_set), allocatable :: colli_data
end type type_molecule_energy_set


type :: type_statistic_equil_params
  integer nitem
  double precision :: RTOL = 1D-6, ATOL = 1D-20
  double precision :: t_max = 1D12, dt_first_step = 1D-4, ratio_tstep = 1.5D0
  real :: max_runtime_allowed = 20.0
  integer n_record
  integer :: &
        NERR, &
        NEQ, &
        ITOL = 1, &
        ITASK = 1, &
        ISTATE = 1, &
        IOPT = 1, &
        LIW = 128, &
        LRW = 4096, &
        MF = 21
  double precision, dimension(4096) :: RWORK
  integer, dimension(128) :: IWORK
end type type_statistic_equil_params

type(type_molecule_energy_set), pointer :: a_molecule_using

type(type_statistic_equil_params) statistic_equil_params

contains


subroutine reset_statistic_equil_params
  statistic_equil_params%NERR = 0
  statistic_equil_params%ITASK = 1
  statistic_equil_params%ISTATE = 1
  statistic_equil_params%IOPT = 1
  statistic_equil_params%IWORK(6) = 5000
end subroutine reset_statistic_equil_params


subroutine calc_cooling_rate
  integer i
  double precision f_up, f_low
  a_molecule_using%cooling_rate_total = 0D0
  do i=1, a_molecule_using%rad_data%n_transition
    associate( &
          n_mol => a_molecule_using%density_mol, &
          beta  => a_molecule_using%rad_data%list(i)%beta, &
          nu    => a_molecule_using%rad_data%list(i)%freq, &
          Aul   => a_molecule_using%rad_data%list(i)%Aul, &
          Bul   => a_molecule_using%rad_data%list(i)%Bul, &
          Blu   => a_molecule_using%rad_data%list(i)%Blu, &
          iup   => a_molecule_using%rad_data%list(i)%iup, &
          ilow  => a_molecule_using%rad_data%list(i)%ilow, &
          J_ave => a_molecule_using%rad_data%list(i)%J_ave)
      a_molecule_using%rad_data%list(i)%cooling_rate = &
        beta * phy_hPlanck_CGS * nu * n_mol * &
        ((Aul + Bul * J_ave) * a_molecule_using%f_occupation(iup) - &
          Blu * J_ave * a_molecule_using%f_occupation(ilow))
      a_molecule_using%cooling_rate_total = a_molecule_using%cooling_rate_total + &
        a_molecule_using%rad_data%list(i)%cooling_rate
    end associate
  end do
end subroutine calc_cooling_rate


subroutine load_moldata_LAMBDA(filename)
  character(len=*) filename
  character(len=512) strtmp
  integer, parameter :: nstr_split = 64
  character(len=32), dimension(nstr_split) :: str_split
  integer i, j, k, iup_, ilow_, fU, nout
  character(len=8), parameter :: strfmt_row = '(A512)'
  character(len=8), parameter :: strfmt_float = '(F16.3)'
  character(len=8), parameter :: strfmt_int = '(I6)'
  double precision, parameter :: freq_conv_factor = 1D9
  !
  integer n_T_, n_transition_
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a free file unit.  In load_moldata_LAMBDA.'
    stop
  end if
  call openFileSequentialRead(fU, filename, 99999)
  ! Get molecule name.
  read(fU,'(A1)') strtmp
  read(fU, strfmt_row) strtmp
  call split_str_by_space(strtmp, str_split, nstr_split, nout)
  a_molecule_using%name_molecule = str_split(1)
  ! Get energy level list
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU,'(I4)') a_molecule_using%n_level
  read(fU,'(A1)') strtmp
  allocate(a_molecule_using%level_list(a_molecule_using%n_level), &
           a_molecule_using%f_occupation(a_molecule_using%n_level))
  do i=1, a_molecule_using%n_level
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_float) a_molecule_using%level_list(i)%energy
    read(str_split(3), strfmt_float) a_molecule_using%level_list(i)%weight
  end do
  ! Convert the unit into Kelvin from cm-1
  a_molecule_using%level_list%energy = a_molecule_using%level_list%energy * phy_cm_1_2K
  !
  ! Get radiative transitions
  allocate(a_molecule_using%rad_data)
  read(fU,'(A1)') strtmp
  read(fU,'(I4)') a_molecule_using%rad_data%n_transition
  read(fU,'(A1)') strtmp
  allocate(a_molecule_using%rad_data%list(a_molecule_using%rad_data%n_transition))
  do i=1, a_molecule_using%rad_data%n_transition
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_int) a_molecule_using%rad_data%list(i)%iup
    read(str_split(3), strfmt_int) a_molecule_using%rad_data%list(i)%ilow
    read(str_split(4), strfmt_float) a_molecule_using%rad_data%list(i)%Aul
    read(str_split(5), strfmt_float) a_molecule_using%rad_data%list(i)%freq
    read(str_split(6), strfmt_float) a_molecule_using%rad_data%list(i)%Eup
  end do
  a_molecule_using%rad_data%list%freq = a_molecule_using%rad_data%list%freq * freq_conv_factor
  a_molecule_using%rad_data%list%Bul = a_molecule_using%rad_data%list%Aul / &
    ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * (a_molecule_using%rad_data%list%freq)**3)
  do i=1, a_molecule_using%rad_data%n_transition
    j = a_molecule_using%rad_data%list(i)%iup
    k = a_molecule_using%rad_data%list(i)%ilow
    a_molecule_using%rad_data%list(i)%Blu = a_molecule_using%rad_data%list(i)%Bul * &
        a_molecule_using%level_list(j)%weight / a_molecule_using%level_list(k)%weight
  end do
  !
  ! Collisional transitions
  allocate(a_molecule_using%colli_data)
  ! Get the number of collisional partners
  read(fU,'(A1)') strtmp
  read(fU,'(I4)') a_molecule_using%colli_data%n_partner
  allocate(a_molecule_using%colli_data%list(a_molecule_using%colli_data%n_partner))
  do i=1, a_molecule_using%colli_data%n_partner
    ! Get the name of partner
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    a_molecule_using%colli_data%list(i)%name_partner = str_split(4)
    ! Get the number of transitions and temperatures
    read(fU,'(A1)') strtmp
    read(fU,'(I4)') a_molecule_using%colli_data%list(i)%n_transition
    read(fU,'(A1)') strtmp
    read(fU,'(I4)') a_molecule_using%colli_data%list(i)%n_T
    !
    ! Name too long...
    n_transition_ = a_molecule_using%colli_data%list(i)%n_transition
    n_T_ = a_molecule_using%colli_data%list(i)%n_T
    if ((n_T_+3) .gt. nstr_split) then
      write(*,*) 'The number of different temperatures is too large!'
      write(*,*) 'nstr_split = ', nstr_split
      write(*,*) 'Change nstr_split of the source code to a higher value.'
      stop
    end if
    !
    allocate(a_molecule_using%colli_data%list(i)%iup(n_transition_), &
             a_molecule_using%colli_data%list(i)%ilow(n_transition_), &
             a_molecule_using%colli_data%list(i)%T_coll(n_T_), &
             a_molecule_using%colli_data%list(i)%Cul(n_T_, n_transition_))
    !
    ! Get the list of temperatures
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    do j=1, n_T_
      read(str_split(j), strfmt_float) a_molecule_using%colli_data%list(i)%T_coll(j)
    end do
    ! Get the collision coefficients
    read(fU,'(A1)') strtmp
    do j=1, n_transition_
      read(fU, strfmt_row) strtmp
      call split_str_by_space(strtmp, str_split, nstr_split, nout)
      read(str_split(2), strfmt_int) a_molecule_using%colli_data%list(i)%iup(j)
      read(str_split(3), strfmt_int) a_molecule_using%colli_data%list(i)%ilow(j)
      do k=1, n_T_
        read(str_split(3+k), strfmt_float) a_molecule_using%colli_data%list(i)%Cul(k, j)
        iup_ = a_molecule_using%colli_data%list(i)%iup(j)
        ilow_ = a_molecule_using%colli_data%list(i)%ilow(j)
      end do
    end do
  end do
  !
  close(fU)
  ! Test the results
  ! write(*,*) a_molecule_using%name_molecule
  ! do i=1, a_molecule_using%n_level
  !   write(*,*) i, a_molecule_using%level_list(i)%energy, a_molecule_using%level_list(i)%weight
  ! end do
  ! do i=1, a_molecule_using%rad_data%n_transition
  !   write(*,*) i, a_molecule_using%rad_data%list(i)%iup, &
  !     a_molecule_using%rad_data%list(i)%ilow, &
  !     a_molecule_using%rad_data%list(i)%Aul, &
  !     a_molecule_using%rad_data%list(i)%Bul, &
  !     a_molecule_using%rad_data%list(i)%Blu, &
  !     a_molecule_using%rad_data%list(i)%freq, &
  !     a_molecule_using%rad_data%list(i)%Eup
  ! end do
  ! do i=1, a_molecule_using%colli_data%n_partner
  !   write(*,*) i, a_molecule_using%colli_data%list(i)%name_partner
  !   write(*,*) a_molecule_using%colli_data%list(i)%T_coll
  !   do j=1, a_molecule_using%colli_data%list(i)%n_transition
  !     write(*,*) i, j, a_molecule_using%colli_data%list(i)%iup(j), &
  !       a_molecule_using%colli_data%list(i)%ilow(j), &
  !       a_molecule_using%colli_data%list(i)%Cul(:, j)
  !   end do
  ! end do
end subroutine load_moldata_LAMBDA



subroutine statistic_equil_solve
  use my_timer
  external stat_equili_ode_f, stat_equili_ode_jac
  integer i, itmp
  double precision t, tout, t_step, t_scale_min
  type(atimer) timer
  real time_thisstep, runtime_thisstep, time_laststep, runtime_laststep
  !
  call reset_statistic_equil_params
  !
  t = 0D0
  tout = statistic_equil_params%dt_first_step
  t_step = statistic_equil_params%dt_first_step
  !
  call timer%init('Stati_equil')
  time_laststep = timer%elapsed_time()
  runtime_laststep = huge(0.0)
  !
  statistic_equil_params%n_record = ceiling( &
    log(statistic_equil_params%t_max / statistic_equil_params%dt_first_step * &
        (statistic_equil_params%ratio_tstep - 1D0) + 1D0) &
    / log(statistic_equil_params%ratio_tstep))
  !
  do i=2, statistic_equil_params%n_record
    !write (*, '(A, 25X, "Solving... ", I5, " (", F5.1, "%)", "  t = ", ES9.2, "  tStep = ", ES9.2)') &
    !  CHAR(27)//'[A', i, real(i*100)/real(statistic_equil_params%n_record), t, t_step
    !
    call DLSODE( &
         stat_equili_ode_f, &
         !
         a_molecule_using%n_level, &
         a_molecule_using%f_occupation, &
         !
         t, &
         tout, &
         !
         statistic_equil_params%ITOL, &
         statistic_equil_params%RTOL, &
         statistic_equil_params%ATOL, &
         statistic_equil_params%ITASK, &
         statistic_equil_params%ISTATE, &
         statistic_equil_params%IOPT, &
         statistic_equil_params%RWORK, &
         statistic_equil_params%LRW, &
         statistic_equil_params%IWORK, &
         statistic_equil_params%LIW, &
         !
         stat_equili_ode_jac, &
         !
         statistic_equil_params%MF)
    !
    time_thisstep = timer%elapsed_time()
    runtime_thisstep = time_thisstep - time_laststep
    if ((runtime_thisstep .gt. max(5.0*runtime_laststep, 0.1*statistic_equil_params%max_runtime_allowed)) &
        .or. &
        (time_thisstep .gt. statistic_equil_params%max_runtime_allowed)) then
      write(*, '(A, ES9.2/)') 'Premature finish: t = ', t
      exit
    end if
    time_laststep = time_thisstep
    runtime_laststep = runtime_thisstep
    !
    if (statistic_equil_params%ISTATE .LT. 0) then
      statistic_equil_params%NERR = statistic_equil_params%NERR + 1
      write(*, '(A, I3/)') 'Error: ', statistic_equil_params%ISTATE
      statistic_equil_params%ISTATE = 3
    end if
    t_step = t_step * statistic_equil_params%ratio_tstep
    tout = t + t_step
  end do
  !
  do i=1, a_molecule_using%n_level
    if (a_molecule_using%f_occupation(i) .lt. 0D0) then
      a_molecule_using%f_occupation(i) = 0D0
    end if
  end do
  a_molecule_using%f_occupation = a_molecule_using%f_occupation / sum(a_molecule_using%f_occupation)
end subroutine statistic_equil_solve


end module statistic_equilibrium




subroutine stat_equili_ode_f(NEQ, t, y, ydot)
  ! sum(y) = 1
  use statistic_equilibrium
  use phy_const
  implicit none
  integer NEQ
  double precision t, y(NEQ), ydot(NEQ)
  integer i, j, itmp, iup, ilow, iL, iR
  double precision nu, J_ave, rtmp, Tkin, Cul, Clu, TL, TR, deltaE, del_nu, alpha, tau, beta
  double precision tmp1
  double precision, parameter :: const_small_num = 1D-8
  ydot = 0D0
  Tkin = a_molecule_using%Tkin
  do i=1, a_molecule_using%rad_data%n_transition
    iup = a_molecule_using%rad_data%list(i)%iup
    ilow = a_molecule_using%rad_data%list(i)%ilow
    nu = a_molecule_using%rad_data%list(i)%freq
    del_nu = nu * a_molecule_using%dv / phy_SpeedOfLight_CGS
    alpha = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_molecule_using%density_mol * &
            (y(ilow) * a_molecule_using%rad_data%list(i)%Blu - &
             y(iup)  * a_molecule_using%rad_data%list(i)%Bul) / del_nu
    tau = alpha * a_molecule_using%length_scale
    if (abs(tau) .le. const_small_num) then
      beta = 1D0
    else
      beta = (1D0 - exp(-3D0*tau)) / (3D0 * tau)
    end if
    tmp1 = phy_hPlanck_CGS*nu / (phy_kBoltzmann_CGS*Tkin)
    if (tmp1 .le. const_small_num) then
      J_ave = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / &
              tmp1 * (1D0 - beta)
    else
      J_ave = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / &
              (exp(tmp1) - 1D0) * (1D0 - beta)
    end if
    !
    a_molecule_using%rad_data%list(i)%beta = beta
    a_molecule_using%rad_data%list(i)%J_ave = J_ave
    !
    rtmp = a_molecule_using%rad_data%list(i)%Aul * y(iup) + &
           a_molecule_using%rad_data%list(i)%Bul * J_ave * y(iup) - &
           a_molecule_using%rad_data%list(i)%Blu * J_ave * y(ilow)
    ydot(iup) = ydot(iup)   - rtmp
    ydot(ilow) = ydot(ilow) + rtmp
  end do
  do i=1, a_molecule_using%colli_data%n_partner
    ! Find the T interval
    itmp = a_molecule_using%colli_data%list(i)%n_T
    if (Tkin .le. a_molecule_using%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_molecule_using%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, a_molecule_using%colli_data%list(i)%n_T
        if ((Tkin .ge. a_molecule_using%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. a_molecule_using%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    do j=1, a_molecule_using%colli_data%list(i)%n_transition
      iup = a_molecule_using%colli_data%list(i)%iup(j)
      ilow = a_molecule_using%colli_data%list(i)%ilow(j)
      deltaE = a_molecule_using%level_list(iup)%energy - a_molecule_using%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = a_molecule_using%colli_data%list(i)%Cul(iL, j)
      else
        TL = a_molecule_using%colli_data%list(i)%T_coll(iL)
        TR = a_molecule_using%colli_data%list(i)%T_coll(iR)
        Cul = (a_molecule_using%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
                a_molecule_using%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_molecule_using%level_list(iup)%weight / &
             a_molecule_using%level_list(ilow)%weight
      rtmp = (Cul * y(iup) - Clu * y(ilow)) * a_molecule_using%colli_data%list(i)%dens_partner
      ydot(iup) = ydot(iup)   - rtmp
      ydot(ilow) = ydot(ilow) + rtmp
    end do
  end do
end subroutine stat_equili_ode_f


subroutine stat_equili_ode_jac(NEQ, t, y, ML, MU, PD, NROWPD)
  use statistic_equilibrium
  use phy_const
  implicit none
  double precision t
  integer ML, MU, NROWPD
  double precision, dimension(NEQ) :: y
  double precision, dimension(NROWPD, *) :: PD
  integer NEQ
  integer i, j, itmp, iup, ilow, iL, iR
  double precision nu, J_ave, &
        Tkin, Cul, Clu, TL, TR, deltaE, del_nu, alpha, tau, beta
  double precision S, dbeta_dtau, dtau_dy_up, dtau_dy_low, &
    dJ_ave_dy_up, dJ_ave_dy_low, drtmp_dy_up, drtmp_dy_low
  double precision tmp1
  double precision, parameter :: const_small_num = 1D-8
  Tkin = a_molecule_using%Tkin
  do i=1, a_molecule_using%rad_data%n_transition
    iup = a_molecule_using%rad_data%list(i)%iup
    ilow = a_molecule_using%rad_data%list(i)%ilow
    nu = a_molecule_using%rad_data%list(i)%freq
    del_nu = nu * a_molecule_using%dv / phy_SpeedOfLight_CGS
    alpha = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_molecule_using%density_mol * &
            (y(ilow) * a_molecule_using%rad_data%list(i)%Blu - &
             y(iup)  * a_molecule_using%rad_data%list(i)%Bul) / del_nu
    tau = alpha * a_molecule_using%length_scale
    if (abs(tau) .le. const_small_num) then
      beta = 1D0
    else
      beta = (1D0 - exp(-3D0*tau)) / (3D0 * tau)
    end if
    tmp1 = phy_hPlanck_CGS*nu / (phy_kBoltzmann_CGS*Tkin)
    if (tmp1 .le. const_small_num) then
      S = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / tmp1
    else
      S = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / (exp(tmp1) - 1D0)
    end if
    J_ave = S * (1D0 - beta)
    !
    dbeta_dtau = exp(-3D0*tau) * (1D0/tau + 1D0/3D0/tau/tau) - 1D0/3D0/tau/tau
    dtau_dy_up = a_molecule_using%length_scale * &
                 phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_molecule_using%density_mol * &
                 (-a_molecule_using%rad_data%list(i)%Bul) / del_nu
    dtau_dy_low = a_molecule_using%length_scale * &
                  phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_molecule_using%density_mol * &
                  (a_molecule_using%rad_data%list(i)%Blu) / del_nu
    dJ_ave_dy_up  = -S * dbeta_dtau * dtau_dy_up
    dJ_ave_dy_low = -S * dbeta_dtau * dtau_dy_low
    !
    drtmp_dy_up = a_molecule_using%rad_data%list(i)%Aul + &
             a_molecule_using%rad_data%list(i)%Bul * J_ave + &
             (a_molecule_using%rad_data%list(i)%Bul * y(iup) - &
              a_molecule_using%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_up
    drtmp_dy_low = -a_molecule_using%rad_data%list(i)%Blu * J_ave + &
             (a_molecule_using%rad_data%list(i)%Bul * y(iup) - &
              a_molecule_using%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_low
    PD(iup,  iup)  = PD(iup,  iup)  - drtmp_dy_up
    PD(ilow, iup)  = PD(ilow, iup)  + drtmp_dy_up
    PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
    PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
  end do
  do i=1, a_molecule_using%colli_data%n_partner
    ! Find the T interval
    itmp = a_molecule_using%colli_data%list(i)%n_T
    if (Tkin .le. a_molecule_using%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_molecule_using%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, a_molecule_using%colli_data%list(i)%n_T
        if ((Tkin .ge. a_molecule_using%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. a_molecule_using%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    do j=1, a_molecule_using%colli_data%list(i)%n_transition
      iup = a_molecule_using%colli_data%list(i)%iup(j)
      ilow = a_molecule_using%colli_data%list(i)%ilow(j)
      deltaE = a_molecule_using%level_list(iup)%energy - a_molecule_using%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = a_molecule_using%colli_data%list(i)%Cul(iL, j)
      else
        TL = a_molecule_using%colli_data%list(i)%T_coll(iL)
        TR = a_molecule_using%colli_data%list(i)%T_coll(iR)
        Cul = (a_molecule_using%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
                a_molecule_using%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_molecule_using%level_list(iup)%weight / &
             a_molecule_using%level_list(ilow)%weight
      drtmp_dy_up  = Cul  * a_molecule_using%colli_data%list(i)%dens_partner
      drtmp_dy_low = -Clu * a_molecule_using%colli_data%list(i)%dens_partner
      PD(iup,  iup)  = PD(iup,  iup) - drtmp_dy_up
      PD(ilow, iup)  = PD(ilow, iup) + drtmp_dy_up
      PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
      PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
    end do
  end do
end subroutine stat_equili_ode_jac
