module statistic_equilibrium
! Purpose: solve the the occupation of energy levels in statistical
!          equilibrium based on LVG.
!
! Fujun Du

use data_struct
use montecarlo
use phy_const

implicit none

type(type_molecule_energy_set), pointer :: mol_sta_sol

type(type_cell), pointer :: current_cell_ptr

type(type_statistic_equil_params), private :: sta_equil_params


contains



subroutine init_statistic_sol(neq, meth)
  ! Purpose: allocate memory to be used by the solver.
  !
  ! Note: if more than one molecule are to be solved, then NEQ may be different
  !       for each molecule, so the allocated momory needs to be enough for
  !       the largest set of possible energy levels.
  !
  integer, intent(in) :: neq, meth
  !
  sta_equil_params%LIW = 50 + neq
  sta_equil_params%LRW = (neq + 13) * neq + 61
  if (.not. allocated(sta_equil_params%IWORK)) then
    allocate(sta_equil_params%IWORK(sta_equil_params%LIW), &
             sta_equil_params%RWORK(sta_equil_params%LRW))
  end if
end subroutine init_statistic_sol



subroutine reset_statistic_equil_params
  sta_equil_params%is_good = .true.
  sta_equil_params%NERR = 0
  sta_equil_params%ITASK = 1
  sta_equil_params%ISTATE = 1
  sta_equil_params%IOPT = 1
  !
  sta_equil_params%RWORK = 0D0
  sta_equil_params%IWORK = 0
  sta_equil_params%IWORK(6) = 5000
end subroutine reset_statistic_equil_params


subroutine calc_cooling_rate
  integer i
  mol_sta_sol%cooling_rate_total = 0D0
  do i=1, mol_sta_sol%rad_data%n_transition
    associate( &
          n_mol => mol_sta_sol%density_mol, &
          beta  => mol_sta_sol%rad_data%list(i)%beta, &
          nu    => mol_sta_sol%rad_data%list(i)%freq, &
          Aul   => mol_sta_sol%rad_data%list(i)%Aul, &
          Bul   => mol_sta_sol%rad_data%list(i)%Bul, &
          Blu   => mol_sta_sol%rad_data%list(i)%Blu, &
          iup   => mol_sta_sol%rad_data%list(i)%iup, &
          ilow  => mol_sta_sol%rad_data%list(i)%ilow, &
          J_ave => mol_sta_sol%rad_data%list(i)%J_ave)
      mol_sta_sol%rad_data%list(i)%cooling_rate = &
        beta * phy_hPlanck_CGS * nu * n_mol * &
        ((Aul + Bul * J_ave) * mol_sta_sol%f_occupation(iup) - &
          Blu * J_ave * mol_sta_sol%f_occupation(ilow))
      mol_sta_sol%cooling_rate_total = mol_sta_sol%cooling_rate_total + &
        mol_sta_sol%rad_data%list(i)%cooling_rate
    end associate
  end do
end subroutine calc_cooling_rate


subroutine statistic_equil_solve
  use my_timer
  external stat_equili_ode_f, stat_equili_ode_jac
  integer i
  double precision t, tout, t_step
  type(atimer) timer
  real time_thisstep, runtime_thisstep, time_laststep, runtime_laststep
  !
  call reset_statistic_equil_params
  !
  t = 0D0
  tout = sta_equil_params%dt_first_step
  t_step = sta_equil_params%dt_first_step
  !
  call timer%init('Stati_equil')
  time_laststep = timer%elapsed_time()
  runtime_laststep = 1D30  ! huge(0.0)
  !
  sta_equil_params%n_record = ceiling( &
    log(sta_equil_params%t_max / sta_equil_params%dt_first_step * &
        (sta_equil_params%ratio_tstep - 1D0) + 1D0) &
    / log(sta_equil_params%ratio_tstep))
  !
  do i=2, sta_equil_params%n_record
    !write (*, '(A, 25X, "StaEq-Solving... ", I5, " (", F5.1, "%)", "  t = ", ES9.2, "  tStep = ", ES9.2)') &
    !  CHAR(27)//'[A', i, real(i*100)/real(sta_equil_params%n_record), t, t_step
    !
    call DLSODE( &
         stat_equili_ode_f, &
         !
         mol_sta_sol%n_level, &
         mol_sta_sol%f_occupation, &
         !
         t, &
         tout, &
         !
         sta_equil_params%ITOL, &
         sta_equil_params%RTOL, &
         sta_equil_params%ATOL, &
         sta_equil_params%ITASK, &
         sta_equil_params%ISTATE, &
         sta_equil_params%IOPT, &
         sta_equil_params%RWORK, &
         sta_equil_params%LRW, &
         sta_equil_params%IWORK, &
         sta_equil_params%LIW, &
         !
         stat_equili_ode_jac, &
         !
         sta_equil_params%MF)
    !
    time_thisstep = timer%elapsed_time()
    runtime_thisstep = time_thisstep - time_laststep
    if ((runtime_thisstep .gt. max(5.0*runtime_laststep, 0.1*sta_equil_params%max_runtime_allowed)) &
        .or. &
        (time_thisstep .gt. sta_equil_params%max_runtime_allowed)) then
      write(*, '(A, ES9.2/)') 'Running too slow... Premature finish: t = ', t
      exit
    end if
    time_laststep = time_thisstep
    runtime_laststep = runtime_thisstep
    !
    if (sta_equil_params%ISTATE .LT. 0) then
      sta_equil_params%NERR = sta_equil_params%NERR + 1
      write(*, '(A, I3/)') 'Error: ', sta_equil_params%ISTATE
      sta_equil_params%ISTATE = 3
    end if
    t_step = t_step * sta_equil_params%ratio_tstep
    tout = t + t_step
  end do
  !
  do i=1, mol_sta_sol%n_level
    if (mol_sta_sol%f_occupation(i) .lt. 0D0) then
      mol_sta_sol%f_occupation(i) = 0D0
    end if
  end do
  mol_sta_sol%f_occupation = mol_sta_sol%f_occupation / sum(mol_sta_sol%f_occupation)
end subroutine statistic_equil_solve


subroutine statistic_equil_solve_Newton
  external stat_equili_fcn, stat_equili_jac
  double precision, dimension(mol_sta_sol%n_level) :: XSCAL
  integer IERR
  integer, dimension(50) :: IOPT
  !
  XSCAL = mol_sta_sol%f_occupation*1D-4 + 1D-20
  sta_equil_params%RTOL = 1D-4
  IOPT = 0
  !
  sta_equil_params%IWORK = 0
  sta_equil_params%RWORK = 0D0
  !
  IOPT(3) = 1 ! JACGEN
  IOPT(11) = 0 ! MPRERR
  IOPT(31) = 4 ! NONLIN
  !
  sta_equil_params%IWORK(1) = 50000 ! NITER
  !
  call NLEQ1( &
             mol_sta_sol%n_level, &
             stat_equili_fcn, &
             stat_equili_jac, &
             mol_sta_sol%f_occupation, &
             XSCAL, &
             sta_equil_params%RTOL, &
             IOPT, &
             IERR, &
             sta_equil_params%LIW, &
             sta_equil_params%IWORK, &
             sta_equil_params%LRW, &
             sta_equil_params%RWORK)
  if (IERR .eq. 10) then
    write(*, '(/A)') 'LIW too small!'
    call error_stop()
  end if
  mol_sta_sol%f_occupation = mol_sta_sol%f_occupation / sum(mol_sta_sol%f_occupation)
end subroutine statistic_equil_solve_Newton


subroutine get_cont_alpha(lam, alp, J)
  double precision, intent(in) :: lam
  double precision, intent(out) :: alp, J
  integer i, imin, imax, imid, k
  integer, parameter :: ITH = 5
  if (current_cell_ptr%cont_lut%n .le. 0) then
    alp = 0D0
    J = 0D0
    return
  end if
  if (lam .lt. dust_0%lam(1)) then
    alp = current_cell_ptr%optical%ext_tot(1)
    J = current_cell_ptr%cont_lut%J(1)
  else if (lam .gt. dust_0%lam(current_cell_ptr%cont_lut%n)) then
    alp = current_cell_ptr%optical%ext_tot(current_cell_ptr%cont_lut%n)
    J = current_cell_ptr%cont_lut%J(current_cell_ptr%cont_lut%n)
  else
    imin = 1
    imax = current_cell_ptr%cont_lut%n
    do i=1, current_cell_ptr%cont_lut%n
      if (imin .ge. imax-ITH) then
        do k=imin, imax-1
          if ((dust_0%lam(k) .le. lam) .and. &
              (dust_0%lam(k+1) .gt. lam)) then
            alp = current_cell_ptr%optical%ext_tot(k)
            J = current_cell_ptr%cont_lut%J(k)
            return
          end if
        end do
        exit
      else
        imid = (imin + imax) / 2
        if (dust_0%lam(imid) .le. lam) then
          imin = imid
        else
          imax = imid
        end if
      end if
    end do
  end if
end subroutine get_cont_alpha


end module statistic_equilibrium


subroutine stat_equili_fcn(N, X, F, IFAIL)
  integer, intent(in) :: N
  double precision, intent(in), dimension(N) :: X
  double precision, intent(out), dimension(N) :: F
  integer, intent(inout) :: IFAIL
  integer i
  call stat_equili_ode_f(N, 0D0, X, F)
  F(N) = 0D0
  do i=1, N
    F(N) = F(N) + X(i)
  end do
  F(N) = F(N) - 1D0
  IFAIL = 0
end subroutine stat_equili_fcn


subroutine stat_equili_jac(N, LDJAC, X, DFDX, IFAIL)
  integer, intent(in) :: N, LDJAC
  double precision, intent(in), dimension(N) :: X
  double precision, intent(out), dimension(LDJAC, N) :: DFDX
  integer, intent(inout) :: IFAIL
  call stat_equili_ode_jac(N, 0D0, X, 0, 0, DFDX, LDJAC)
  DFDX(LDJAC, 1:N) = 1D0
  IFAIL = 0
end subroutine stat_equili_jac


subroutine stat_equili_ode_f(NEQ, t, y, ydot)
  ! sum(y) = 1
  use statistic_equilibrium
  use phy_const
  implicit none
  integer NEQ
  double precision t, y(NEQ), ydot(NEQ)
  integer i, j, itmp, iup, ilow, iL, iR
  double precision nu, J_ave, rtmp, Tkin, Cul, Clu, TL, TR, deltaE, del_nu, alpha, tau, beta
  double precision lambda, cont_alpha, cont_J
  !double precision tmp1
  double precision jnu, knu, S
  double precision, parameter :: const_small_num = 1D-6
  double precision, parameter :: const_big_num = 100D0
  double precision t1
  ydot = 0D0
  Tkin = mol_sta_sol%Tkin
  do i=1, mol_sta_sol%rad_data%n_transition
    iup = mol_sta_sol%rad_data%list(i)%iup
    ilow = mol_sta_sol%rad_data%list(i)%ilow
    nu = mol_sta_sol%rad_data%list(i)%freq
    lambda = mol_sta_sol%rad_data%list(i)%lambda
    del_nu = nu * mol_sta_sol%dv / phy_SpeedOfLight_CGS
    call get_cont_alpha(lambda, cont_alpha, cont_J)
    !
    t1 = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol / del_nu
    jnu = y(iup) *  mol_sta_sol%rad_data%list(i)%Aul
    knu = y(ilow) * mol_sta_sol%rad_data%list(i)%Blu - &
          y(iup)  * mol_sta_sol%rad_data%list(i)%Bul
    alpha = t1 * knu + cont_alpha
    !alpha = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol * &
    !        (y(ilow) * mol_sta_sol%rad_data%list(i)%Blu - &
    !         y(iup)  * mol_sta_sol%rad_data%list(i)%Bul) / del_nu + &
    !         cont_alpha
    tau = alpha * mol_sta_sol%length_scale
    if (isnan(tau)) then
      write(*, *) 'tau is NaN'
      write(*, *) 'In stat_equili_ode_f'
      call error_stop()
    end if
    if (abs(tau) .le. const_small_num) then
      beta = 1D0
    else if (tau .ge. const_big_num) then
      beta = 1D0  / (3D0 * tau)
    else if (tau .lt. 0D0) then
#ifdef DIAGNOSIS_TRACK_FUNC_CALL_DETAILED
      write(*, *) 'tau is very negative: ', tau
      write(*, *) 'Tkin, density_mol, ilow, iup, y(ilow), y(iup), Blu, Bul, f, cont_alpha, length_scale:'
      write(*, *) Tkin, mol_sta_sol%density_mol, ilow, iup, y(ilow), y(iup), &
        mol_sta_sol%rad_data%list(i)%Blu, mol_sta_sol%rad_data%list(i)%Bul, &
        mol_sta_sol%rad_data%list(i)%freq, &
        cont_alpha, mol_sta_sol%length_scale
      write(*, *) 'In stat_equili_ode_f'
#endif
      beta = 1D0 - 1.5D0 * tau  ! Special treatment for negative tau
    else
      beta = (1D0 - exp(-3D0*tau)) / (3D0 * tau)
    end if
    !tmp1 = phy_hPlanck_CGS*nu / (phy_kBoltzmann_CGS*Tkin)
    !if (tmp1 .le. const_small_num) then
    !  J_ave = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / &
    !          tmp1 * (1D0 - beta)
    !else
    !  J_ave = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / &
    !          (exp(tmp1) - 1D0) * (1D0 - beta)
    !end if
    !
    if ((knu .gt. 1D-30) .or. (knu .lt. -1D-30)) then
      S = jnu / knu
    !else if ((knu .ge. 0D0) .and. (knu .le. 1D-50)) then
    !  J_ave = jnu / (knu + 1D-50)
    !else
    !  J_ave = jnu / (knu - 1D-50)
    else
      S = jnu * mol_sta_sol%length_scale * t1
    end if
    !
    J_ave = S * (1D0 - beta) + cont_J * beta
    !
    mol_sta_sol%rad_data%list(i)%beta = beta
    mol_sta_sol%rad_data%list(i)%J_ave = J_ave
    !
    rtmp = mol_sta_sol%rad_data%list(i)%Aul * y(iup) + &
           mol_sta_sol%rad_data%list(i)%Bul * J_ave * y(iup) - &
           mol_sta_sol%rad_data%list(i)%Blu * J_ave * y(ilow)
    ydot(iup) = ydot(iup)   - rtmp
    ydot(ilow) = ydot(ilow) + rtmp
  end do
  do i=1, mol_sta_sol%colli_data%n_partner
    ! Find the T interval
    iL = -1; iR = -1
    itmp = mol_sta_sol%colli_data%list(i)%n_T
    if (Tkin .le. mol_sta_sol%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. mol_sta_sol%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, mol_sta_sol%colli_data%list(i)%n_T
        if ((Tkin .ge. mol_sta_sol%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. mol_sta_sol%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    if ((iL .le. 0) .or. (iR .le. 0)) then
      write(*, '(A, 2I6)') 'Invalid iL or iR index: ', iL, iR
      call error_stop()
    end if
    do j=1, mol_sta_sol%colli_data%list(i)%n_transition
      iup = mol_sta_sol%colli_data%list(i)%iup(j)
      ilow = mol_sta_sol%colli_data%list(i)%ilow(j)
      deltaE = mol_sta_sol%level_list(iup)%energy - mol_sta_sol%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = mol_sta_sol%colli_data%list(i)%Cul(iL, j)
      else
        TL = mol_sta_sol%colli_data%list(i)%T_coll(iL)
        TR = mol_sta_sol%colli_data%list(i)%T_coll(iR)
        Cul = (mol_sta_sol%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
                mol_sta_sol%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             mol_sta_sol%level_list(iup)%weight / &
             mol_sta_sol%level_list(ilow)%weight
      rtmp = (Cul * y(iup) - Clu * y(ilow)) &
             * mol_sta_sol%colli_data%list(i)%dens_partner
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
  double precision lambda, cont_alpha, cont_J
  double precision dbeta_dtau, dtau_dy_up, dtau_dy_low, &
    dJ_ave_dy_up, dJ_ave_dy_low, drtmp_dy_up, drtmp_dy_low, &
    dS_dy_up, dS_dy_low
  !double precision tmp1
  double precision jnu, knu, S
  double precision t1
  double precision, parameter :: const_small_num = 1D-6
  double precision, parameter :: const_big_num = 100D0
  !
  PD(1:NROWPD, 1:NEQ) = 0D0
  !
  Tkin = mol_sta_sol%Tkin
  do i=1, mol_sta_sol%rad_data%n_transition
    iup = mol_sta_sol%rad_data%list(i)%iup
    ilow = mol_sta_sol%rad_data%list(i)%ilow
    nu = mol_sta_sol%rad_data%list(i)%freq
    lambda = mol_sta_sol%rad_data%list(i)%lambda
    del_nu = nu * mol_sta_sol%dv / phy_SpeedOfLight_CGS
    call get_cont_alpha(lambda, cont_alpha, cont_J)
    !
    t1 = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol / del_nu
    jnu = y(iup) *  mol_sta_sol%rad_data%list(i)%Aul
    knu = y(ilow) * mol_sta_sol%rad_data%list(i)%Blu - &
          y(iup)  * mol_sta_sol%rad_data%list(i)%Bul
    alpha = t1 * knu + cont_alpha
    !alpha = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol * &
    !        (y(ilow) * mol_sta_sol%rad_data%list(i)%Blu - &
    !         y(iup)  * mol_sta_sol%rad_data%list(i)%Bul) / del_nu + &
    !         cont_alpha
    tau = alpha * mol_sta_sol%length_scale
    if (abs(tau) .le. const_small_num) then
      beta = 1D0
      dbeta_dtau = -1.5D0
    else if (tau .ge. const_big_num) then
      beta = 1D0  / (3D0 * tau)
      dbeta_dtau =  -1D0/3D0/tau/tau
    else if (tau .lt. 0D0) then
#ifdef DIAGNOSIS_TRACK_FUNC_CALL_DETAILED
      write(*, *) 'tau is very negative: ', tau
      write(*, *) 'Tkin, density_mol, ilow, iup, y(ilow), y(iup), Blu, Bul, f, cont_alpha, length_scale:'
      write(*, *) Tkin, mol_sta_sol%density_mol, ilow, iup, y(ilow), y(iup), &
        mol_sta_sol%rad_data%list(i)%Blu, mol_sta_sol%rad_data%list(i)%Bul, &
        mol_sta_sol%rad_data%list(i)%freq, &
        cont_alpha, mol_sta_sol%length_scale
      write(*, *) 'In stat_equili_ode_f'
#endif
      beta = 1D0 - 1.5D0 * tau
      dbeta_dtau = -1.5D0
    else
      beta = (1D0 - exp(-3D0*tau)) / (3D0 * tau)
      dbeta_dtau = exp(-3D0*tau) * (1D0/tau + 1D0/3D0/tau/tau) - 1D0/3D0/tau/tau
    end if
    !tmp1 = phy_hPlanck_CGS*nu / (phy_kBoltzmann_CGS*Tkin)
    !if (tmp1 .le. const_small_num) then
    !  S = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / tmp1
    !else
    !  S = 2D0 * phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / (exp(tmp1) - 1D0)
    !end if
    !
    !if ((knu .ge. 0D0) .and. (knu .le. 1D-50)) then
    !  knu = knu + 1D-50
    !else if ((knu .lt. 0D0) .and. (knu .ge. -1D-50)) then
    !  knu = knu - 1D-50
    !end if
    if ((knu .gt. 1D-30) .or. (knu .lt. -1D-30)) then
      S = jnu / knu
      dS_dy_up = (mol_sta_sol%rad_data%list(i)%Aul + &
                  S * mol_sta_sol%rad_data%list(i)%Bul) / knu
      dS_dy_low = -S * mol_sta_sol%rad_data%list(i)%Blu / knu
    else
      S = jnu * mol_sta_sol%length_scale * t1
      dS_dy_up = mol_sta_sol%rad_data%list(i)%Aul * mol_sta_sol%length_scale * t1
      dS_dy_low = 0D0
    end if
    !
    J_ave = S * (1D0 - beta) + cont_J * beta
    !
    dtau_dy_up = mol_sta_sol%length_scale * &
                 phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol * &
                 (-mol_sta_sol%rad_data%list(i)%Bul) / del_nu
    dtau_dy_low = mol_sta_sol%length_scale * &
                  phy_hPlanck_CGS * nu / (4D0*phy_Pi) * mol_sta_sol%density_mol * &
                  (mol_sta_sol%rad_data%list(i)%Blu) / del_nu
    !if (knu .eq. 0D0) then
    !  dS_dy_up = 0D0
    !  dS_dy_low = 0D0
    !else
    !  dS_dy_up = (mol_sta_sol%rad_data%list(i)%Aul + &
    !              S * mol_sta_sol%rad_data%list(i)%Bul) / knu
    !  dS_dy_low = -S * mol_sta_sol%rad_data%list(i)%Blu / knu
    !end if
    dJ_ave_dy_up  = (cont_J-S) * dbeta_dtau * dtau_dy_up  + dS_dy_up * (1D0 - beta)
    dJ_ave_dy_low = (cont_J-S) * dbeta_dtau * dtau_dy_low + dS_dy_low * (1D0 - beta)
    !
    drtmp_dy_up = mol_sta_sol%rad_data%list(i)%Aul + &
             mol_sta_sol%rad_data%list(i)%Bul * J_ave + &
             (mol_sta_sol%rad_data%list(i)%Bul * y(iup) - &
              mol_sta_sol%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_up
    drtmp_dy_low = -mol_sta_sol%rad_data%list(i)%Blu * J_ave + &
             (mol_sta_sol%rad_data%list(i)%Bul * y(iup) - &
              mol_sta_sol%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_low
    PD(iup,  iup)  = PD(iup,  iup)  - drtmp_dy_up
    PD(ilow, iup)  = PD(ilow, iup)  + drtmp_dy_up
    PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
    PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
  end do
  do i=1, mol_sta_sol%colli_data%n_partner
    ! Find the T interval
    iL = -1; iR = -1
    itmp = mol_sta_sol%colli_data%list(i)%n_T
    if (Tkin .le. mol_sta_sol%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. mol_sta_sol%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, mol_sta_sol%colli_data%list(i)%n_T
        if ((Tkin .ge. mol_sta_sol%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. mol_sta_sol%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    if ((iL .le. 0) .or. (iR .le. 0)) then
      write(*, '(A, 2I6)') 'Invalid iL or iR index: ', iL, iR
      call error_stop()
    end if
    do j=1, mol_sta_sol%colli_data%list(i)%n_transition
      iup = mol_sta_sol%colli_data%list(i)%iup(j)
      ilow = mol_sta_sol%colli_data%list(i)%ilow(j)
      deltaE = mol_sta_sol%level_list(iup)%energy - mol_sta_sol%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = mol_sta_sol%colli_data%list(i)%Cul(iL, j)
      else
        TL = mol_sta_sol%colli_data%list(i)%T_coll(iL)
        TR = mol_sta_sol%colli_data%list(i)%T_coll(iR)
        Cul = (mol_sta_sol%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
                mol_sta_sol%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             mol_sta_sol%level_list(iup)%weight / &
             mol_sta_sol%level_list(ilow)%weight
      drtmp_dy_up  = Cul  * mol_sta_sol%colli_data%list(i)%dens_partner
      drtmp_dy_low = -Clu * mol_sta_sol%colli_data%list(i)%dens_partner
      PD(iup,  iup)  = PD(iup,  iup) - drtmp_dy_up
      PD(ilow, iup)  = PD(ilow, iup) + drtmp_dy_up
      PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
      PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
    end do
  end do
end subroutine stat_equili_ode_jac
