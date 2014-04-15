module statistic_equilibrium

use data_struct
use phy_const

implicit none

type(type_molecule_energy_set), pointer :: mol_sta_sol

type(type_statistic_equil_params) statistic_equil_params

type(type_continuum_lut), pointer :: cont_lut_ptr


contains


subroutine reset_statistic_equil_params
  statistic_equil_params%is_good = .true.
  statistic_equil_params%NERR = 0
  statistic_equil_params%ITASK = 1
  statistic_equil_params%ISTATE = 1
  statistic_equil_params%IOPT = 1
  !
  statistic_equil_params%RWORK = 0D0
  statistic_equil_params%IWORK = 0
  statistic_equil_params%IWORK(6) = 5000
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
         mol_sta_sol%n_level, &
         mol_sta_sol%f_occupation, &
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
  do i=1, mol_sta_sol%n_level
    if (mol_sta_sol%f_occupation(i) .lt. 0D0) then
      mol_sta_sol%f_occupation(i) = 0D0
    end if
  end do
  mol_sta_sol%f_occupation = mol_sta_sol%f_occupation / sum(mol_sta_sol%f_occupation)
end subroutine statistic_equil_solve


subroutine get_cont_alpha(lam, alp, J)
  double precision, intent(in) :: lam
  double precision, intent(out) :: alp, J
  integer i, imin, imax, imid, k
  integer, parameter :: ITH = 5
  if (cont_lut_ptr%n .le. 0) then
    alp = 0D0
    J = 0D0
    return
  end if
  if (lam .lt. cont_lut_ptr%lam(1)) then
    alp = cont_lut_ptr%alpha(1)
    J = cont_lut_ptr%J(1)
  else if (lam .gt. cont_lut_ptr%lam(cont_lut_ptr%n)) then
    alp = cont_lut_ptr%alpha(cont_lut_ptr%n)
    J = cont_lut_ptr%J(cont_lut_ptr%n)
  else
    imin = 1
    imax = cont_lut_ptr%n
    do i=1, cont_lut_ptr%n
      if (imin .ge. imax-ITH) then
        do k=imin, imax-1
          if ((cont_lut_ptr%lam(k) .le. lam) .and. &
              (cont_lut_ptr%lam(k+1) .gt. lam)) then
            alp = cont_lut_ptr%alpha(k)
            J = cont_lut_ptr%J(k)
            return
          end if
        end do
        exit
      else
        imid = (imin + imax) / 2
        if (cont_lut_ptr%lam(imid) .le. lam) then
          imin = imid
        else
          imax = imid
        end if
      end if
    end do
  end if
end subroutine get_cont_alpha


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
  double precision lambda, cont_alpha, cont_J
  !double precision tmp1
  double precision jnu, knu
  double precision, parameter :: const_small_num = 1D-6
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
    if (abs(tau) .le. const_small_num) then
      beta = 1D0
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
      J_ave = jnu / knu
    !else if ((knu .ge. 0D0) .and. (knu .le. 1D-50)) then
    !  J_ave = jnu / (knu + 1D-50)
    !else
    !  J_ave = jnu / (knu - 1D-50)
    else
      J_ave = jnu * mol_sta_sol%length_scale * t1
    end if
    !
    J_ave = J_ave * (1D0 - beta) + cont_J
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
      rtmp = (Cul * y(iup) - Clu * y(ilow)) * mol_sta_sol%colli_data%list(i)%dens_partner
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
  double precision S, dbeta_dtau, dtau_dy_up, dtau_dy_low, &
    dJ_ave_dy_up, dJ_ave_dy_low, drtmp_dy_up, drtmp_dy_low, &
    dS_dy_up, dS_dy_low
  !double precision tmp1
  double precision jnu, knu
  double precision t1
  double precision, parameter :: const_small_num = 1D-6
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
    J_ave = S * (1D0 - beta) + cont_J
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
    dJ_ave_dy_up  = -S * dbeta_dtau * dtau_dy_up + dS_dy_up * (1D0 - beta)
    dJ_ave_dy_low = -S * dbeta_dtau * dtau_dy_low + dS_dy_low * (1D0 - beta)
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
