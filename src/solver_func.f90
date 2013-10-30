subroutine chem_ode_f(NEQ, t, y, ydot)
  use chemistry
  use heating_cooling
  implicit none
  integer NEQ, i, j, i1
  double precision t, y(NEQ), ydot(NEQ), rtmp
  ydot = 0D0
  !
  if (chem_solver_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (y(chem_species%nSpecies+1) .ne. chem_params%Tgas) then
      chem_params%Tgas = y(chem_species%nSpecies+1)
      call chem_cal_rates
    end if
  end if
  !
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
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          if (chem_solver_params%H2_form_use_moeq) then
            i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
            rtmp = chem_net%rates(i) * y(i1) * y(chem_net%reac(1, i))
            ydot(i1) = ydot(i1) - rtmp ! It's like H + gH -> gH2. So dt(H) -= rtmp, dt(gH) += rtmp
            ydot(chem_net%reac(1, i)) = ydot(chem_net%reac(1, i)) + rtmp
          else
            rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
          end if
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
  !
  if (chem_solver_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    call realtime_heating_cooling_rate(ydot(chem_species%nSpecies+1), NEQ, y)
  else
    ydot(chem_species%nSpecies+1) = 0D0
  end if
  !
end subroutine chem_ode_f




subroutine realtime_heating_cooling_rate(r, NEQ, y)
  use chemistry
  use heating_cooling
  double precision, intent(out) :: r
  integer, intent(in) :: NEQ
  double precision, dimension(NEQ), intent(in) :: y
  heating_cooling_params%Tgas    = y(chem_species%nSpecies+1)
  heating_cooling_params%X_H2    = y(chem_idx_some_spe%i_H2)
  heating_cooling_params%X_HI    = y(chem_idx_some_spe%i_HI)
  heating_cooling_params%X_CI    = y(chem_idx_some_spe%i_CI)
  heating_cooling_params%X_Cplus = y(chem_idx_some_spe%i_Cplus)
  heating_cooling_params%X_OI    = y(chem_idx_some_spe%i_OI)
  heating_cooling_params%X_CO    = y(chem_idx_some_spe%i_CO)
  heating_cooling_params%X_H2O   = y(chem_idx_some_spe%i_H2O)
  heating_cooling_params%X_OH    = y(chem_idx_some_spe%i_OH)
  heating_cooling_params%X_E     = y(chem_idx_some_spe%i_E)
  heating_cooling_params%X_Hplus = y(chem_idx_some_spe%i_Hplus)
  heating_cooling_params%X_gH    = y(chem_idx_some_spe%i_gH)
  heating_cooling_params%R_H2_form_rate_coeff = chem_params%R_H2_form_rate_coeff
  if (chem_solver_params%H2_form_use_moeq) then
    heating_cooling_params%R_H2_form_rate = &
      heating_cooling_params%R_H2_form_rate_coeff * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%X_HI * &
      heating_cooling_params%n_gas
  else
    heating_cooling_params%R_H2_form_rate = &
      heating_cooling_params%R_H2_form_rate_coeff * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%X_gH * &
      heating_cooling_params%n_gas
  end if
  hc_Tgas = y(chem_species%nSpecies+1)
  hc_Tdust = heating_cooling_params%Tdust
  r = &
    heating_minus_cooling() * phy_SecondsPerYear / (chem_params%n_gas * phy_kBoltzmann_CGS)
end subroutine realtime_heating_cooling_rate




subroutine chem_ode_jac(NEQ, t, y, j, ian, jan, pdj)
  use chemistry
  use heating_cooling
  use trivials
  implicit none
  double precision t, rtmp
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, i, j, k, i1
  double precision dT_dt_1, dT_dt_2, del_ratio, del_0, delta_y
  double precision, dimension(NEQ) :: ydot1, ydot2
  del_ratio = 1D-3
  del_0 = 1D-16
  pdj = 0D0
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
        if (chem_net%reac_names(1, i) .eq. 'gH') then
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
  end do
  !
  if (chem_solver_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (is_in_list_int(j, chem_idx_some_spe%nItem, chem_idx_some_spe%idx)) then
      call realtime_heating_cooling_rate(dT_dt_1, NEQ, y)
      delta_y = y(j) * del_ratio + del_0
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call realtime_heating_cooling_rate(dT_dt_2, NEQ, y)
      pdj(chem_species%nSpecies+1) = (dT_dt_2 - dT_dt_1) / delta_y
      y(j) = rtmp
    else if (j .eq. (chem_species%nSpecies+1)) then
      call chem_ode_f(NEQ, t, y, ydot1)
      delta_y = y(j) * del_ratio * 1D1 + del_0
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call chem_ode_f(NEQ, t, y, ydot2)
      pdj = (ydot2 - ydot1) / delta_y
      y(j) = rtmp
    end if
  else
    pdj(chem_species%nSpecies+1) = 0D0
  end if
end subroutine chem_ode_jac
