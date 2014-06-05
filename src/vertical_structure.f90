module vertical

use data_struct
use grid
use ray_propagating
use phy_const

implicit none

contains


subroutine vertical_pressure_gravity_balance_alt(mstar, useTdust, &
    Tdust_lowerlimit, ngas_lowerlimit, ndust_lowerlimit, fix_dust_struct)
  double precision, intent(in) :: mstar
  logical, intent(in), optional :: useTdust, fix_dust_struct
  double precision, intent(in), optional :: Tdust_lowerlimit, &
    ngas_lowerlimit, ndust_lowerlimit
  logical useTd, fix_d
  integer ic, ir
  type(type_cell), pointer :: c1, c2
  double precision Sig0, Sig1, fac1, fac2, fac, r1, r2, z0, z1, z2, T1, T2
  double precision, dimension(MaxNumOfDustComponents) :: SigD0, SigD1, facD
  double precision minfac, maxfac, mulfac
  double precision Td_low, nd_low, ng_low
  double precision, parameter :: min_d2g_ratio=1D-20, max_d2g_ratio=1D-5
  !
  if (present(useTdust)) then
    useTd = useTdust
    if (useTd) then
      if (present(Tdust_lowerlimit)) then
        Td_low = Tdust_lowerlimit
      else
        Td_low = 5D0
      end if
    end if
  else
    useTd = .false.
  end if
  if (present(fix_dust_struct)) then
    fix_d = fix_dust_struct
  else
    fix_d = .true.
  end if
  if (present(ndust_lowerlimit)) then
    nd_low = ndust_lowerlimit
  else
    nd_low = 1D-20
  end if
  if (present(ngas_lowerlimit)) then
    ng_low = ngas_lowerlimit
  else
    ng_low = 1D-2
  end if
  !
  do ic=1, bott_cells%nlen
    minfac = 1D100
    maxfac = 0D0
    mulfac = 1D0
    Sig0 = 0D0
    SigD0 = 0D0
    do ir=1, columns(ic)%nlen
      c1 => columns(ic)%list(ir)%p
      if (c1%using) then
        Sig0 = Sig0 + (c1%ymax - c1%ymin) * &
          c1%par%n_gas*c1%par%MeanMolWeight*phy_mProton_CGS
        SigD0 = SigD0 + (c1%ymax - c1%ymin) * c1%par%rho_dusts
      end if
    end do
    !
    do ir=2, columns(ic)%nlen
      c1 => columns(ic)%list(ir-1)%p
      c2 => columns(ic)%list(ir)%p
      if (.not. c2%using) then
        exit
      end if
      r1 = sqrt((c1%xmin+c1%xmax)**2 + (c1%ymin+c1%ymax)**2) * 0.5D0 * phy_AU2cm
      r2 = sqrt((c2%xmin+c2%xmax)**2 + (c2%ymin+c2%ymax)**2) * 0.5D0 * phy_AU2cm
      z0 = 0.5D0 * (c1%ymax + c1%ymin) * phy_AU2cm
      z1 = c1%ymax * phy_AU2cm
      z2 = 0.5D0 * (c2%ymax + c2%ymin) * phy_AU2cm
      if (useTd) then
        T1 = c1%par%Tdust
        T2 = c2%par%Tdust
        if ((T1 .le. Td_low) .or. (T2 .le. Td_low)) then
          cycle
        end if
      else
        T1 = c1%par%Tgas
        T2 = c2%par%Tgas
      end if
      fac1 = phy_GravitationConst_CGS * mstar * phy_Msun_CGS &
            * c1%par%MeanMolWeight * phy_mProton_CGS &
            / 2D0 / r1**3 / (phy_kBoltzmann_CGS * T1) &
            * (z1-z0)*(z1+z0)
      fac2 = phy_GravitationConst_CGS * mstar * phy_Msun_CGS &
            * c2%par%MeanMolWeight * phy_mProton_CGS &
            / 2D0 / r2**3 / (phy_kBoltzmann_CGS * T2) &
            * (z2-z1)*(z2+z1)
      !
      fac = exp(-fac1-fac2) * T1 / T2
      minfac = min(minfac, fac)
      maxfac = max(maxfac, fac)
      mulfac = mulfac * fac
      !
      !write(*, '(3ES16.6)') fac1, fac2, fac
      !
      c2%par%n_gas = c1%par%n_gas * fac
      if (.not. fix_d) then
        c2%par%rho_dusts = c1%par%rho_dusts * fac
      end if
    end do
    !
    Sig1 = 0D0
    SigD1 = 0D0
    do ir=1, columns(ic)%nlen
      c1 => columns(ic)%list(ir)%p
      if (c1%using) then
        Sig1 = Sig1 + (c1%ymax - c1%ymin) * &
          c1%par%n_gas*c1%par%MeanMolWeight*phy_mProton_CGS
        SigD1 = SigD1 + (c1%ymax - c1%ymin) * c1%par%rho_dusts
      end if
    end do
    !
    fac = Sig0/(Sig1+1D-100)
    facD = SigD0/(SigD1+1D-100)
    !
    do ir=1, columns(ic)%nlen
      c1 => columns(ic)%list(ir)%p
      if (c1%using) then
        !
        c1%par%n_gas = c1%par%n_gas * fac
        !
        if (.not. fix_d) then
          !
          c1%par%rho_dusts = c1%par%rho_dusts * facD
          !
          call calc_dustgas_struct_snippet1(c1)
          !
        end if
        !
        call calc_dustgas_struct_snippet2(c1)
        !
        if ((c1%par%ndust_tot .le. nd_low) .or. &
            (c1%par%n_gas .le. ng_low) .or. &
            (c1%par%n_gas .le. c1%par%ndust_tot/max_d2g_ratio) .or. &
            (c1%par%ndust_tot .le. c1%par%n_gas * min_d2g_ratio)) then
          c1%using = .false.
        end if
      end if
    end do
    !
    write(*, '(A, I6, 4X, A, 4ES16.6)') 'Column: ', ic, &
        'min,max,mul,fac: ', minfac, maxfac, mulfac, fac
  end do
end subroutine vertical_pressure_gravity_balance_alt


subroutine calc_dustgas_struct_snippet1(c)
  type(type_cell), intent(inout) :: c
  c%par%n_dusts = c%par%rho_dusts / (c%par%mp_dusts + 1D-100)
  c%par%mdusts_cell = c%par%rho_dusts * c%par%volume
  c%par%mdust_tot = sum(c%par%mdusts_cell)
  c%par%ndust_tot = sum(c%par%n_dusts)
  c%par%sigdust_ave = sum(c%par%n_dusts * c%par%sig_dusts) / &
                      (c%par%ndust_tot+1D-100)
  !
  c%par%GrainRadius_CGS = sqrt(c%par%sigdust_ave / phy_Pi)
  c%par%SitesPerGrain = 4D0 * c%par%sigdust_ave * const_SitesDensity_CGS
end subroutine calc_dustgas_struct_snippet1


subroutine calc_dustgas_struct_snippet2(c)
  type(type_cell), intent(inout) :: c
  c%par%mgas_cell = c%par%n_gas * c%par%volume * &
                    (phy_mProton_CGS * c%par%MeanMolWeight)
  !
  c%par%ratioDust2GasMass = c%par%mdust_tot / (c%par%mgas_cell + 1D-100)
  c%par%ratioDust2HnucNum = c%par%ndust_tot / (c%par%n_gas + 1D-100)
  c%par%dust_depletion = c%par%ratioDust2GasMass / (phy_ratioDust2GasMass_ISM + 1D-100)
  !
  c%par%abso_wei = c%par%n_dusts*c%par%sig_dusts / &
                   (sum(c%par%n_dusts*c%par%sig_dusts) + 1D-100)
  c%val(1) = c%par%n_gas
end subroutine calc_dustgas_struct_snippet2


subroutine vertical_pressure_gravity_balance(frescale_max, frescale_min, useTdust, max_dz)
  integer i
  double precision dznew, pold, pnew
  double precision, intent(out), optional :: frescale_max, frescale_min
  double precision, intent(in),optional :: max_dz
  logical, intent(in), optional :: useTdust
  logical useTd
  double precision :: frescale, fr_max, fr_min, maxdz
  !
  if (present(useTdust)) then
    useTd = useTdust
  else
    useTd = .false.
  end if
  !
  if (.not. grid_config%columnwise) then
    write(*, '(A)') &
      'At the present I will not do vertical structure calculation if I am not in column mode.'
    return
  end if
  !
  fr_max = 0D0
  fr_min = 1D99
  !
  ! Calculate the new density and size of each cell
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (useTd) then
        if (c%par%Tdust .le. 5D0) then
          cycle
        end if
        c%par%pressure_thermal = &
          c%par%n_gas * c%par%Tdust * phy_kBoltzmann_CGS
      else
        c%par%pressure_thermal = &
          c%par%n_gas * c%par%Tgas  * phy_kBoltzmann_CGS
      end if
      !
      pold = c%par%pressure_thermal
      !
      pnew = -c%par%gravity_acc_z / c%par%area_T
      !
      ! Not change too much in one go if use gas temperature to calculate pressure.
      pnew = sqrt(sqrt(pnew**3 * pold))
      !
      pnew = max(min(pnew, pold*1D2), pold*1D-2)
      if (present(max_dz)) then
        maxdz = max_dz
      else
        maxdz = 0.25D0 * (c%xmin + c%xmax + c%ymin + c%ymax) + 1D0 * root%ymax
      end if
      !
      frescale = max(pnew / pold, (c%ymax-c%ymin)/maxdz)
      !
      fr_max = max(frescale, fr_max)
      fr_min = min(frescale, fr_min)
      !
      c%par%n_gas     = c%par%n_gas     * frescale
      c%par%n_dusts   = c%par%n_dusts   * frescale
      c%par%ndust_tot = c%par%ndust_tot * frescale
      c%par%rho_dusts = c%par%rho_dusts * frescale
      c%par%volume    = c%par%volume    / frescale
      c%par%dz        = c%par%dz        / frescale
      c%par%area_I    = c%par%area_I    / frescale
      c%par%area_O    = c%par%area_O    / frescale
      c%val(1) = c%par%n_gas
      c%par%pressure_thermal = pnew
      !
      c%par%surf_area = c%par%area_T + c%par%area_B + c%par%area_I + c%par%area_O
    end associate
  end do
  !
  if (present(frescale_max)) then
    frescale_max = fr_max
  end if
  if (present(frescale_min)) then
    frescale_min = fr_min
  end if
  !
  ! Move the cells column by column from bottom to top.
  call shift_and_scale_above
  !
  ! Reset the bounds of all the background cells.
  call reset_cell_bounds(root)
  !
  call clean_peculiar_cells
  !
end subroutine vertical_pressure_gravity_balance


subroutine clean_peculiar_cells
  integer i
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (c%ymax .ge. grid_config%zmax*1D1) then
        c%using = .false.
      end if
    end associate
  end do
end subroutine clean_peculiar_cells


subroutine shift_and_scale_above
  type(type_cell), pointer :: cthis
  double precision ybelow, frescale
  integer i, ic
  !
  do ic=1, bott_cells%nlen ! Loop over the columns
    ybelow = columns(ic)%list(1)%p%ymin
    do i=1, columns(ic)%nlen
      cthis => columns(ic)%list(i)%p
      if (associated(cthis%par)) then
        cthis%ymin = ybelow
        cthis%ymax = ybelow + cthis%par%dz
        !
        cthis%par%zmin = cthis%ymin
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        ! dz is already set
      else
        cthis%ymax = ybelow + (cthis%ymax - cthis%ymin)
        cthis%ymin = ybelow
      end if
      !
      root%ymax = max(root%ymax, cthis%ymax)
      !
      ybelow = cthis%ymax
      !
    end do
  end do
  !
  ! Align the top edge to the domain upper boundary
  do ic=1, bott_cells%nlen
    i = columns(ic)%nlen
    ! The top cell of a column
    cthis => columns(ic)%list(i)%p
    if (cthis%ymax .le. root%ymax) then
      ! Rescale the density
      frescale = (cthis%ymax - cthis%ymin) / (root%ymax - cthis%ymin)
      cthis%val(1) = cthis%val(1) * frescale
      cthis%ymax = root%ymax
      !
      if (associated(cthis%par)) then
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        !
        cthis%par%n_gas     = cthis%par%n_gas   * frescale
        cthis%par%n_dusts   = cthis%par%n_dusts * frescale
        cthis%par%ndust_tot = cthis%par%ndust_tot * frescale
        cthis%par%rho_dusts = cthis%par%rho_dusts * frescale
        cthis%par%volume    = cthis%par%volume / frescale
        cthis%par%dz        = cthis%par%dz     / frescale
        !
        cthis%par%area_I = cthis%par%area_I / frescale
        cthis%par%area_O = cthis%par%area_O / frescale
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
      end if
    else
      write(*,'(A)') 'Should not have this case:'
      write(*,'(A/)') 'in shift_and_scale_above.'
      write(*,'(A, 2ES16.6)') 'ymax of cthis and root:', cthis%ymax, root%ymax
    end if
  end do
end subroutine shift_and_scale_above



subroutine vertical_pressure_gravity_balance_simple
  integer i
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      if (c%using) then
        c%par%n_gas = c%par%n_gas * &
          abs(c%par%gravity_acc_z / c%par%area_T) / c%par%pressure_thermal
      end if
    end associate
  end do
end subroutine vertical_pressure_gravity_balance_simple

end module vertical
