module vertical_structure

use data_struct
use grid
use ray_propagating
use phy_const

implicit none


double precision, dimension(:), allocatable, private :: dz_s

contains


subroutine vertical_pressure_gravity_balance_alt(mstar, useTdust, &
    Tdust_lowerlimit, ngas_lowerlimit, ndust_lowerlimit, fix_dust_struct, &
    disk_gas_mass_preset, &
    maxfac, minfac)
  double precision, intent(in) :: mstar
  logical, intent(in), optional :: useTdust, fix_dust_struct
  double precision, intent(in), optional :: Tdust_lowerlimit, &
    ngas_lowerlimit, ndust_lowerlimit, disk_gas_mass_preset
  double precision, intent(out), optional :: maxfac, minfac
  logical useTd, fix_d
  integer ic, ir
  type(type_cell), pointer :: c1, c2
  double precision Sig0, Sig1, fac1, fac2, fac, fac_ch, r1, r2, z0, z1, z2, T1, T2
  double precision, dimension(MaxNumOfDustComponents) :: SigD0, SigD1, facD
  double precision Td_low, nd_low, ng_low
  double precision, parameter :: min_d2g_ratio=1D-30, max_d2g_ratio=1D-3
  double precision disk_gas_mass_actual, f_resc_global
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
    ng_low = 1D-4
  end if
  if (present(maxfac)) then
    maxfac = 0D0
  end if
  if (present(minfac)) then
    minfac = 1D100
  end if
  !
  if (present(disk_gas_mass_preset)) then
    disk_gas_mass_actual = calc_disk_gas_mass()
    if ((disk_gas_mass_actual .le. 0D0) .or. isnan(disk_gas_mass_actual))  then
      write(*,*) 'In vertical_pressure_gravity_balance_alt:'
      write(*,*) 'Disk mass <= 0 or is nan: ', disk_gas_mass_actual
      call error_stop()
    end if
    f_resc_global = disk_gas_mass_preset / disk_gas_mass_actual
  else
    f_resc_global = 1D0
  end if
  !
#ifdef DIAGNOSIS_TRACK_FUNC_CALL
  write(*,*) 'Running vertical_pressure_gravity_balance_alt'
  write(*,*) 'f_resc_global = ', f_resc_global
#endif
  do ic=1, bott_cells%nlen
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
      fac = min(exp(-fac1-fac2) * T1 / T2, 1D0)
      !
      fac_ch = c1%par%n_gas * fac / (c2%par%n_gas + 1D-100)
      !
      c2%par%n_gas = c1%par%n_gas * fac
      if (.not. fix_d) then
        c2%par%rho_dusts = c1%par%rho_dusts * fac
      end if
      !
      if (present(maxfac) .and. (c1%par%n_gas .ge. ng_low)) then
        maxfac = max(maxfac, fac_ch)
      end if
      if (present(minfac) .and. (c1%par%n_gas .ge. ng_low)) then
        minfac = min(minfac, fac_ch)
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
    fac  = f_resc_global * Sig0/(Sig1+1D-100)
    facD = f_resc_global * SigD0/(SigD1+1D-100)
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
            (c1%par%n_gas*max_d2g_ratio .le. c1%par%ndust_tot) .or. &
            (c1%par%ndust_tot .le. c1%par%n_gas * min_d2g_ratio)) then
          c1%using = .false.
        end if
      end if
    end do
  end do
end subroutine vertical_pressure_gravity_balance_alt



subroutine calc_dustgas_struct_snippet1(c)
  type(type_cell), intent(inout) :: c
  integer i
  c%par%mdust_tot = 0D0
  c%par%ndust_tot = 0D0
  c%par%sigdust_ave = 0D0
  do i=1, c%par%ndustcompo
    c%par%n_dusts(i) = c%par%rho_dusts(i) / c%par%mp_dusts(i)
    c%par%mdusts_cell(i) = c%par%rho_dusts(i) * c%par%volume
    c%par%mdust_tot = c%par%mdust_tot + c%par%mdusts_cell(i)
    c%par%ndust_tot = c%par%ndust_tot + c%par%n_dusts(i)
    c%par%sigdust_ave = c%par%sigdust_ave + c%par%n_dusts(i) * c%par%sig_dusts(i)
  end do
  if (c%par%ndust_tot .le. 1D-100) then
    c%par%sigdust_ave = 0D0
  else
    c%par%sigdust_ave = c%par%sigdust_ave / c%par%ndust_tot
  end if
  !
  c%par%GrainRadius_CGS = sqrt(c%par%sigdust_ave / phy_Pi)
  c%par%SitesPerGrain = 4D0 * c%par%sigdust_ave * const_SitesDensity_CGS
end subroutine calc_dustgas_struct_snippet1


subroutine calc_dustgas_struct_snippet2(c)
  type(type_cell), intent(inout) :: c
  c%par%mgas_cell = c%par%n_gas * c%par%volume * &
                    (phy_mProton_CGS * c%par%MeanMolWeight)
  !
  c%par%ratioDust2GasMass = c%par%mdust_tot / c%par%mgas_cell
  c%par%ratioDust2HnucNum = c%par%ndust_tot / c%par%n_gas
  c%par%dust_depletion = c%par%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
  !
  c%par%abso_wei = c%par%n_dusts*c%par%sig_dusts / &
                   (sum(c%par%n_dusts*c%par%sig_dusts) + 1D-100)
  c%val(1) = c%par%n_gas
end subroutine calc_dustgas_struct_snippet2


function calc_disk_gas_mass() result(m)
  double precision m, vol
  type(type_cell), pointer :: c
  integer ic
  m = 0D0
  do ic=1, leaves%nlen
    c => leaves%list(ic)%p
    if (c%using) then
      vol = phy_Pi * (c%xmax + c%xmin) * &
            (c%xmax - c%xmin) * (c%ymax-c%ymin) * phy_AU2cm**3
      m = m + vol * c%par%n_gas * &
          (phy_mProton_CGS * c%par%MeanMolWeight)
    end if
  end do
  m = m * 2D0 / phy_Msun_CGS  ! Two sides of the disk.
end function calc_disk_gas_mass


subroutine get_ndiv(c, n_div)
  type(type_cell), pointer, intent(in) :: c
  integer, intent(out) :: n_div
  !
  double precision, parameter :: mindens_refine = 1D3
  double precision, parameter :: minTdust_refine = 5.0D0
  double precision, parameter :: minG0_UV_refine = 1.0D-3
  double precision, parameter :: maxdz_ratio = 0.1D0
  double precision, parameter :: mindz_ratio = 0.005D0
  double precision, parameter :: dTdust_ratio_max = 0.3D0
  double precision, parameter :: dens_ratio_max = 1.5D0
  integer, parameter :: max_n_div = 10
  !
  double precision :: maxdens, mindens, maxTdust, minTdust
  double precision maxdz_here, mindz_here
  integer i, i0
  !
  n_div = 0
  !
  if (.not. c%using) then
    return
  end if
  !
  if ((c%ymax-c%ymin) .le. grid_config%smallest_cell_size) then
    return
  end if
  !
  if (c%par%n_gas .le. mindens_refine) then
    return
  end if
  !
  if (c%par%Tdust .le. minTdust_refine) then
    return
  end if
  !
  if (c%par%G0_UV_toStar_photoDesorb .le. minG0_UV_refine) then
    return
  end if
  !
  maxdens = c%par%n_gas
  mindens = c%par%n_gas
  maxTdust = c%par%Tdust
  minTdust = c%par%Tdust
  do i=1, c%above%n
    i0 = c%above%idx(i)
    if (leaves%list(i0)%p%using) then
      maxdens = max(maxdens, leaves%list(i0)%p%par%n_gas)
      mindens = min(mindens, leaves%list(i0)%p%par%n_gas)
      maxTdust = max(maxTdust, leaves%list(i0)%p%par%Tdust)
      minTdust = min(minTdust, leaves%list(i0)%p%par%Tdust)
    end if
  end do
  do i=1, c%below%n
    i0 = c%below%idx(i)
    if (leaves%list(i0)%p%using) then
      maxdens = max(maxdens, leaves%list(i0)%p%par%n_gas)
      mindens = min(mindens, leaves%list(i0)%p%par%n_gas)
      maxTdust = max(maxTdust, leaves%list(i0)%p%par%Tdust)
      minTdust = min(minTdust, leaves%list(i0)%p%par%Tdust)
    end if
  end do
  !
  minTdust = max(minTdust, minTdust_refine)
  mindens = max(mindens, mindens_refine)
  !
  maxdz_here = (c%xmax+c%xmin) * 0.5D0 * maxdz_ratio
  mindz_here = (c%xmax+c%xmin) * 0.5D0 * mindz_ratio
  !
  if (((c%par%Av_toISM .ge. 1D-3) .and. (c%par%Av_toISM .le. 2D0)) .or. &
      ((c%par%Av_toStar .ge. 1D-3) .and. (c%par%Av_toStar .le. 2D0))) then
    maxdz_here = min(maxdz_here, &
                     2D-1 * min(c%par%Av_toISM, c%par%Av_toStar) &
                     / (c%par%sigdust_ave * c%par%ndust_tot) &
                     / phy_AU2cm)
  end if
  !
  n_div = ceiling((c%ymax - c%ymin) / maxdz_here)
  !
  n_div = max(n_div, &
              ceiling((maxTdust/c%par%Tdust - 1D0) / dTdust_ratio_max))
  n_div = max(n_div, &
              ceiling((c%par%Tdust/minTdust - 1D0) / dTdust_ratio_max))
  n_div = max(n_div, ceiling(log(maxdens/c%par%n_gas) / log(dens_ratio_max)))
  n_div = max(n_div, ceiling(log(c%par%n_gas/mindens) / log(dens_ratio_max)))
  !
  if      ((n_div .le. 1) .and. (associated(c%inner))) then
    if (c%inner%n .ge. 4) then
      n_div = c%inner%n / 2
    end if
  else if ((n_div .le. 1) .and. (associated(c%outer))) then
    if (c%outer%n .ge. 4) then
      n_div = c%outer%n / 2
    end if
  end if
  !
  n_div = min(n_div, ceiling((c%ymax - c%ymin) / mindz_here), max_n_div)
end subroutine get_ndiv


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
  allocate(dz_s(leaves%nlen))
  !
  do i=1, leaves%nlen
    dz_s(i) = leaves%list(i)%p%ymax - leaves%list(i)%p%ymin
  end do
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
      dz_s(i)         = dz_s(i)         / frescale
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
  deallocate(dz_s)
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
        cthis%ymax = ybelow + dz_s(ic)
        !
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
        !
        cthis%par%n_gas     = cthis%par%n_gas   * frescale
        cthis%par%n_dusts   = cthis%par%n_dusts * frescale
        cthis%par%ndust_tot = cthis%par%ndust_tot * frescale
        cthis%par%rho_dusts = cthis%par%rho_dusts * frescale
        cthis%par%volume    = cthis%par%volume / frescale
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

end module vertical_structure
