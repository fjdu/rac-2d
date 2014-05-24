module vertical

use data_struct
use grid
use ray_propagating
use phy_const


implicit none

contains

subroutine vertical_pressure_gravity_balance(frescale_max, frescale_min, useTdust)
  integer i
  double precision dznew, pold, pnew
  double precision, intent(out), optional :: frescale_max, frescale_min
  logical, intent(in), optional :: useTdust
  logical useTd
  double precision :: frescale, fr_max, fr_min
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
      if (.not. useTd) then
        ! Not change too much in one go if use gas temperature to calculate pressure.
        pnew = sqrt(sqrt(pnew**3 * pold))
      end if
      !
      pnew = max(min(pnew, pold*1D2), pold*1D-2)
      !
      frescale = pnew / pold
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


end module vertical
