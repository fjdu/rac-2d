module vertical

use data_struct
use grid
use chemistry


implicit none

contains

subroutine vertical_pressure_gravity_balance
  integer i
  double precision dznew, pnew, nnew, vnew
  ! For displaying some text to the screen
  !
  if (.not. grid_config%columnwise) then
    write(*, '(A)') &
      'It is difficult to do vertical mechanical balance if not in column mode.'
    return
  end if
  !
  ! Calculate the new density and size of each cell
  do i=1, leaves%nlen
    pnew = -leaves%list(i)%p%par%gravity_acc_z / leaves%list(i)%p%par%area_T
    ! Not change too much in one go
    pnew = sqrt(pnew * leaves%list(i)%p%par%pressure_thermal)
    !
    nnew = pnew / &
      (leaves%list(i)%p%par%Tgas * phy_kBoltzmann_CGS * &
       (leaves%list(i)%p%abundances(chem_idx_some_spe%i_HI) + &
        leaves%list(i)%p%abundances(chem_idx_some_spe%i_H2)))
    vnew = leaves%list(i)%p%par%mgas_cell / nnew / &
      (phy_mProton_CGS * leaves%list(i)%p%par%MeanMolWeight)
    dznew = vnew / leaves%list(i)%p%par%area_T / phy_AU2cm
    if (isnan(dznew)) then
      !write(str_disp, '(A)') '! dznew is NaN!'
      !call display_string_both(str_disp, a_book_keeping%fU)
      !write(str_disp, '(A, 3ES16.9)') 'pnew,nnew,vnew:', pnew, nnew, vnew
      !call display_string_both(str_disp, a_book_keeping%fU)
      cycle
    end if
    !
    leaves%list(i)%p%par%n_dusts = leaves%list(i)%p%par%n_dusts * nnew / &
                                   leaves%list(i)%p%par%n_gas
    leaves%list(i)%p%par%n_gas = nnew
    leaves%list(i)%p%val(1) = nnew
    leaves%list(i)%p%par%volume = vnew
    leaves%list(i)%p%par%dz = dznew
    leaves%list(i)%p%par%pressure_thermal = pnew
  end do
  !
  ! Move the cells column by column from bottom to top.
  call shift_and_scale_above
  !
  ! Remake neigbor information
  call grid_make_neighbors
  !
  call grid_make_surf_bott
  !
end subroutine vertical_pressure_gravity_balance



subroutine shift_and_scale_above
  type(type_cell), pointer :: cthis
  double precision ybelow
  integer i, ic
  !
  do ic=1, bott_cells%nlen ! Loop over the columns
    ybelow = columns(ic)%list(1)%p%ymin
    do i=1, columns(ic)%nlen
      cthis => columns(ic)%list(i)%p
      if (cthis%using) then
        cthis%ymin = ybelow
        cthis%ymax = ybelow + cthis%par%dz
        !
        cthis%par%zmin = cthis%ymin
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        ! dz is already set
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
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
    cthis => columns(ic)%list(i)%p
    if (cthis%ymax .le. root%ymax) then
      ! Rescale the density
      cthis%val(1) = cthis%val(1) * &
        (cthis%ymax - cthis%ymin) / (root%ymax - cthis%ymin)
      cthis%ymax = root%ymax
      !
      if (cthis%using) then
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        cthis%par%dz = cthis%par%zmax - cthis%par%zmin
        !
        cthis%par%n_dusts = cthis%par%n_dusts * cthis%val(1) / cthis%par%n_gas
        cthis%par%n_gas  = cthis%val(1)
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
      end if
    else
      write(*,'(A)') 'Should not have this case:'
      write(*,'(A/)') 'in shift_and_scale_above.'
    end if
  end do
end subroutine shift_and_scale_above


end module vertical
