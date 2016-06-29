module cont_lookuptable

use grid
use montecarlo

implicit none

contains

subroutine allocate_local_cont_lut(c)
  type(type_cell), intent(inout), pointer :: c
  !
  if (.not. allocated(c%cont_lut)) then
    allocate(c%cont_lut)
  end if
  !
  if (.not. allocated(c%cont_lut%J)) then
    c%cont_lut%n = dust_0%n
    allocate(c%cont_lut%J(dust_0%n))
  end if
end subroutine allocate_local_cont_lut



subroutine make_local_cont_lut(c)
  type(type_cell), intent(inout), pointer :: c
  integer i
  double precision dlam, lam
  !
  do i=1, c%cont_lut%n
    if (i .lt. c%cont_lut%n) then
      dlam = dust_0%lam(i+1) - dust_0%lam(i)
      lam = (dust_0%lam(i+1) + dust_0%lam(i)) * 0.5D0
      ! Energy per unit area per unit frequency per second per sqradian
      c%cont_lut%J(i) = c%optical%flux(i) &
        / dlam * lam * lam * phy_Angstrom2cm / phy_SpeedOfLight_CGS &
        / (4D0 * phy_Pi)
    else
      c%cont_lut%J(i) = c%cont_lut%J(i-1)
    end if
  end do
end subroutine make_local_cont_lut




end module cont_lookuptable
