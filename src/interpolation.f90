module spline_1d_2d

implicit none


type :: type_spline_1D
  integer n
  integer :: itype=0
  double precision, dimension(:), allocatable :: xi, yi, ddyi
end type type_spline_1D


type :: type_spline_2d
  integer :: nx, ny, itype=0
  double precision, dimension(:), allocatable :: xi, yi
  double precision, dimension(:,:), allocatable :: vi
  type(type_spline_1D), dimension(:), allocatable :: sp1d
  type(type_spline_1D) sp1d_using
end type type_spline_2d


contains


function spline2d_interpol(x, y, sp2d, extrapolate, stat1, stat2)
  double precision spline2d_interpol
  double precision, intent(in) :: x, y
  type(type_spline_2d), intent(inout) :: sp2d
  logical, intent(in), optional :: extrapolate
  integer, intent(out), optional :: stat1, stat2
  logical extrap
  integer j, sta1, sta2
  !if ((x .lt. sp2d%xi(1)) .or. &
  !    (x .gt. sp2d%xi(sp2d%nx)) .or. &
  !    (y .lt. sp2d%yi(1)) .or. &
  !    (y .gt. sp2d%yi(sp2d%ny))) then
  !  spline2d_interpol = 0D0
  !  return
  !end if
  if (present(extrapolate)) then
    extrap = extrapolate
  else
    extrap = .false.
  end if
  associate(sp => sp2d%sp1d_using)
    do j=1, sp2d%nx
      sp%yi(j) = spline1d_interpol(y, sp2d%sp1d(j), extrapolate=extrap, stat=sta2)
    end do
    call spline1d_prepare(sp)
    spline2d_interpol = spline1d_interpol(x, sp, extrapolate=extrap, stat=sta1)
  end associate
  if (present(stat1)) then
    stat1 = sta1
  end if
  if (present(stat2)) then
    stat2 = sta2
  end if
end function spline2d_interpol


pure subroutine spline2d_prepare(sp2d)
  type(type_spline_2d), intent(inout) :: sp2d
  integer i
  allocate(sp2d%sp1d_using%xi(sp2d%nx), &
           sp2d%sp1d_using%yi(sp2d%nx), &
           sp2d%sp1d_using%ddyi(sp2d%nx))
  sp2d%sp1d_using%n = sp2d%nx
  sp2d%sp1d_using%itype = sp2d%itype
  sp2d%sp1d_using%xi = sp2d%xi
  allocate(sp2d%sp1d(sp2d%nx))
  do i=1, sp2d%nx
    associate(sp => sp2d%sp1d(i))
      sp%n = sp2d%ny
      sp%itype = sp2d%itype
      allocate(sp%xi(sp2d%ny), &
               sp%yi(sp2d%ny), &
               sp%ddyi(sp2d%ny))
      sp%xi = sp2d%yi
      sp%yi = sp2d%vi(i, :)
      call spline1d_prepare(sp)
    end associate
  end do
end subroutine spline2d_prepare


pure subroutine spline1d_prepare(sp1d)
  ! Not-a-knot end condition is taken from:
  ! http://www.mathworks.com/matlabcentral/newsreader/view_thread/172988
  ! It makes the third derives to be continuous at the left and right bars.
  type(type_spline_1D), intent(inout) :: sp1d
  double precision, dimension(:), allocatable :: a, b, c, r
  integer j
  sp1d%ddyi(1) = 0D0
  sp1d%ddyi(sp1d%n) = 0D0
  allocate( a(sp1d%n-2), &
            b(sp1d%n-2), &
            c(sp1d%n-2), &
            r(sp1d%n-2))
  do j=2, sp1d%n-1
    c(j-1) = (sp1d%xi(j+1) - sp1d%xi(j)) / 6D0
    b(j-1) = (sp1d%xi(j+1) - sp1d%xi(j-1)) / 3D0
    a(j-1) = (sp1d%xi(j) - sp1d%xi(j-1)) / 6D0
    r(j-1) = (sp1d%yi(j+1) - sp1d%yi(j)) / (sp1d%xi(j+1) - sp1d%xi(j)) - &
             (sp1d%yi(j) - sp1d%yi(j-1)) / (sp1d%xi(j) - sp1d%xi(j-1))
  end do
  !
  if (sp1d%itype .eq. 0) then ! Linear interpolation
    sp1d%ddyi = 0D0
  else if (sp1d%itype .eq. 1) then ! Set dy/dx = dely/delx at the boundary
    b(1) = b(1) - (sp1d%xi(2) - sp1d%xi(1)) / 12D0
    b(sp1d%n-2) = b(sp1d%n-2) - (sp1d%xi(sp1d%n) - sp1d%xi(sp1d%n-1)) / 12D0
    call tridiagonal_solve(sp1d%n-2, a, b, c, r, sp1d%ddyi(2:(sp1d%n-1)))
    sp1d%ddyi(1) = -0.5D0 * sp1d%ddyi(2)
    sp1d%ddyi(sp1d%n) = -0.5D0 * sp1d%ddyi(sp1d%n-1)
  else if (sp1d%itype .eq. 2) then ! Set ddy/ddx = 0 at the boundary
    call tridiagonal_solve(sp1d%n-2, a, b, c, r, sp1d%ddyi(2:(sp1d%n-1)))
  else if (sp1d%itype .eq. 3) then ! Not-a-knot boundary condition
    associate(&
      x1 => sp1d%xi(1), x2 => sp1d%xi(2), x3 => sp1d%xi(3), &
      xn => sp1d%xi(sp1d%n), xn1 => sp1d%xi(sp1d%n-1), xn2 => sp1d%xi(sp1d%n-2))
      b(1) = b(1) + (x2-x1)/6D0 * (1D0 + (x2-x1)/(x3-x2))
      c(1) = c(1) - (x2-x1)/6D0 * (x2-x1)/(x3-x2)
      b(sp1d%n-2) = b(sp1d%n-2) + (xn-xn1)/6D0 * (1D0 + (xn-xn1)/(xn1-xn2))
      a(sp1d%n-2) = a(sp1d%n-2) - (xn-xn1)/6D0 * (xn-xn1)/(xn1-xn2)
      call tridiagonal_solve(sp1d%n-2, a, b, c, r, sp1d%ddyi(2:(sp1d%n-1)))
      sp1d%ddyi(1) = sp1d%ddyi(2) - (x2-x1)/(x3-x2) * (sp1d%ddyi(3) - sp1d%ddyi(2))
      sp1d%ddyi(sp1d%n) = sp1d%ddyi(sp1d%n-1) + (xn-xn1)/(xn1-xn2)*(sp1d%ddyi(sp1d%n-1)-sp1d%ddyi(sp1d%n-2))
    end associate
  end if
  deallocate(a, b, c, r)
end subroutine spline1d_prepare


function spline1d_interpol(x, sp1d, extrapolate, stat)
  double precision spline1d_interpol
  double precision, intent(in) :: x
  type(type_spline_1D), intent(in) :: sp1d
  logical, intent(in), optional :: extrapolate
  integer, intent(out), optional :: stat
  double precision A, B, C, D
  logical extrap
  integer j, sta
  if (present(extrapolate)) then
    extrap = extrapolate
  else
    extrap = .false.
  end if
  if (sp1d%xi(1) .gt. x) then
    sta = -1
    if (.not. extrap) then
      spline1d_interpol = sp1d%yi(1)
      if (present(stat)) then
        stat = sta
      end if
      !
      return
      !
    end if
    j = 1
  else if (sp1d%xi(sp1d%n) .lt. x) then
    sta =  1
    if (.not. extrap) then
      spline1d_interpol = sp1d%yi(sp1d%n)
      if (present(stat)) then
        stat = sta
      end if
      !
      return
      !
    end if
    j = sp1d%n - 1
  else
    sta = 0
    do j=1, sp1d%n-1
      if ((sp1d%xi(j) .le. x) .and. (sp1d%xi(j+1) .ge. x)) then
        exit
      end if
    end do
  end if
  A = (sp1d%xi(j+1) - x) / (sp1d%xi(j+1) - sp1d%xi(j))
  B = 1D0 - A
  C = (A**3 - A) * (sp1d%xi(j+1) - sp1d%xi(j))**2 / 6D0
  D = (B**3 - B) * (sp1d%xi(j+1) - sp1d%xi(j))**2 / 6D0
  spline1d_interpol = &
    A * sp1d%yi(j)  + B * sp1d%yi(j+1) + &
    C * sp1d%ddyi(j) + D * sp1d%ddyi(j+1)
  if (present(stat)) then
    stat = sta
  end if
end function spline1d_interpol


pure subroutine tridiagonal_solve(n, a, b, c, r, u, msg)
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: a, b, c, r
  double precision, dimension(n), intent(out) :: u
  double precision, dimension(n) :: gam
  integer, intent(out), optional :: msg
  double precision bet
  integer j
  if (present(msg)) then
    msg = 0
  end if
  if (b(1) .eq. 0D0) then
    !write(*,*) 'Error 1!  In tridiagonal_solve.'
    if (present(msg)) then
      msg = 1
    end if
  end if
  u(1) = r(1) / b(1)
  bet = b(1)
  do j=2, n
    gam(j) = c(j-1) / bet
    bet = b(j) - a(j) * gam(j)
    if (bet .eq. 0D0) then
      !write(*,*) 'Error 2!  In tridiagonal_solve.'
      if (present(msg)) then
        msg = msg + 2
      end if
    end if
    u(j) = (r(j) - a(j)*u(j-1)) / bet
  end do
  do j=n-1, 1, -1
    u(j) = u(j) - gam(j+1) * u(j+1)
  end do
end subroutine tridiagonal_solve

end module spline_1d_2d




module barycentric_1d_2d

implicit none

type :: type_barycentric_1d
  integer n, d
  double precision, dimension(:), allocatable :: xi, yi, wi
end type type_barycentric_1d


type :: type_barycentric_2d
  integer :: nx, ny, d
  double precision, dimension(:), allocatable :: xi, yi
  double precision, dimension(:,:), allocatable :: vi
  type(type_barycentric_1d), dimension(:), allocatable :: bary1d
  type(type_barycentric_1d) bary1d_using
end type type_barycentric_2d


contains


subroutine barycentric2d_prepare(bary2d)
  type(type_barycentric_2d) :: bary2d
  integer i
  bary2d%bary1d_using%n = bary2d%nx
  bary2d%bary1d_using%d = bary2d%d
  allocate(bary2d%bary1d_using%xi(bary2d%nx), &
           bary2d%bary1d_using%yi(bary2d%nx), &
           bary2d%bary1d_using%wi(bary2d%nx))
  bary2d%bary1d_using%xi = bary2d%xi
  allocate(bary2d%bary1d(bary2d%nx))
  do i=1, bary2d%nx
    associate(br => bary2d%bary1d(i))
      br%n = bary2d%ny
      br%d = bary2d%d
      allocate(br%xi(bary2d%ny), &
               br%yi(bary2d%ny), &
               br%wi(bary2d%ny))
      br%xi = bary2d%yi
      br%yi = bary2d%vi(i, :)
      call barycentric1d_prepare(br)
    end associate
  end do
end subroutine barycentric2d_prepare


function barycentric2d_interpol(x, y, bary2d)
  double precision barycentric2d_interpol
  double precision, intent(in) :: x, y
  type(type_barycentric_2d), intent(inout) :: bary2d
  integer j
  associate(br => bary2d%bary1d_using)
    do j=1, bary2d%nx
      br%yi(j) = barycentric1d_interpol(y, bary2d%bary1d(j))
    end do
    call barycentric1d_prepare(br)
    barycentric2d_interpol = barycentric1d_interpol(x, br)
  end associate
end function barycentric2d_interpol


subroutine barycentric1d_prepare(a_barycentric_1d)
  type(type_barycentric_1d), intent(inout) :: a_barycentric_1d
  integer i, j, k, imin
  double precision :: s, p
  if (a_barycentric_1d%n .le. a_barycentric_1d%d) then
    write(*,*) 'd too large.'
    stop
  end if
  if (a_barycentric_1d%d .eq. -1) then
    return
  end if
  if (a_barycentric_1d%d .eq. 0) then
    do i=1, a_barycentric_1d%n
      a_barycentric_1d%wi(i) = dble((-1)**i)
    end do
    return
  end if
  s = 1D0
  do i=1, a_barycentric_1d%n
    a_barycentric_1d%wi(i) = 0D0
    do j=max(i-a_barycentric_1d%d, 1), min(i, a_barycentric_1d%n-a_barycentric_1d%d)
      p = 1D0
      do k=j, j+a_barycentric_1d%d
        if (k .ne. i) then
          p = p * abs(a_barycentric_1d%xi(i) - a_barycentric_1d%xi(k))
        end if
        a_barycentric_1d%wi(i) = a_barycentric_1d%wi(i) + p
      end do
    end do
    a_barycentric_1d%wi(i) = a_barycentric_1d%wi(i) * s
    s = -s
  end do
end subroutine barycentric1d_prepare


function barycentric1d_interpol(x, a_barycentric_1d)
  double precision barycentric1d_interpol
  double precision, intent(in) :: x
  type(type_barycentric_1d), intent(in) :: a_barycentric_1d
  integer i
  double precision a, b, h
  a = 0D0
  b = 0D0
  do i=1, a_barycentric_1d%n
    h = x - a_barycentric_1d%xi(i)
    if (abs(h) .lt. tiny(0D0)*1D2) then
      barycentric1d_interpol = a_barycentric_1d%yi(i)
      return
    else
      if (a_barycentric_1d%d .eq. -1) then
        a = a + 1D0 / (h*h) * a_barycentric_1d%yi(i)
        b = b + 1D0 / (h*h)
      else
        a = a + a_barycentric_1d%wi(i) / h * a_barycentric_1d%yi(i)
        b = b + a_barycentric_1d%wi(i) / h
      end if
    end if
  end do
  barycentric1d_interpol = a / b
end function barycentric1d_interpol

end module barycentric_1d_2d
