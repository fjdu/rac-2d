module grid

use data_struct

implicit none

type :: type_grid_bb
  double precision rmin, rmax, zmin, zmax
end type type_grid_bb

type :: type_refinement_data
  integer :: nlen=0
  double precision, dimension(:,:), allocatable :: xyv
  integer n_idx_incell
  integer, dimension(:), allocatable :: idx_incell
  double precision max_val, min_val, ave_val
  double precision, dimension(2) :: xy_weighted
end type type_refinement_data

type :: type_cell_ptr
  type(type_cell), pointer :: p
end type type_cell_ptr

type :: type_leaves
  integer :: nlen = 0
  type(type_cell_ptr), dimension(:), allocatable :: list
end type type_leaves

double precision :: max_ratio_to_be_uniform = 2.0D0, min_val_considered = 1D-8
double precision :: very_small_len = 1D-7
double precision :: max_cell_size_frac = 0.5D0
double precision :: max_cell_size_0 = 0.5D0
double precision :: aspect_ratio_max = 20.0D0

double precision :: dx_0_frac = 1D-4
double precision :: dy_0_frac = 1D-4
double precision :: dx_0_min_frac = 1D-1
double precision :: dy_0_min_frac = 1D-1
double precision :: dx_ratio = 1.05D0
double precision :: dy_ratio = 1.05D0

type(type_cell), pointer :: root

type(type_leaves) :: cell_leaves

type(type_grid_bb) grid_bb_init


!private :: refinement_data
type(type_refinement_data) refinement_data


contains


! Make grid
subroutine make_grid
  call grid_init
  call load_data_for_refinement
  call grid_refine(root)
  call grid_make_leaves
end subroutine make_grid


subroutine grid_make_leaves
  integer :: idx = 0
  allocate(cell_leaves%list(cell_leaves%nlen))
  call grid_add_leaves(root, idx)
end subroutine grid_make_leaves


recursive subroutine grid_add_leaves(c, idx)
  integer i, idx
  type(type_cell), target :: c
  if (c%nOffspring .EQ. 0) then
    idx = idx + 1
    cell_leaves%list(idx)%p => c
    ! 
    allocate(c%par)
    c%par%n_gas = get_ave_val_analytic(c)
    !
    return
  else
    do i=1, c%nChildren
      call grid_add_leaves(c%children(i), idx)
    end do
  end if
end subroutine grid_add_leaves


function get_ave_val_analytic(c)
  double precision get_ave_val_analytic
  type(type_cell), target :: c
  integer :: i, j, nx=40, ny=50
  double precision dx, dy, x, y
  dx = (c%xmax - c%xmin) / dble(nx-1)
  dy = (c%ymax - c%ymin) / dble(ny-1)
  get_ave_val_analytic = 0D0
  x = c%xmin
  do i=1, nx
    y = c%ymin
    do j=1, ny
      get_ave_val_analytic = get_ave_val_analytic + &
        get_density_analytic(x, y)
      y = y + dy
    end do
    x = x + dx
  end do
  get_ave_val_analytic = get_ave_val_analytic / dble(nx*ny)
end function get_ave_val_analytic


! Initialization
subroutine grid_init
  allocate(root)
  root%xmin = grid_bb_init%rmin
  root%xmax = grid_bb_init%rmax
  root%ymin = grid_bb_init%zmin
  root%ymax = grid_bb_init%zmax
end subroutine grid_init


recursive subroutine grid_refine(c)
  type(type_cell), target :: c
  integer i
  write(*,*) 'XX', c%order, c%nChildren
  if (.not. is_uniform(c)) then
    call sub_divide(c)
  else
    !call grid_decorate(c)
  end if
  !
  do i=1, c%nChildren
    write(*,*) 'YY', c%children(i)%order
    call grid_refine(c%children(i))
  end do
  !
  if (c%nChildren .GT. 0) then
    do i=1, c%nChildren
      c%nOffspring = c%nOffspring + c%children(i)%nOffspring + 1
    end do
  else
    cell_leaves%nlen = cell_leaves%nlen + 1
  end if
end subroutine grid_refine


subroutine sub_divide(c)
  type(type_cell), target :: c
  double precision xmid, ymid, del_x_1, del_x_2, del_y_1, del_y_2
  double precision aspect_1, aspect_2
  integer i
  !
  xmid = refinement_data%xy_weighted(1)
  ymid = refinement_data%xy_weighted(2)
  del_x_1 = xmid - c%xmin
  del_x_2 = c%xmax - xmid
  del_y_1 = ymid - c%ymin
  del_y_2 = c%ymax - ymid
  !
  if ((min(del_x_1, del_x_2) .LE. very_small_len) .AND. &
      (min(del_y_1, del_y_2) .LE. very_small_len)) then
    write(*, '(A, 10ES10.3, 3I5)') 'Find ill-defined cell! (1)', &
      del_x_1, del_x_2, del_y_1, del_y_2, &
      c%xmin, c%xmax, xmid, c%ymin, c%ymax, ymid, refinement_data%n_idx_incell, c%order, c%nChildren
    return
  end if
  !
  c%nChildren = 2
  !call init_children(c, c%nChildren)
  allocate(c%children(c%nChildren))
  do i=1, c%nChildren
    allocate(c%children(i))
  end do
  do i=1, c%nChildren
    c%children(i)%order = c%order + 1
    c%children(i)%parent => c
  write(*,*) 'ZZ', c%children(1)%order
  end do
  !
  aspect_1 = max(max((c%ymax - c%ymin), del_x_1) / min((c%ymax - c%ymin), del_x_1), &
                 max((c%ymax - c%ymin), del_x_2) / min((c%ymax - c%ymin), del_x_2))
  aspect_2 = max(max((c%xmax - c%xmin), del_y_1) / min((c%xmax - c%xmin), del_y_1), &
                 max((c%xmax - c%xmin), del_y_2) / min((c%xmax - c%xmin), del_y_2))
  if ((aspect_1 .LE. aspect_2) .AND. (min(del_x_1, del_x_2) .GT. very_small_len)) then
    c%children(1)%xmin = c%xmin
    c%children(1)%xmax = xmid
    c%children(1)%ymin = c%ymin
    c%children(1)%ymax = c%ymax
    c%children(2)%xmin = xmid
    c%children(2)%xmax = c%xmax
    c%children(2)%ymin = c%ymin
    c%children(2)%ymax = c%ymax
  else
    c%children(1)%xmin = c%xmin
    c%children(1)%xmax = c%xmax
    c%children(1)%ymin = c%ymin
    c%children(1)%ymax = ymid
    c%children(2)%xmin = c%xmin
    c%children(2)%xmax = c%xmax
    c%children(2)%ymin = ymid
    c%children(2)%ymax = c%ymax
  end if
end subroutine sub_divide


! Initialize a cell
recursive subroutine cell_init(c, p, nChildren)
  type(type_cell), target :: c
  type(type_cell), target :: p
  integer nChildren
  integer i
  c%order = p%order + 1
  c%nChildren = nChildren
  c%parent => p
  if (nChildren .gt. 0) then
    allocate(c%children(nChildren))
    do i=1, nChildren
      allocate(c%children(i))
      call cell_init(c%children(i), c, 0)
    end do
  end if
end subroutine cell_init


subroutine init_children(c, nChildren)
  type(type_cell), target :: c
  integer i, nChildren
  if (nChildren .gt. 0) then
    allocate(c%children(nChildren))
    do i=1, nChildren
      allocate(c%children(i))
      !call cell_init(c%children(i), c, 0)
      c%children(i)%order = c%order + 1
      c%children(i)%parent => c
    end do
  end if
end subroutine init_children


subroutine grid_decorate(c)
  ! Esthetic refinement
  type(type_cell), target :: c
  double precision aspect_ratio, delta
  integer i, j, iSub, nsplit
  double precision del_max_x, del_max_y
  double precision del_x, del_y, del1_x, del1_y
  integer ndiv_x, ndiv_y
  ! Max size allowed at this position
  del_max_x = max_cell_size_frac * sqrt(c%xmin*c%xmin + c%ymin*c%ymin) + max_cell_size_0
  del_max_y = del_max_x
  del_x = c%xmax - c%xmin
  del_y = c%ymax - c%ymin
  ndiv_x = ceiling(del_x / del_max_x)
  ndiv_y = ceiling(del_y / del_max_y)
  del1_x = del_x / dble(ndiv_x)
  del1_y = del_y / dble(ndiv_y)
  if (del1_x .ge. del1_y) then
    aspect_ratio = del1_x / del1_y
    ndiv_x = ndiv_x * min(floor(aspect_ratio), ceiling(aspect_ratio/aspect_ratio_max))
  else
    aspect_ratio = del1_y / del1_x
    ndiv_y = ndiv_y * min(floor(aspect_ratio), ceiling(aspect_ratio/aspect_ratio_max))
  end if
  nsplit = ndiv_x * ndiv_y
  if (nsplit .EQ. 1) then
    return
  end if
  del_x = del_x / dble(ndiv_x)
  del_y = del_y / dble(ndiv_y)
  c%nChildren = nsplit
  allocate(c%children(nsplit))
  iSub = 0
  do i=0, ndiv_x-1
    do j=0, ndiv_y-1
      iSub = iSub + 1
      call cell_init(c%children(iSub), c, 0)
      c%children(iSub)%xmin = c%xmin + del_x * dble(i)
      c%children(iSub)%xmax = c%children(iSub)%xmin + del_x
      c%children(iSub)%ymin = c%ymin + del_y * dble(j)
      c%children(iSub)%ymax = c%children(iSub)%ymin + del_y
      ! Get rid of the rounding errors.
      if (i .EQ. 0)          c%children(iSub)%xmin = c%xmin
      if (i .EQ. (ndiv_x-1)) c%children(iSub)%xmax = c%xmax
      if (j .EQ. 0)          c%children(iSub)%ymin = c%ymin
      if (j .EQ. (ndiv_y-1)) c%children(iSub)%ymax = c%ymax
    end do
  end do
end subroutine grid_decorate


function is_uniform(c)
  logical is_uniform
  type(type_cell), target :: c
  !if (c%order .lt. 9955) then
    is_uniform = is_uniform_based_on_data(c)
  !else
  !  is_uniform = is_uniform_analytic(c)
  !end if
end function is_uniform


function is_uniform_analytic(c)
  logical is_uniform_analytic
  type(type_cell), target :: c
  integer, parameter :: nx=4, ny=5
  double precision, dimension(nx, ny) :: vals
  integer i, j
  double precision dx, dy, x, y
  dx = (c%xmax - c%xmin) / dble(nx-1)
  dy = (c%ymax - c%ymin) / dble(ny-1)
  associate(d    => refinement_data, &
            n    => refinement_data%n_idx_incell, &
            idx  => refinement_data%idx_incell, &
            maxv => refinement_data%max_val, &
            minv => refinement_data%min_val, &
            avev => refinement_data%ave_val, &
            xyw  => refinement_data%xy_weighted)
    n = -1
    idx = -1
    xyw = 0D0
    avev = 0D0
    x = c%xmin
    do i=1, nx
      y = c%ymin
      do j=1, ny
        vals(i, j) = get_density_analytic(x, y)
        xyw(1) = xyw(1) + vals(i, j) * x
        xyw(2) = xyw(2) + vals(i, j) * y
        avev = avev + vals(i, j)
        y = y + dy
      end do
      x = x + dx
    end do
    maxv = maxval(vals)
    minv = minval(vals)
    xyw  = xyw / avev
    avev = avev / dble(nx*ny)
    if (maxv .gt. min_val_considered) then
      if (maxv / (minv + tiny(0D0)) .gt. max_ratio_to_be_uniform) then
        is_uniform_analytic = .false.
      else
        is_uniform_analytic = .true.
      end if
    else
      is_uniform_analytic = .true.
    end if
    !write(*,*) xyw, avev, dx, dy, is_uniform_analytic
    write(*,*) c%order, c%xmin, c%xmax, c%ymin, c%ymax
  end associate
end function is_uniform_analytic


function is_uniform_based_on_data(c)
  logical is_uniform_based_on_data
  type(type_cell), target :: c
  integer i, n_in
  associate(d    => refinement_data, &
            n    => refinement_data%n_idx_incell, &
            idx  => refinement_data%idx_incell, &
            maxv => refinement_data%max_val, &
            minv => refinement_data%min_val, &
            avev => refinement_data%ave_val, &
            xyw  => refinement_data%xy_weighted)
    n = 0
    maxv = 0D0
    minv = huge(0D0)
    avev = 0D0
    xyw  = 0D0
    do i=1, d%nlen
      if (is_inside_rect(d%xyv(1:2, i), c%xmin, c%xmax, c%ymin, c%ymax)) then
        n = n + 1
        idx(n) = i
        maxv = max(maxv, d%xyv(3, i))
        minv = min(minv, d%xyv(3, i))
        avev = avev + d%xyv(3, i)
        xyw  = xyw + d%xyv(1:2, i) * d%xyv(3, i)
      end if
    end do
    xyw  = xyw / avev
    avev = avev / dble(n)
    if (n .gt. 0) then
      if (maxv .gt. min_val_considered) then
        if (maxv / (minv + tiny(0D0)) .gt. (max_ratio_to_be_uniform + max(0D0,-log10(maxv)/3D0))) then
          is_uniform_based_on_data = .false.
        else
          is_uniform_based_on_data = .true.
        end if
      else
        is_uniform_based_on_data = .true.
      end if
    else
      is_uniform_based_on_data = .true.
    end if
  end associate
end function is_uniform_based_on_data


function is_inside_rect(v, xmin, xmax, ymin, ymax)
  logical is_inside_rect
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision, dimension(2) :: v
  if ( (v(1) .ge. xmin) .and. (v(1) .le. xmax) .and. &
       (v(2) .ge. ymin) .and. (v(2) .le. ymax)) then
    is_inside_rect = .true.
  else
    is_inside_rect = .false.
  end if
end function is_inside_rect


subroutine load_data_for_refinement
  call load_data_for_refinement_analytic
end subroutine load_data_for_refinement


subroutine load_data_for_refinement_analytic
  use quick_sort
  integer nx, ny
  double precision, dimension(:,:), allocatable :: x, y, z
  integer i, j
  double precision xmin, xmax, ymin, ymax
  double precision dx, dx_0
  double precision dy, dy_0
  !
  xmin = root%xmin
  xmax = root%xmax
  ymin = root%ymin
  ymax = root%ymax
  !
  dx_0 = max(xmin * dx_0_min_frac, (xmax - xmin) * dx_0_frac)
  dy_0 = (ymax - ymin) * dy_0_frac
  !
  nx = ceiling(log( &
         (xmax-xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (ymax-ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  write(*,*) 'nx, ny ', nx, ny
  !
  allocate(x(nx, ny), y(nx, ny), z(nx, ny))
  !
  x(1, :) = xmin
  dx = dx_0
  do i=1, nx
    if (i .gt. 1) then
      x(i, :) = x(i-1, :) + dx
      dx = dx * dx_ratio
    end if
    y(i, 1) = 0D0
    dy = dy_0 * min(1D0, dy_0_min_frac + x(i, 1) / xmax)
    do j=1, ny
      if (j .gt. 1) then
        y(i, j) = y(i, j-1) + dy
        dy = dy * dy_ratio
      end if
      z(i, j) = get_density_analytic(x(i,j), y(i,j))
    end do
  end do
  associate(d => refinement_data)
    d%nlen = nx * ny
    allocate(d%xyv(3, d%nlen), d%idx_incell(d%nlen))
    d%xyv(1, :) = reshape(x, (/d%nlen/))
    d%xyv(2, :) = reshape(y, (/d%nlen/))
    d%xyv(3, :) = reshape(z, (/d%nlen/))
    ! Sort data
    !call quick_sort_array(d%xyv, 3, d%nlen, 2, (/1, 2/))
  end associate
  deallocate(x, y, z)
end subroutine load_data_for_refinement_analytic


function get_density_analytic(x, y)
  double precision get_density_analytic
  double precision, intent(in) :: x, y
  get_density_analytic = density_analytic_Hayashi(x, y)
end function get_density_analytic
! Sub-divide method

! Merge method

! Neighborhood relationship


function density_analytic_Hayashi(r, z)
  ! Hayashi 1981, equations 2.7, 2.8, 2.9
  use phy_const
  double precision, intent(in) :: r, z ! in AU
  double precision density_analytic_Hayashi ! in number cm-3
  double precision, parameter :: rho0 = 1.4D-9 ! g cm-3
  double precision, parameter :: n0 = rho0 / (1.4D0 * phy_mProton_CGS)
  double precision z0
  z0 = 0.0472D0 * r**1.25D0
  density_analytic_Hayashi = n0 * r**(-2.75D0) * exp(-(z/z0)**2)
end function density_analytic_Hayashi


end module grid
