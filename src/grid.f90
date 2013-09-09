module grid

use data_struct
use phy_const
use trivials
use spline_1d_2d
use barycentric_1d_2d

implicit none

type :: type_refinement_data
  integer :: nlen=0, ncol=0, nx=0, ny=0
  double precision, dimension(:,:), allocatable :: xyv
  integer n_idx_incell
  integer, dimension(:), allocatable :: idx_incell
  double precision max_val, min_val
  double precision, dimension(:), allocatable :: ave_val
  double precision, dimension(2) :: xy_weighted
end type type_refinement_data

type :: type_cell_ptr
  type(type_cell), pointer :: p
end type type_cell_ptr

type :: type_leaves
  integer :: nlen = 0
  type(type_cell_ptr), dimension(:), allocatable :: list
end type type_leaves

type :: type_grid_config
  double precision rmin, rmax, zmin, zmax
  logical :: use_data_file_input = .false.
  logical :: columnwise = .true.
  integer :: ncol
  character(len=128) :: data_dir = './'
  character(len=128) :: data_filename = ''
  character(len=16) :: analytical_to_use = 'Andrew'
  character(len=16) :: interpolation_method = 'barycentric'
  double precision :: max_ratio_to_be_uniform = 5.0D0
  double precision :: min_val_considered = 5D1
  double precision :: very_small_len = 1D-3
  double precision :: small_len_frac = 1D-2
end type type_grid_config

type :: type_leaves_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: idx
end type type_leaves_list


type(type_leaves) :: cell_leaves

type(type_leaves_list) :: surf_cells, bott_cells

type(type_grid_config) :: grid_config

private :: refinement_data
type(type_refinement_data) refinement_data

type(type_cell), pointer :: root

type(type_spline_2D), allocatable :: n_spline2d, T_spline2d
type(type_barycentric_2d), allocatable :: n_bary2d, T_bary2d

double precision :: MeanMolWeight = 1.4D0
double precision :: RADMC_gas2dust_mass_ratio = 1D2 ! Todo

double precision, parameter, private :: const_uniform_a = 0.1D0, const_uniform_b = 9D0

namelist /grid_configure/ &
  grid_config


contains


subroutine make_grid
  integer i
  if (grid_config%use_data_file_input) then
    call load_data_for_refinement
  end if
  !
  allocate(root)
  if (grid_config%columnwise) then
    call grid_init_columnwise(root)
    do i=1, root%nChildren
      call grid_refine(root%children(i)%p)
      root%nOffspring = root%nOffspring + root%children(i)%p%nOffspring + 1
      root%nleaves = root%nleaves + root%children(i)%p%nleaves
    end do
  else
    call grid_init(root)
    call grid_refine(root)
  end if
  cell_leaves%nlen = root%nleaves
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
end subroutine make_grid


subroutine grid_make_surf_bott
  integer i
  integer, parameter :: nSurfBott_max = 1024
  integer, dimension(nSurfBott_max) :: idxSurf, idxBott
  double precision, parameter :: fra_tot_thre = 0.99D0
  do i=1, cell_leaves%nlen
    !if (cell_leaves%list(i)%p%above%n .eq. 0) then
    if (cell_leaves%list(i)%p%above%fra_tot .lt. fra_tot_thre) then
      surf_cells%nlen = surf_cells%nlen + 1
      idxSurf(surf_cells%nlen) = i
    end if
    !if (cell_leaves%list(i)%p%below%n .eq. 0) then
    if (cell_leaves%list(i)%p%below%fra_tot .lt. fra_tot_thre) then
      bott_cells%nlen = bott_cells%nlen + 1
      idxBott(bott_cells%nlen) = i
    end if
  end do
  if (surf_cells%nlen .gt. 0) then
    if (allocated(surf_cells%idx)) then
      deallocate(surf_cells%idx)
    end if
    allocate(surf_cells%idx(surf_cells%nlen))
    surf_cells%idx = idxSurf(1:surf_cells%nlen)
  end if
  if (bott_cells%nlen .gt. 0) then
    if (allocated(bott_cells%idx)) then
      deallocate(bott_cells%idx)
    end if
    allocate(bott_cells%idx(bott_cells%nlen))
    bott_cells%idx = idxBott(1:bott_cells%nlen)
  end if
end subroutine grid_make_surf_bott


subroutine grid_make_leaves(c)
  type(type_cell), target, intent(in) :: c
  integer :: idx = 0
  if (allocated(cell_leaves%list)) then
    deallocate(cell_leaves%list)
  end if
  allocate(cell_leaves%list(cell_leaves%nlen))
  call grid_add_leaves(c, idx)
end subroutine grid_make_leaves


recursive subroutine get_number_of_leaves(c)
  integer i
  type(type_cell), target :: c
  if (c%using) then
    c%nleaves = 1
    c%nOffspring = 0
    return
  else
    c%nleaves = 0
    c%nOffspring = 0
    do i=1, c%nChildren
      call get_number_of_leaves(c%children(i)%p)
      c%nleaves = c%nleaves + c%children(i)%p%nleaves
      c%nOffspring = c%nOffspring + c%children(i)%p%nOffspring + 1
    end do
  end if
end subroutine get_number_of_leaves


recursive subroutine grid_add_leaves(c, idx)
  integer i, idx
  type(type_cell), target :: c
  if (c%using) then
    idx = idx + 1
    cell_leaves%list(idx)%p => c
    return
  else
    do i=1, c%nChildren
      call grid_add_leaves(c%children(i)%p, idx)
    end do
  end if
end subroutine grid_add_leaves


subroutine grid_make_neighbors
  integer i
  do i=1, cell_leaves%nlen
    call make_neighbors(cell_leaves%list(i)%p)
  end do
end subroutine grid_make_neighbors


function get_ave_val_analytic(xmin, xmax, ymin, ymax)
  double precision get_ave_val_analytic
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  integer :: i, j, nx, ny
  double precision dx, dy, dx0, dy0, x, y, area, dely
  double precision, parameter :: dx0_frac = 1D-3, dy0_frac = 1D-4, dx_ratio = 1.5D0, dy_ratio=1.5D0
  dx0 = max((xmax - xmin) * dx0_frac, grid_config%very_small_len*0.01D0)
  dy0 = max((ymax - ymin) * dy0_frac, grid_config%very_small_len*0.01D0)
  nx = ceiling(log( &
         (xmax-xmin)/dx0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (ymax-ymin)/dy0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  get_ave_val_analytic = 0D0
  area = 0D0
  x = xmin
  dx = dx0
  do i=1, nx
    if (x .gt. xmax) then
      exit
    end if
    !
    get_ave_val_analytic = get_ave_val_analytic + &
      get_int_val_along_y(x+0.5D0*dx, ymin, ymax, dely) * dx
    area = area + dx * dely
    x = x + dx
    dx = dx * dx_ratio
  end do
  if (area .eq. 0D0) then
    write(*,*) 'Area  = 0!', xmin, xmax, dx0, dely
  end if
  get_ave_val_analytic = get_ave_val_analytic / area
end function get_ave_val_analytic


function get_ave_val_based_on_data(xmin, xmax, ymin, ymax)
  double precision, dimension(:), allocatable :: get_ave_val_based_on_data
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision, dimension(4) :: vals
  allocate(get_ave_val_based_on_data(refinement_data%nlen-2))
  vals(1) =  max(0D0, get_RADMC_n(xmin, ymin))
  vals(2) =  max(0D0, get_RADMC_n(xmin, ymax))
  vals(3) =  max(0D0, get_RADMC_n(xmax, ymin))
  vals(4) =  max(0D0, get_RADMC_n(xmax, ymax))
  get_ave_val_based_on_data(1) = sum(vals) * 0.25D0
  vals(1) =  max(0D0, get_RADMC_T(xmin, ymin))
  vals(2) =  max(0D0, get_RADMC_T(xmin, ymax))
  vals(3) =  max(0D0, get_RADMC_T(xmax, ymin))
  vals(4) =  max(0D0, get_RADMC_T(xmax, ymax))
  get_ave_val_based_on_data(2) = sum(vals) * 0.25D0
end function get_ave_val_based_on_data


subroutine grid_init(c)
  type(type_cell), target :: c
  c%xmin = grid_config%rmin
  c%xmax = grid_config%rmax
  c%ymin = grid_config%zmin
  c%ymax = grid_config%zmax
end subroutine grid_init


subroutine grid_init_columnwise(c)
  type(type_cell), target :: c
  double precision dx0, del_ratio
  double precision x, dx, del, tmp
  integer i
  c%xmin = grid_config%rmin
  c%xmax = grid_config%rmax
  c%ymin = grid_config%zmin
  c%ymax = grid_config%zmax
  c%nChildren = grid_config%ncol
  !dx0 = max(1D-2, (c%xmax - c%xmin) * 2D-4)
  tmp = log(c%xmax/c%xmin) / dble(c%nChildren)
  if (tmp .le. 1D-4) then
    dx0 = tmp * c%xmin
  else
    dx0 = (exp(tmp) - 1D0) * c%xmin
  end if
  !del_ratio = get_ratio_of_interval_log(c%xmin, c%xmax, dx0, c%nChildren)
  del_ratio = exp(tmp)
  call init_children(c, c%nChildren)
  dx = dx0
  x = c%xmin
  do i=1, c%nChildren
    c%children(i)%p%xmin = x
    c%children(i)%p%xmax = x + dx
    c%children(i)%p%ymin = c%ymin
    c%children(i)%p%ymax = get_ymax_here(x + 0.5D0 * dx, c%ymin, c%ymax)
    x = x + dx
    dx = dx * del_ratio
  end do
  c%children(c%nChildren)%p%xmax = c%xmax
end subroutine grid_init_columnwise


function get_ymax_here(x, y0, y1)
  double precision get_ymax_here
  double precision, intent(in) :: x, y0, y1
  integer :: i, n = 100
  double precision y, dy, del_ratio
  dy = (y1 - y0) * 1D-4
  del_ratio = get_ratio_of_interval_log(y0, y1, dy, n)
  dy = dy * del_ratio**(n-1)
  y = y1
  do i=1, n
    if (get_RADMC_n(x, y) .ge. grid_config%min_val_considered) then
      get_ymax_here = y
      return
    end if
    y = y - dy
    dy = dy / del_ratio
  end do
end function get_ymax_here



recursive subroutine grid_refine(c)
  type(type_cell), target :: c
  integer i
  if (.not. is_uniform(c)) then
    call sub_divide(c)
  else
    !call grid_decorate(c)
  end if
  !
  do i=1, c%nChildren
    call grid_refine(c%children(i)%p)
  end do
  !
  if (c%nChildren .gt. 0) then
    do i=1, c%nChildren
      c%nOffspring = c%nOffspring + c%children(i)%p%nOffspring + 1
      c%nleaves = c%nleaves + c%children(i)%p%nleaves
    end do
  else
    call set_cell_par_preliminary(c)
    if (c%using) then
      c%nleaves = 1
    end if
  end if
end subroutine grid_refine


subroutine set_cell_par_preliminary(c)
  type(type_cell), target :: c
  if (grid_config%use_data_file_input) then
    c%val = get_ave_val_based_on_data(c%xmin, c%xmax, c%ymin, c%ymax)
  else
    c%val(1) = get_ave_val_analytic(c%xmin, c%xmax, c%ymin, c%ymax)
  end if
  c%using = c%val(1) .gt. grid_config%min_val_considered
end subroutine set_cell_par_preliminary


subroutine make_neighbors(c)
  type(type_cell), target :: c
  integer i
  integer n, pos
  integer, parameter :: nnei_max = 32
  integer, dimension(nnei_max) :: idx_inner, idx_outer, idx_below, idx_above
  double precision, dimension(nnei_max) :: fra_inner, fra_outer, fra_below, fra_above
  double precision frac
  if (.not. associated(c%inner)) then
    allocate(c%inner, c%outer, c%below, c%above, c%around)
  else
    if (allocated(c%inner%idx)) then
      deallocate(c%inner%idx, c%inner%fra)
    end if
    if (allocated(c%outer%idx)) then
      deallocate(c%outer%idx, c%outer%fra)
    end if
    if (allocated(c%above%idx)) then
      deallocate(c%above%idx, c%above%fra)
    end if
    if (allocated(c%below%idx)) then
      deallocate(c%below%idx, c%below%fra)
    end if
    if (allocated(c%around%idx)) then
      deallocate(c%around%idx, c%around%fra)
    end if
  end if
  do i=1, cell_leaves%nlen
    if (is_neighbor(c, cell_leaves%list(i)%p, 0.1D0*grid_config%very_small_len, pos, frac)) then
      select case(pos)
        case(1)
          c%inner%n = c%inner%n + 1
          idx_inner(c%inner%n) = i
          fra_inner(c%inner%n) = frac
        case(2)
          c%outer%n = c%outer%n + 1
          idx_outer(c%outer%n) = i
          fra_outer(c%outer%n) = frac
        case(3)
          c%below%n = c%below%n + 1
          idx_below(c%below%n) = i
          fra_below(c%below%n) = frac
        case(4)
          c%above%n = c%above%n + 1
          idx_above(c%above%n) = i
          fra_above(c%above%n) = frac
      end select
    end if
  end do
  c%around%n = c%inner%n + c%outer%n + c%below%n + c%above%n
  if (c%around%n .gt. 0) then
    allocate(c%around%idx(c%around%n), c%around%fra(c%around%n))
  end if
  i = 0
  if (c%above%n .gt. 0) then
    allocate(c%above%idx(c%above%n), c%above%fra(c%above%n))
    c%above%idx = idx_above(1:c%above%n)
    c%above%fra = fra_above(1:c%above%n)
    c%above%fra_tot = sum(fra_above(1:c%above%n))
    c%around%idx(i+1:i+c%above%n) = c%above%idx
    c%around%fra(i+1:i+c%above%n) = c%above%fra
    i = i + c%above%n
  end if
  if (c%inner%n .gt. 0) then
    allocate(c%inner%idx(c%inner%n), c%inner%fra(c%inner%n))
    c%inner%idx = idx_inner(1:c%inner%n)
    c%inner%fra = fra_inner(1:c%inner%n)
    c%inner%fra_tot = sum(fra_inner(1:c%inner%n))
    c%around%idx(i+1:i+c%inner%n) = c%inner%idx
    c%around%fra(i+1:i+c%inner%n) = c%inner%fra
    i = i + c%inner%n
  end if
  if (c%outer%n .gt. 0) then
    allocate(c%outer%idx(c%outer%n), c%outer%fra(c%outer%n))
    c%outer%idx = idx_outer(1:c%outer%n)
    c%outer%fra = fra_outer(1:c%outer%n)
    c%outer%fra_tot = sum(fra_outer(1:c%outer%n))
    c%around%idx(i+1:i+c%outer%n) = c%outer%idx
    c%around%fra(i+1:i+c%outer%n) = c%outer%fra
    i = i + c%outer%n
  end if
  if (c%below%n .gt. 0) then
    allocate(c%below%idx(c%below%n), c%below%fra(c%below%n))
    c%below%idx = idx_below(1:c%below%n)
    c%below%fra = fra_below(1:c%below%n)
    c%below%fra_tot = sum(fra_below(1:c%below%n))
    c%around%idx(i+1:i+c%below%n) = c%below%idx
    c%around%fra(i+1:i+c%below%n) = c%below%fra
    i = i + c%below%n
  end if
end subroutine make_neighbors


function is_neighbor(c1, c2, tol, pos, frac)
  logical is_neighbor
  type(type_cell), target, intent(in) :: c1, c2
  double precision, intent(in) :: tol
  integer, intent(out) :: pos
  double precision, intent(out) :: frac
  is_neighbor = .false.
  if ((abs(c2%xmax - c1%xmin) .le. tol) .and. &
      Two_lines_has_common_section(c1%ymin, c1%ymax, c2%ymin, c2%ymax, frac)) then
    is_neighbor = .true.
    pos = 1
    frac = frac / (c1%ymax - c1%ymin)
    return
  end if
  if ((abs(c2%xmin - c1%xmax) .le. tol) .and. &
      Two_lines_has_common_section(c1%ymin, c1%ymax, c2%ymin, c2%ymax, frac)) then
    is_neighbor = .true.
    pos = 2
    frac = frac / (c1%ymax - c1%ymin)
    return
  end if
  if ((abs(c2%ymax - c1%ymin) .le. tol) .and. &
      Two_lines_has_common_section(c1%xmin, c1%xmax, c2%xmin, c2%xmax, frac)) then
    is_neighbor = .true.
    pos = 3
    frac = frac / (c1%xmax - c1%xmin)
    return
  end if
  if ((abs(c2%ymin - c1%ymax) .le. tol) .and. &
      Two_lines_has_common_section(c1%xmin, c1%xmax, c2%xmin, c2%xmax, frac)) then
    is_neighbor = .true.
    pos = 4
    frac = frac / (c1%xmax - c1%xmin)
    return
  end if
end function is_neighbor


function Two_lines_has_common_section(xmin, xmax, ymin, ymax, len_sec)
  logical Two_lines_has_common_section
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision, optional, intent(out) :: len_sec
  if ((xmax .le. ymin) .or. (ymax .le. xmin)) then
    Two_lines_has_common_section = .false.
  else
    Two_lines_has_common_section = .true.
    if (present(len_sec)) then
      len_sec = min(xmax, ymax) - max(xmin, ymin)
    end if
  end if
end function Two_lines_has_common_section


function is_uniform(c)
  logical is_uniform
  type(type_cell), target :: c
  if (grid_config%columnwise) then
    is_uniform = test_uniformity_columnwise(c%xmin, c%xmax, c%ymin, c%ymax)
  else
    is_uniform = test_uniformity(c%xmin, c%xmax, c%ymin, c%ymax)
  end if
end function is_uniform


subroutine sub_divide(c)
  type(type_cell), target :: c
  if (grid_config%columnwise) then
    call find_mid_columnwise(c)
    call sub_divide_columnwise(c)
  else
    call find_mid_quadtree(c)
    call sub_divide_8cases(c)
  end if
end subroutine sub_divide


subroutine sub_divide_8cases(c)
  type(type_cell), target :: c
  double precision xmid, ymid, del_x_1, del_x_2, del_y_1, del_y_2
  double precision min_del_x, min_del_y
  integer i, icase
  logical unif_L, unif_R, unif_B, unif_T
  double precision small_len_this
  !
  xmid = refinement_data%xy_weighted(1)
  ymid = refinement_data%xy_weighted(2)
  del_x_1 = xmid - c%xmin
  del_x_2 = c%xmax - xmid
  del_y_1 = ymid - c%ymin
  del_y_2 = c%ymax - ymid
  !
  small_len_this = max(sqrt(xmid*xmid + ymid*ymid) * grid_config%small_len_frac, &
                      grid_config%very_small_len)
  !
  min_del_x = min(del_x_1, del_x_2)
  min_del_y = min(del_y_1, del_y_2)
  if ((min_del_x .LE. small_len_this) .and. &
      (min_del_y .LE. small_len_this)) then
    icase = 1
  else if ((min_del_x .gt. small_len_this) .and. &
           (min_del_y .le. small_len_this)) then
    icase = 2
  else if ((min_del_x .le. small_len_this) .and. &
           (min_del_y .gt. small_len_this)) then
    icase = 3
  else
    unif_L = test_uniformity(c%xmin, xmid, c%ymin, c%ymax)
    unif_R = test_uniformity(xmid, c%xmax, c%ymin, c%ymax)
    unif_B = test_uniformity(c%xmin, c%xmax, c%ymin, ymid)
    unif_T = test_uniformity(c%xmin, c%xmax, ymid, c%ymax)
    if (unif_L .and. unif_R) then
      icase = 2
    else if (unif_B .and. unif_T) then
      icase = 3
    else if ((.not. unif_L) .and. unif_R) then
      icase = 4
    else if (unif_L .and. (.not. unif_R)) then
      icase = 5
    else if ((.not. unif_B) .and. unif_T) then
      icase = 6
    else if (unif_B .and. (.not. unif_T)) then
      icase = 7
    else if ((.not. unif_B) .and. (.not. unif_T)) then
      icase = 8
    end if
  end if
  select case(icase)
    case(1)
      return
    case(2)
      c%nChildren = 2
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = xmid
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = c%ymax
      c%children(2)%p%xmin = xmid
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = c%ymin
      c%children(2)%p%ymax = c%ymax
    case(3)
      c%nChildren = 2
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = c%xmax
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = ymid
      c%children(2)%p%xmin = c%xmin
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = ymid
      c%children(2)%p%ymax = c%ymax
    case(4)
      c%nChildren = 3
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = xmid
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = ymid
      c%children(2)%p%xmin = c%xmin
      c%children(2)%p%xmax = xmid
      c%children(2)%p%ymin = ymid
      c%children(2)%p%ymax = c%ymax
      c%children(3)%p%xmin = xmid
      c%children(3)%p%xmax = c%xmax
      c%children(3)%p%ymin = c%ymin
      c%children(3)%p%ymax = c%ymax
    case(5)
      c%nChildren = 3
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = xmid
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = c%ymax
      c%children(2)%p%xmin = xmid
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = c%ymin
      c%children(2)%p%ymax = ymid
      c%children(3)%p%xmin = xmid
      c%children(3)%p%xmax = c%xmax
      c%children(3)%p%ymin = ymid
      c%children(3)%p%ymax = c%ymax
    case(6)
      c%nChildren = 3
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = xmid
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = ymid
      c%children(2)%p%xmin = xmid
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = c%ymin
      c%children(2)%p%ymax = ymid
      c%children(3)%p%xmin = c%xmin
      c%children(3)%p%xmax = c%xmax
      c%children(3)%p%ymin = ymid
      c%children(3)%p%ymax = c%ymax
    case(7)
      c%nChildren = 3
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = c%xmax
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = ymid
      c%children(2)%p%xmin = xmid
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = ymid
      c%children(2)%p%ymax = c%ymax
      c%children(3)%p%xmin = c%xmin
      c%children(3)%p%xmax = xmid
      c%children(3)%p%ymin = ymid
      c%children(3)%p%ymax = c%ymax
    case(8)
      c%nChildren = 4
      call init_children(c, c%nChildren)
      c%children(1)%p%xmin = c%xmin
      c%children(1)%p%xmax = xmid
      c%children(1)%p%ymin = c%ymin
      c%children(1)%p%ymax = ymid
      c%children(2)%p%xmin = xmid
      c%children(2)%p%xmax = c%xmax
      c%children(2)%p%ymin = c%ymin
      c%children(2)%p%ymax = ymid
      c%children(3)%p%xmin = xmid
      c%children(3)%p%xmax = c%xmax
      c%children(3)%p%ymin = ymid
      c%children(3)%p%ymax = c%ymax
      c%children(4)%p%xmin = c%xmin
      c%children(4)%p%xmax = xmid
      c%children(4)%p%ymin = ymid
      c%children(4)%p%ymax = c%ymax
  end select
end subroutine sub_divide_8cases


subroutine sub_divide_columnwise(c)
  type(type_cell), target :: c
  double precision xmid, ymid, del_y_1, del_y_2
  double precision min_del_y
  double precision small_len_this
  !
  xmid = refinement_data%xy_weighted(1)
  ymid = refinement_data%xy_weighted(2)
  del_y_1 = ymid - c%ymin
  del_y_2 = c%ymax - ymid
  !
  small_len_this = max(sqrt(xmid*xmid + ymid*ymid) * grid_config%small_len_frac, &
                      grid_config%very_small_len)
  !
  min_del_y = min(del_y_1, del_y_2)
  if (min_del_y .LE. small_len_this) then
    return
  end if
  c%nChildren = 2
  call init_children(c, c%nChildren)
  c%children(1)%p%xmin = c%xmin
  c%children(1)%p%xmax = c%xmax
  c%children(1)%p%ymin = c%ymin
  c%children(1)%p%ymax = ymid
  c%children(2)%p%xmin = c%xmin
  c%children(2)%p%xmax = c%xmax
  c%children(2)%p%ymin = ymid
  c%children(2)%p%ymax = c%ymax
end subroutine sub_divide_columnwise


subroutine find_mid_quadtree(c)
  type(type_cell), target :: c
  associate( &
            xmid => refinement_data%xy_weighted(1), &
            ymid => refinement_data%xy_weighted(2))
    xmid = 0.5D0 * (c%xmin + c%xmax)
    ymid = 0.5D0 * (c%ymin + c%ymax)
  end associate
end subroutine find_mid_quadtree


subroutine find_mid_columnwise(c)
  use trivials
  type(type_cell), target :: c
  associate( &
            xmid => refinement_data%xy_weighted(1), &
            ymid => refinement_data%xy_weighted(2))
    xmid = dblNaN() ! Undefined
    ymid = 0.5D0 * (c%ymin + c%ymax)
  end associate
end subroutine find_mid_columnwise


function test_uniformity(xmin, xmax, ymin, ymax)
  logical test_uniformity
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  if (grid_config%use_data_file_input) then
    test_uniformity = test_uniformity_based_on_data(xmin, xmax, ymin, ymax)
  else
    test_uniformity = test_uniformity_simple_analytic(xmin, xmax, ymin, ymax)
  end if
end function test_uniformity


function test_uniformity_columnwise(xmin, xmax, ymin, ymax)
  logical test_uniformity_columnwise
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  if (grid_config%use_data_file_input) then
    test_uniformity_columnwise = test_uniformity_based_on_data_columnwise(xmin, xmax, ymin, ymax)
  else
    test_uniformity_columnwise = test_uniformity_simple_analytic_columnwise(xmin, xmax, ymin, ymax)
  end if
end function test_uniformity_columnwise


function test_uniformity_simple_analytic(xmin, xmax, ymin, ymax)
  logical test_uniformity_simple_analytic
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  integer, parameter :: n = 4
  double precision, dimension(n) :: vals
  double precision minv, maxv, max_ratio_to_be_uniform_here!, Av_vert_0
  integer i
  double precision, parameter :: const_small_num = 1D-100
  vals(1) = get_density_analytic(xmin, ymin)
  vals(2) = get_density_analytic(xmax, ymin)
  vals(3) = get_density_analytic(xmax, ymax)
  vals(4) = get_density_analytic(xmin, ymax)
  maxv = maxval(vals)
  minv = minval(vals)
  max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
    const_uniform_a * (log10(maxv) - const_uniform_b)**2
  if ((maxv .le. grid_config%min_val_considered)) then
    test_uniformity_simple_analytic = .true.
  else if (maxv / (minv + const_small_num) .le. max_ratio_to_be_uniform_here) then
    test_uniformity_simple_analytic = .true.
  else
    test_uniformity_simple_analytic = .false.
  end if
end function test_uniformity_simple_analytic


function test_uniformity_simple_analytic_columnwise(xmin, xmax, ymin, ymax)
  logical test_uniformity_simple_analytic_columnwise
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  integer, parameter :: n = 3
  double precision, dimension(n) :: vals
  double precision minv, maxv, max_ratio_to_be_uniform_here!, Av_vert_0
  integer i
  double precision, parameter :: const_small_num = 1D-100
  vals(1) = get_density_analytic(0.5D0*(xmin+xmax), ymin)
  vals(2) = get_density_analytic(0.5D0*(xmin+xmax), 0.5D0*(ymin+ymax))
  vals(3) = get_density_analytic(0.5D0*(xmin+xmax), ymax)
  maxv = maxval(vals)
  minv = minval(vals)
  max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
    const_uniform_a * (log10(maxv) - const_uniform_b)**2
  if ((maxv .le. grid_config%min_val_considered)) then
    test_uniformity_simple_analytic_columnwise = .true.
  else if (maxv / (minv + const_small_num) .le. max_ratio_to_be_uniform_here) then
    test_uniformity_simple_analytic_columnwise = .true.
  else
    test_uniformity_simple_analytic_columnwise = .false.
  end if
end function test_uniformity_simple_analytic_columnwise


function test_uniformity_based_on_data(xmin, xmax, ymin, ymax)
  logical test_uniformity_based_on_data
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  integer i, n_in
  double precision max_ratio_to_be_uniform_here
  double precision, parameter :: const_small_num = 1D-100
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
    !xyw  = 0D0
    do i=1, d%nlen
      !if (is_inside_rect(d%xyv(1:2, i), xmin, xmax, ymin, ymax)) then
      if (is_close_to_rect(d%xyv(1:2, i), xmin, xmax, ymin, ymax)) then
        n = n + 1
        idx(n) = i
        maxv = max(maxv, d%xyv(3, i))
        minv = min(minv, d%xyv(3, i))
        avev = avev + d%xyv(3, i)
        !xyw  = xyw + d%xyv(1:2, i) * d%xyv(3, i)
      end if
    end do
    !xyw  = xyw / avev
    avev = avev / dble(n)
    if (n .gt. 0) then
      max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
        const_uniform_a * (log10(maxv) - const_uniform_b)**2
      if ((maxv .le. grid_config%min_val_considered)) then
        test_uniformity_based_on_data = .true.
      else if (maxv / (minv + const_small_num) .le. max_ratio_to_be_uniform_here) then
        test_uniformity_based_on_data = .true.
      else
        test_uniformity_based_on_data = .false.
      end if
    else
      test_uniformity_based_on_data = .true.
    end if
  end associate
end function test_uniformity_based_on_data


function test_uniformity_based_on_data_columnwise(xmin, xmax, ymin, ymax)
  logical test_uniformity_based_on_data_columnwise
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision x0
  double precision max_ratio_to_be_uniform_here
  integer, parameter :: n = 10
  integer i
  double precision dely, y
  double precision, dimension(n) :: vals
  double precision, parameter :: const_small_num = 1D-100
  associate(d    => refinement_data, &
            maxv => refinement_data%max_val, &
            minv => refinement_data%min_val)
    x0 = 0.5D0 * (xmin + xmax)
    dely = (ymax - ymin) / dble(n-1)
    y = ymin
    do i=1, n
      vals(i) = get_RADMC_n(x0, y)
      !if (vals(i) .lt. 0D0) then
      !  write(*,*) x0, y, vals(i)
      !end if
      y = y + dely
    end do
    maxv = maxval(vals)
    minv = minval(vals)
    if (maxv .le. grid_config%min_val_considered) then
      test_uniformity_based_on_data_columnwise = .true.
    else
      max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
        const_uniform_a * (log10(maxv) - const_uniform_b)**2
      if (maxv / (abs(minv) + const_small_num) .le. max_ratio_to_be_uniform_here) then
        test_uniformity_based_on_data_columnwise = .true.
      else
        test_uniformity_based_on_data_columnwise = .false.
      end if
    end if
    !write(*,*) xmin, xmax, ymin, ymax, maxv, minv, test_uniformity_based_on_data_columnwise
  end associate
end function test_uniformity_based_on_data_columnwise


function get_int_val_along_y(x, lmin, lmax, del_span)
  double precision get_int_val_along_y
  double precision, intent(in) :: x, lmin, lmax
  double precision, intent(out), optional :: del_span
  double precision del0, del, y
  double precision, parameter :: del0_frac = 1D-4, del_ratio = 1.5D0
  integer i, n
  del0 = max((lmax - lmin) * del0_frac, grid_config%very_small_len*0.01D0)
  n = ceiling(log( &
         (lmax-lmin)/del0 * (del_ratio - 1D0) + 1D0) / log(del_ratio))
  get_int_val_along_y = 0D0
  y = lmin
  del = del0
  if (present(del_span)) then
    del_span = 0D0
  end if
  do i=1, n
    if (y .gt. lmax) then
      exit
    end if
    get_int_val_along_y = get_int_val_along_y + &
      get_density_analytic(x, y+0.5D0*del) * del
    if (present(del_span)) then
      del_span = del_span + del
    end if
    y = y + del
    del = del * del_ratio
  end do
end function get_int_val_along_y


function get_int_val_along_x(y, lmin, lmax, del_span)
  double precision get_int_val_along_x
  double precision, intent(in) :: y, lmin, lmax
  double precision, intent(out), optional :: del_span
  double precision del0, del, x
  double precision, parameter :: del0_frac = 1D-3, del_ratio = 2.0D0
  integer i, n
  del0 = max((lmax - lmin) * del0_frac, grid_config%very_small_len*0.01D0)
  n = ceiling(log( &
         (lmax-lmin)/del0 * (del_ratio - 1D0) + 1D0) / log(del_ratio))
  get_int_val_along_x = 0D0
  x = lmin
  del = del0
  if (present(del_span)) then
    del_span = 0D0
  end if
  do i=1, n
    if (x .gt. lmax) then
      exit
    end if
    get_int_val_along_x = get_int_val_along_x + &
      get_density_analytic(x+0.5D0*del, y) * del
    if (present(del_span)) then
      del_span = del_span + del
    end if
    x = x + del
    del = del * del_ratio
  end do
end function get_int_val_along_x


! Initialize a cell
recursive subroutine cell_init(c, parent, nChildren)
  type(type_cell), target :: c
  type(type_cell), target :: parent
  integer nChildren
  integer i
  if (.not. allocated(c%val)) then
    allocate(c%val(max(refinement_data%ncol - 2, 2)))
  end if
  c%order = parent%order + 1
  c%nChildren = nChildren
  c%parent => parent
  if (nChildren .gt. 0) then
    allocate(c%children(nChildren))
    do i=1, nChildren
      allocate(c%children(i)%p)
      call cell_init(c%children(i)%p, c, 0)
    end do
  end if
end subroutine cell_init


subroutine init_children(c, nChildren)
  type(type_cell), target :: c
  integer i, nChildren
  if (nChildren .gt. 0) then
    allocate(c%children(nChildren))
    do i=1, nChildren
      allocate(c%children(i)%p)
      call cell_init(c%children(i)%p, c, 0)
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
  double precision :: max_cell_size_frac = 0.1D0
  double precision :: max_cell_size_0 = 0.1D0
  double precision :: aspect_ratio_max = 20.0D0
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
  call init_children(c, c%nChildren)
  iSub = 0
  do i=0, ndiv_x-1
    do j=0, ndiv_y-1
      iSub = iSub + 1
      c%children(iSub)%p%xmin = c%xmin + del_x * dble(i)
      c%children(iSub)%p%xmax = c%children(iSub)%p%xmin + del_x
      c%children(iSub)%p%ymin = c%ymin + del_y * dble(j)
      c%children(iSub)%p%ymax = c%children(iSub)%p%ymin + del_y
      ! Get rid of the rounding errors.
      if (i .EQ. 0)          c%children(iSub)%p%xmin = c%xmin
      if (i .EQ. (ndiv_x-1)) c%children(iSub)%p%xmax = c%xmax
      if (j .EQ. 0)          c%children(iSub)%p%ymin = c%ymin
      if (j .EQ. (ndiv_y-1)) c%children(iSub)%p%ymax = c%ymax
    end do
  end do
end subroutine grid_decorate


function is_close_to_rect(v, xmin, xmax, ymin, ymax)
  logical is_close_to_rect
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  double precision, dimension(2) :: v
  double precision xtol, ytol
  xtol = 0.1D0 * (xmax - xmin)
  ytol = 0.1D0 * (ymax - ymin)
  if ( (v(1) .ge. xmin - xtol) .and. (v(1) .le. xmax + xtol) .and. &
       (v(2) .ge. ymin - ytol) .and. (v(2) .le. ymax + ytol)) then
    is_close_to_rect = .true.
  else
    is_close_to_rect = .false.
  end if
end function is_close_to_rect


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
  allocate(refinement_data%ave_val(refinement_data%ncol - 2))
  !call load_data_for_refinement_analytic
  call load_data_from_RADMC
  call init_interpolation
end subroutine load_data_for_refinement


subroutine load_data_for_refinement_analytic
  integer nx, ny
  double precision, dimension(:,:), allocatable :: x, y, z
  integer i, j
  double precision xmin, xmax, ymin, ymax
  double precision dx, dx_0
  double precision dy, dy_0
  !
  double precision :: dx_0_frac = 1D-4
  double precision :: dy_0_frac = 1D-4
  double precision :: dx_0_min_frac = 1D-1
  double precision :: dy_0_min_frac = 1D-1
  double precision :: dx_ratio = 1.05D0
  double precision :: dy_ratio = 1.05D0
  !
  xmin = 0D0
  xmax = 50D0
  ymin = 0D0
  ymax = 50D0
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
    allocate(d%xyv(2+1, d%nlen), d%idx_incell(d%nlen))
    d%xyv(1, :) = reshape(x, (/d%nlen/))
    d%xyv(2, :) = reshape(y, (/d%nlen/))
    d%xyv(3, :) = reshape(z, (/d%nlen/))
  end associate
  deallocate(x, y, z)
end subroutine load_data_for_refinement_analytic


function get_density_analytic(x, y)
  double precision get_density_analytic
  double precision, intent(in) :: x, y
  select case(grid_config%analytical_to_use)
    case ('Hayashi')
      get_density_analytic = density_analytic_Hayashi(x, y)
    case ('Andrew')
      get_density_analytic = density_analytic_Andrew(x, y)
  end select
end function get_density_analytic


! Merge method

! Neighborhood relationship


function density_analytic_Hayashi(r, z)
  ! Hayashi 1981, equations 2.7, 2.8, 2.9
  double precision, intent(in) :: r, z ! in AU
  double precision density_analytic_Hayashi ! in number cm-3
  double precision, parameter :: rho0 = 1.4D-9 ! g cm-3
  double precision, parameter :: n0 = rho0 / (1.4D0 * phy_mProton_CGS)
  double precision z0
  z0 = 0.0472D0 * r**1.25D0
  density_analytic_Hayashi = n0 * r**(-2.75D0) * exp(-(z/z0)**2)
end function density_analytic_Hayashi


function density_analytic_Andrew(r, z)
  ! Andrew 2009, equations 1, 2, ...
  ! The value of these parameters are chosen arbitrarily based on his table.
  use phy_const
  double precision, intent(in) :: r, z ! in AU
  double precision density_analytic_Andrew ! in number cm-3
  double precision theta, R_sph, sigma, h, RRc
  double precision, parameter :: MeanMolWeight = 1.4D0
  double precision, parameter :: Rc=150D0, H100=5D0, gam=1.0D0, psi=0.1D0, Md=0.1D0
  double precision, parameter :: sigma_c = (2D0 - gam) * Md / (phy_2Pi * Rc**2)
  double precision, parameter :: hc = H100 / (100D0 * (100D0/Rc)**psi)
  theta = phy_Pi_2 - atan2(z, r)
  R_sph = sqrt(r*r + z*z)
  RRc = R_sph / Rc
  h = hc * RRc**psi
  sigma = sigma_c * (RRc)**(-gam) * exp(-RRc**(2D0-gam))
  density_analytic_Andrew = sigma / (phy_sqrt2Pi * R_sph * h) * &
    exp(-0.5D0 * ((phy_Pi*0.5D0-theta)/h)**2) &
    * phy_Msun_CGS / (MeanMolWeight * phy_mProton_CGS) / (phy_AU2cm)**3
end function density_analytic_Andrew


function get_num_of_interval_log(xmin, xmax, dx0, del_ratio)
  integer get_num_of_interval_log
  double precision, intent(in) :: xmin, xmax, dx0, del_ratio
  get_num_of_interval_log = ceiling(log( &
         (xmax-xmin)/dx0 * (del_ratio - 1D0) + 1D0) / log(del_ratio))
end function get_num_of_interval_log


function get_ratio_of_interval_log(xmin, xmax, dx0, n, iout)
  double precision get_ratio_of_interval_log
  double precision, intent(in) :: xmin, xmax, dx0
  integer, intent(in) :: n
  integer, intent(out), optional :: iout
  integer :: i, imax = 100
  double precision :: frac_precision = 1D-3
  double precision k, p
  double precision tmp
  k = (xmax - xmin) / dx0
  p = 1D0 / dble(n)
  get_ratio_of_interval_log = &
    exp(p*log(k + 1D0)) - 1D0
  do i=1, imax
    tmp = get_ratio_of_interval_log
    get_ratio_of_interval_log = &
      exp(p * log(get_ratio_of_interval_log*k + 1D0)) - 1D0
    if (abs(tmp - get_ratio_of_interval_log) .le. get_ratio_of_interval_log * frac_precision) then
      exit
    end if
  end do
  get_ratio_of_interval_log = get_ratio_of_interval_log + 1D0
  if (present(iout)) then
    iout = i
  end if
end function get_ratio_of_interval_log


subroutine load_data_from_RADMC
  integer nx, ny, i, j
  character(len=256) pathname
  character(len=32) :: commentstr = ''
  pathname = combine_dir_filename(grid_config%data_dir, grid_config%data_filename)
  associate(d => refinement_data)
    call load_array_from_txt(pathname, d%xyv, d%ncol, d%nlen, d%nx, d%ny, commentstr)
  end associate
end subroutine load_data_from_RADMC


subroutine init_interpolation
  if (grid_config%interpolation_method .eq. 'barycentric') then
    call init_interpol_barycentric
  else if (grid_config%interpolation_method .eq. 'spline') then
    call init_interpol_spline
  else
    write(*,*) 'Unknown interpolation method.'
    stop
  end if
end subroutine init_interpolation


subroutine init_interpol_barycentric
  associate(d => refinement_data)
    allocate(n_bary2d, T_bary2d)
    n_bary2d%nx = d%nx
    n_bary2d%ny = d%ny
    n_bary2d%d  = 0
    allocate(n_bary2d%xi(n_bary2d%nx), &
             n_bary2d%yi(n_bary2d%ny), &
             n_bary2d%vi(n_bary2d%nx, n_bary2d%ny))
    n_bary2d%xi = d%xyv(1, 1:d%nlen:n_bary2d%ny) / phy_AU2cm
    n_bary2d%yi = d%xyv(2, 1:n_bary2d%ny)
    n_bary2d%vi = transpose(reshape(d%xyv(3, :), (/n_bary2d%ny, n_bary2d%nx/))) &
                  * RADMC_gas2dust_mass_ratio / (phy_mProton_CGS*MeanMolWeight)
    !
    T_bary2d%nx = n_bary2d%nx
    T_bary2d%ny = n_bary2d%ny
    T_bary2d%d = 0
    allocate(T_bary2d%xi(T_bary2d%nx), &
             T_bary2d%yi(T_bary2d%ny), &
             T_bary2d%vi(T_bary2d%nx, T_bary2d%ny))
    T_bary2d%xi = n_bary2d%xi
    T_bary2d%yi = n_bary2d%yi
    T_bary2d%vi = transpose(reshape(d%xyv(4, :), (/T_bary2d%ny, T_bary2d%nx/)))
    !
    write(*, '("Max n and T in RADMC", 2ES14.4)') maxval(n_bary2d%vi), maxval(T_bary2d%vi)
    !
    call barycentric2d_prepare(n_bary2d)
    call barycentric2d_prepare(T_bary2d)
  end associate
end subroutine init_interpol_barycentric


subroutine init_interpol_spline
  associate(d => refinement_data)
    allocate(n_spline2d, T_spline2d)
    n_spline2d%nx = d%nx
    n_spline2d%ny = d%ny
    n_spline2d%itype = 0
    allocate(n_spline2d%xi(n_spline2d%nx), &
             n_spline2d%yi(n_spline2d%ny), &
             n_spline2d%vi(n_spline2d%nx, n_spline2d%ny))
    n_spline2d%xi = d%xyv(1, 1:d%nlen:n_spline2d%ny) / phy_AU2cm
    n_spline2d%yi = d%xyv(2, 1:n_spline2d%ny)
    n_spline2d%vi = transpose(reshape(d%xyv(3, :), (/n_spline2d%ny, n_spline2d%nx/))) &
                    * RADMC_gas2dust_mass_ratio / (phy_mProton_CGS*MeanMolWeight)
    !
    T_spline2d%nx = n_spline2d%nx
    T_spline2d%ny = n_spline2d%ny
    T_spline2d%itype = 0
    allocate(T_spline2d%xi(T_spline2d%nx), &
             T_spline2d%yi(T_spline2d%ny), &
             T_spline2d%vi(T_spline2d%nx, T_spline2d%ny))
    T_spline2d%xi = n_spline2d%xi
    T_spline2d%yi = n_spline2d%yi
    T_spline2d%vi = transpose(reshape(d%xyv(4, :), (/T_spline2d%ny, T_spline2d%nx/)))
    !
    write(*, '("Max n and T in RADMC", 2ES14.4)') maxval(n_spline2d%vi), maxval(T_spline2d%vi)
    !
    call spline2d_prepare(n_spline2d)
    call spline2d_prepare(T_spline2d)
  end associate
end subroutine init_interpol_spline


function get_RADMC_n(x, y)
  double precision get_RADMC_n
  double precision, intent(in) :: x, y
  double precision y_theta
  y_theta = phy_Pi_2 - atan2(y, x)
  if (grid_config%interpolation_method .eq. 'barycentric') then
    if (y_theta .gt. n_bary2d%yi(n_bary2d%ny)) then
      y_theta = n_bary2d%yi(n_bary2d%ny)
    end if
    get_RADMC_n = barycentric2d_interpol(x, y_theta, n_bary2d)
  else if (grid_config%interpolation_method .eq. 'spline') then
    if (y_theta .gt. n_spline2d%yi(n_spline2d%ny)) then
      y_theta = n_spline2d%yi(n_spline2d%ny)
    end if
    get_RADMC_n = spline2d_interpol(x, y_theta, n_spline2d)
  else
    write(*,*) 'Unknown interpolation method.'
    stop
  end if
end function get_RADMC_n


function get_RADMC_T(x, y)
  double precision get_RADMC_T
  double precision, intent(in) :: x, y
  double precision y_theta
  y_theta = phy_Pi_2 - atan2(y, x)
  if (grid_config%interpolation_method .eq. 'barycentric') then
    if (y_theta .gt. T_bary2d%yi(T_bary2d%ny)) then
      y_theta = T_bary2d%yi(T_bary2d%ny)
    end if
    get_RADMC_T = barycentric2d_interpol(x, y_theta, T_bary2d)
  else if (grid_config%interpolation_method .eq. 'spline') then
    if (y_theta .gt. T_spline2d%yi(T_spline2d%ny)) then
      y_theta = T_spline2d%yi(T_spline2d%ny)
    end if
    get_RADMC_T = spline2d_interpol(x, y_theta, T_spline2d)
  else
    write(*,*) 'Unknown interpolation method.'
    stop
  end if
end function get_RADMC_T



end module grid
