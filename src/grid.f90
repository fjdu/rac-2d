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

type :: type_grid_config
  double precision rmin, rmax, zmin, zmax
  double precision rin, rout
  double precision :: dr0 = 0D0
  double precision :: ratiotatio = 1D0
  integer :: nratio = 1
  logical :: use_data_file_input = .false.
  logical :: columnwise = .true.
  logical :: refine_at_r0_in_exp = .true.
  integer :: ncol
  character(len=128) :: data_dir = './'
  character(len=128) :: data_filename = ''
  character(len=16) :: analytical_to_use = 'Andrews'
  character(len=16) :: interpolation_method = 'barycentric'
  double precision :: max_ratio_to_be_uniform = 2.0D0
  double precision :: density_log_range = 5D0, density_scale = 14.0D0
  double precision :: max_val_considered = 1D19
  double precision :: min_val_considered = 5D1
  double precision :: min_val_considered_use = 5D1
  double precision :: very_small_len = 1D-4
  double precision :: smallest_cell_size = 1D-2
  double precision :: largest_cell_size = 1D3
  double precision :: small_len_frac = 1D-2
end type type_grid_config

type :: type_leaves_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: idx
end type type_leaves_list

type :: type_Andrews_disk_p
  type(type_Andrews_disk), pointer :: p
end type type_Andrews_disk_p

type :: type_Andrews_disk_p_list
  integer n
  type(type_Andrews_disk_p), dimension(:), allocatable :: list
end type type_Andrews_disk_p_list

type(type_leaves) :: leaves

type(type_leaves) :: all_leaves

type(type_leaves_list) :: surf_cells, bott_cells

! Columns of cells
type(type_leaves), dimension(:), allocatable :: columns
type(type_simple_integer_list), dimension(:), allocatable :: columns_idx

type(type_grid_config) :: grid_config

private :: refinement_data
type(type_refinement_data) refinement_data

type(type_cell), pointer :: root

type(type_spline_2D), allocatable :: n_spline2d, T_spline2d
type(type_barycentric_2d), allocatable :: n_bary2d, T_bary2d

double precision, parameter, private :: MeanMolWeight = 1.4D0
double precision, parameter, private :: RADMC_gas2dust_mass_ratio = 1D2 ! Todo

type(type_Andrews_disk_p_list) :: andrews_list

type(type_Andrews_disk) a_andrews_4ini

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
    call grid_init_columnwise_new(root)
    !call grid_init_columnwise(root)
    !call grid_init_columnwise_alt(root)
    do i=1, root%nChildren
      call grid_refine(root%children(i)%p)
      root%nOffspring = root%nOffspring + root%children(i)%p%nOffspring + 1
      root%nleaves = root%nleaves + root%children(i)%p%nleaves
    end do
  else
    call grid_init(root)
    call grid_refine(root)
  end if
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
  call grid_make_post
end subroutine make_grid


subroutine grid_make_post
  if (allocated(refinement_data%xyv)) then
    deallocate(refinement_data%xyv, &
        refinement_data%idx_incell, refinement_data%ave_val)
  end if
  !
  if (allocated(n_spline2d)) then
    deallocate(n_spline2d%xi, n_spline2d%yi, n_spline2d%vi)
    deallocate(T_spline2d%xi, T_spline2d%yi, T_spline2d%vi)
    deallocate(n_spline2d, T_spline2d)
  end if
  !
  if (allocated(n_bary2d)) then
    deallocate(n_bary2d%xi, n_bary2d%yi, n_bary2d%vi)
    deallocate(T_bary2d%xi, T_bary2d%yi, T_bary2d%vi)
    deallocate(n_bary2d, T_bary2d)
  end if
end subroutine grid_make_post


subroutine make_all_leaves(c)
  type(type_cell), target, intent(in) :: c
  integer idx
  all_leaves%nlen = c%nOffspring
  allocate(all_leaves%list(all_leaves%nlen))
  idx = 0
  call add_leaves_all(c, idx)
end subroutine make_all_leaves


recursive subroutine add_leaves_all(c, idx)
  integer i
  integer, intent(inout) :: idx
  type(type_cell), target :: c
  do i=1, c%nChildren
    idx = idx + 1
    all_leaves%list(idx)%p => c%children(i)%p
    call add_leaves_all(c%children(i)%p, idx)
  end do
end subroutine add_leaves_all


subroutine grid_make_surf_bott
  integer i, j, i1, j1
  integer, parameter :: nSurfBott_max = 1024
  integer, dimension(nSurfBott_max) :: idxSurf, idxBott
  !double precision, parameter :: fra_tot_thre = 0.9D0
  surf_cells%nlen = 0
  bott_cells%nlen = 0
  do i=1, leaves%nlen
    if (leaves%list(i)%p%above%n .eq. 0) then
    !if (leaves%list(i)%p%above%fra_tot .lt. fra_tot_thre) then
      surf_cells%nlen = surf_cells%nlen + 1
      idxSurf(surf_cells%nlen) = i
    end if
    if (leaves%list(i)%p%below%n .eq. 0) then
    !if (leaves%list(i)%p%below%fra_tot .lt. fra_tot_thre) then
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
  !
  ! Sort: in -> out
  do i=1, surf_cells%nlen
    do j=1, i-1
      i1 = surf_cells%idx(i)
      j1 = surf_cells%idx(j)
      if (leaves%list(j1)%p%xmin .gt. leaves%list(i1)%p%xmin) then
        surf_cells%idx(i) = j1
        surf_cells%idx(j) = i1
      end if
    end do
  end do
  do i=1, bott_cells%nlen
    do j=1, i-1
      i1 = bott_cells%idx(i)
      j1 = bott_cells%idx(j)
      if (leaves%list(j1)%p%xmin .gt. leaves%list(i1)%p%xmin) then
        bott_cells%idx(i) = j1
        bott_cells%idx(j) = i1
      end if
    end do
  end do
end subroutine grid_make_surf_bott


subroutine grid_make_leaves(croot)
  type(type_cell), target, intent(in) :: croot
  integer i, idx
  !
  if (allocated(leaves%list)) then
    do i=1, leaves%nlen
      if (associated(leaves%list(i)%p)) then
        leaves%list(i)%p%id = -1
      end if
      nullify(leaves%list(i)%p)
    end do
    deallocate(leaves%list)
  end if
  !
  leaves%nlen = croot%nleaves
  !
  allocate(leaves%list(leaves%nlen))
  do i=1, leaves%nlen
    nullify(leaves%list(i)%p)
  end do
  !
  idx = 0
  call grid_add_leaves(croot, idx)
  !
end subroutine grid_make_leaves


recursive subroutine get_number_of_leaves(c)
  integer i
  type(type_cell), target :: c
  !
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
  integer i
  integer, intent(inout) :: idx
  type(type_cell), target :: c
  if (c%using) then
    idx = idx + 1
    leaves%list(idx)%p => c
    leaves%list(idx)%p%id = idx
    return
  else
    do i=1, c%nChildren
      !write(*,*) i, c%nChildren, size(c%children), c%order, associated(c%children)
      !write(*,*) 'x', associated(c%children(i)%p)
      call grid_add_leaves(c%children(i)%p, idx)
    end do
  end if
end subroutine grid_add_leaves


subroutine grid_make_neighbors
  integer i
  integer nnei_max, idx
  nnei_max = 0
  idx = 1
  do i=1, leaves%nlen
    call make_neighbors(i)
    if (associated(leaves%list(i)%p%around)) then
      if (nnei_max .lt. leaves%list(i)%p%around%n) then
        nnei_max = leaves%list(i)%p%around%n
        idx = i
      end if
    end if
  end do
  write(*, '(A, I6)') 'Max number of neighbors:', nnei_max
  write(*, '(A, 4ES12.4)') 'owned by this one: xmin,xmax,ymin,ymax', &
    leaves%list(idx)%p%xmin, leaves%list(idx)%p%xmax, leaves%list(idx)%p%ymin, leaves%list(idx)%p%ymax
end subroutine grid_make_neighbors


function get_min_cell_size()
  double precision get_min_cell_size
  double precision tmp
  integer i
  get_min_cell_size = 1D99
  do i=1, leaves%nlen
    tmp = min(leaves%list(i)%p%xmax - leaves%list(i)%p%xmin, &
              leaves%list(i)%p%ymax - leaves%list(i)%p%ymin)
    if (tmp .lt. get_min_cell_size) then
      get_min_cell_size = tmp
    end if
  end do
end function get_min_cell_size


function get_ave_val_analytic(xmin, xmax, ymin, ymax, andrews)
  double precision get_ave_val_analytic
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  type(type_Andrews_disk), intent(in), optional :: andrews
  integer :: i, j, nx
  double precision dx, dx0, x, xx, area, dely
  double precision, parameter :: dx0_frac = 1D-5, dy0_frac = 1D-5, &
    dx_ratio = 1.2D0, dy_ratio=1.2D0
  dx0 = min((xmax - xmin) * 0.1D0, &
            max((xmax - xmin) * dx0_frac, &
                grid_config%very_small_len*0.01D0))
  nx = ceiling(log( &
         (xmax-xmin)/dx0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio)) + 5
  get_ave_val_analytic = 0D0
  area = 0D0
  x = xmin
  dx = dx0
  do i=1, nx
    xx = x + 0.5D0 * dx
    if (xx .gt. xmax) then
      exit
    end if
    !
    if (present(andrews)) then
      get_ave_val_analytic = get_ave_val_analytic + &
        get_int_val_along_y(xx, ymin, ymax, dely, andrews) * dx
    else
      get_ave_val_analytic = get_ave_val_analytic + &
        get_int_val_along_y(xx, ymin, ymax, dely) * dx
    end if
    area = area + dx * dely
    x = x + dx
    dx = dx * dx_ratio
  end do
  if (area .eq. 0D0) then
    write(*, '(A)') 'In get_ave_val_analytic:'
    write(*, '(A, 4F16.10)') 'Area  = 0!', xmin, xmax, dx0, dely
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
  double precision x, dx, del, tmp, ymaxtmp
  integer i
  integer nncol
  double precision, parameter :: frac_div_min = 1D-4
  !
  nncol = ceiling(real(grid_config%ncol) / real(grid_config%nratio))
  !
  if (grid_config%ncol .ne. (grid_config%nratio * nncol)) then
    grid_config%ncol = grid_config%nratio * nncol
    write(*, '(A, I5)') 'Changing number of columns into ', grid_config%ncol
  end if
  !
  c%xmin = grid_config%rmin
  c%xmax = grid_config%rmax
  c%ymin = grid_config%zmin
  c%ymax = grid_config%zmax
  c%nChildren = grid_config%ncol * 2
  !dx0 = max(1D-2, (c%xmax - c%xmin) * 2D-4)
  if (grid_config%dr0 .gt. 0D0) then
    dx0 = grid_config%dr0
    del_ratio = get_ratio_of_interval_log(c%xmin, c%xmax, &
                                          dx0, grid_config%ncol)
  else
    tmp = log(c%xmax/c%xmin) / dble(grid_config%ncol)
    if (tmp .le. 1D-4) then
      dx0 = tmp * c%xmin
    else
      dx0 = (exp(tmp) - 1D0) * c%xmin
    end if
    del_ratio = exp(tmp)
  end if
  call init_children(c, c%nChildren)
  dx = dx0
  x = c%xmin
  do i=1, grid_config%ncol
    c%children(i)%p%xmin = x
    c%children(i)%p%xmax = x + dx
    c%children(i)%p%ymin = c%ymin
    ymaxtmp = get_ymax_here(x + 0.5D0 * dx, c%ymin, c%ymax, frac_div_min)
    if ((ymaxtmp - c%ymin) .le. grid_config%smallest_cell_size) then
      ymaxtmp = get_ymax_here(x + 0.5D0 * dx, c%ymin, &
                c%ymin + 5D0*frac_div_min * (c%ymax-c%ymin), frac_div_min)
    end if
    if (ymaxtmp .ge. c%ymax/1.5D0) then
      c%children(i)%p%ymax = c%ymax/1.5D0
    else
      c%children(i)%p%ymax = ymaxtmp
    end if
    if (c%children(i)%p%ymax - c%children(i)%p%ymin .lt. &
        grid_config%smallest_cell_size) then
      write(*, '(A, 4F16.10)') 'Short column:', &
        c%children(i)%p%xmin, c%children(i)%p%xmax, &
        c%children(i)%p%ymin, c%children(i)%p%ymax
      c%children(i)%p%ymax = c%children(i)%p%ymin + &
        grid_config%smallest_cell_size * 4D0
      write(*, '(A, F16.10)') 'Set ymax to:', c%children(i)%p%ymax
    end if
    x = x + dx
    dx = dx * del_ratio
  end do
  c%children(grid_config%ncol)%p%xmax = c%xmax
  !
  do i=1+grid_config%ncol, grid_config%ncol*2
    c%children(i)%p%xmin = c%children(i-grid_config%ncol)%p%xmin
    c%children(i)%p%xmax = c%children(i-grid_config%ncol)%p%xmax
    c%children(i)%p%ymin = c%children(i-grid_config%ncol)%p%ymax
    c%children(i)%p%ymax = c%ymax
  end do
  !
  !c%children(1)%p%xmin = c%xmin
  !c%children(1)%p%xmax = c%xmin + 1D-2
  !c%children(2)%p%xmin = c%children(1)%p%xmax
  !c%children(2)%p%xmax = c%children(2)%p%xmin + 1D-2
  !c%children(3)%p%xmin = c%children(2)%p%xmax
  !c%children(grid_config%ncol+1)%p%xmin = c%children(1)%p%xmin
  !c%children(grid_config%ncol+1)%p%xmax = c%children(1)%p%xmax
  !c%children(grid_config%ncol+2)%p%xmin = c%children(2)%p%xmin
  !c%children(grid_config%ncol+2)%p%xmax = c%children(2)%p%xmax
  !c%children(grid_config%ncol+3)%p%xmin = c%children(3)%p%xmin
  !
end subroutine grid_init_columnwise


subroutine grid_init_columnwise_new(c)
  type(type_cell), target :: c
  double precision dx0, del_ratio
  double precision x, dx, ymaxtmp
  integer i
  integer nncol
  double precision, dimension(:), allocatable :: locs
  double precision, parameter :: frac_div_min = 1D-4
  !
  nncol = ceiling(real(grid_config%ncol) / real(grid_config%nratio))
  !
  if (grid_config%ncol .ne. (grid_config%nratio * nncol)) then
    grid_config%ncol = grid_config%nratio * nncol
    write(*, '(A, I5)') 'Changing number of columns into ', grid_config%ncol
  end if
  !
  call cell_init(c)
  !
  c%xmin = grid_config%rmin
  c%xmax = grid_config%rmax
  c%ymin = grid_config%zmin
  c%ymax = grid_config%zmax
  c%nChildren = grid_config%ncol * 2
  !
  call init_children(c, c%nChildren)
  !
  allocate(locs(1+grid_config%ncol))
  locs = 0D0
  call get_column_locations(1+grid_config%ncol, locs)
  !
  do i=1, grid_config%ncol
    c%children(i)%p%xmin = locs(i)
    c%children(i)%p%xmax = locs(i+1)
    !write(*, '(A, I6, 2ES16.6)') 'Column', i, c%children(i)%p%xmin, c%children(i)%p%xmax
    !
    x = locs(i)
    dx = locs(i+1) - locs(i)
    !
    c%children(i)%p%ymin = c%ymin
    !
    ymaxtmp = get_ymax_here(x + 0.5D0 * dx, c%ymin, c%ymax, frac_div_min)
    if ((ymaxtmp - c%ymin) .le. grid_config%smallest_cell_size) then
      ymaxtmp = get_ymax_here(x + 0.5D0 * dx, c%ymin, &
                c%ymin + 5D0*frac_div_min * (c%ymax-c%ymin), frac_div_min)
    end if
    if (ymaxtmp .ge. c%ymax/1.5D0) then
      c%children(i)%p%ymax = c%ymax/1.5D0
    else
      c%children(i)%p%ymax = ymaxtmp
    end if
    if (c%children(i)%p%ymax - c%children(i)%p%ymin .lt. &
        grid_config%smallest_cell_size) then
      write(*, '(A, 4F16.10)') 'Short column:', &
        c%children(i)%p%xmin, c%children(i)%p%xmax, &
        c%children(i)%p%ymin, c%children(i)%p%ymax
      c%children(i)%p%ymax = c%children(i)%p%ymin + &
        grid_config%smallest_cell_size * 4D0
      write(*, '(A, F16.10)') 'Set ymax to:', c%children(i)%p%ymax
    end if
  end do
  c%children(grid_config%ncol)%p%xmax = c%xmax
  !
  do i=1+grid_config%ncol, grid_config%ncol*2
    c%children(i)%p%xmin = c%children(i-grid_config%ncol)%p%xmin
    c%children(i)%p%xmax = c%children(i-grid_config%ncol)%p%xmax
    c%children(i)%p%ymin = c%children(i-grid_config%ncol)%p%ymax
    c%children(i)%p%ymax = c%ymax
  end do
  !
  if (allocated(locs)) then
    deallocate(locs)
  end if
  !
  !do i=1, c%nChildren
  !  write(*, '(I5, 4ES12.4)') &
  !    i, c%children(i)%p%xmin, c%children(i)%p%xmax, c%children(i)%p%ymin, c%children(i)%p%ymax
  !end do
  !
end subroutine grid_init_columnwise_new


subroutine grid_init_columnwise_alt(c)
  type(type_cell), target :: c
  double precision x, dx, dx0, del, ratio
  integer, parameter :: N = 500
  double precision, dimension(N) :: xs, xx, dens, ncol
  integer i, j, ic, idx
  double precision acc, Ncol_threshold
  !
  dx0 = 1D-4
  Ncol_threshold = 1D23
  !
  c%xmin = grid_config%rmin
  c%xmax = grid_config%rmax
  c%ymin = grid_config%zmin
  c%ymax = grid_config%zmax
  !
  ratio = get_ratio_of_interval_log(c%xmin, c%xmax, dx0, N)
  xs(1) = c%xmin
  dens(1) = get_RADMC_n(xs(1), c%ymin)
  ncol(1) = 0D0
  dx = dx0
  do i=2, N
    xs(i) = xs(i-1) + dx
    dens(i) = get_RADMC_n(xs(i), c%ymin)
    ncol(i) = dx * phy_AU2cm * (dens(i-1) + dens(i)) * 0.5D0
    dx = dx * ratio
  end do
  if (xs(N-1) .lt. c%xmax) then
    xs(N) = c%xmax
  end if
  idx = 1
  ic = 0
  do i=1, N
    acc = 0D0
    do j=idx, N
      acc = acc + ncol(j)
      if (acc .ge. Ncol_threshold) then
        ic = ic + 1
        xx(ic) = xs(j)
        idx = j + 1
        exit
      end if
    end do
    if (idx .gt. N) then
      exit
    end if
  end do
  c%nChildren = ic
  call init_children(c, c%nChildren)
  do i=1, c%nChildren
    if (i .gt. 1) then
      c%children(i)%p%xmin = c%children(i-1)%p%xmax
    else
      c%children(i)%p%xmin = c%xmin
    end if
    c%children(i)%p%xmax = xx(i)
    c%children(i)%p%ymin = c%ymin
    c%children(i)%p%ymax = c%ymax
  end do
end subroutine grid_init_columnwise_alt



subroutine get_column_locations(n, locs)
  integer, intent(in) :: n
  double precision, dimension(n), intent(out) :: locs
  double precision r0, tmp, delr ,delr1
  integer n1, n2, n3
  r0 = a_andrews_4ini%r0_in_exp
  if (((grid_config%rmin .ge. r0) .or. &
       (grid_config%rmax .le. r0)) .or. &
      (.not. grid_config%refine_at_r0_in_exp)) then
    !
    call logspace(locs, log10(grid_config%rmin), log10(grid_config%rmax), &
        grid_config%ncol+1, 10D0)
    !write(*,*) locs(1:4), locs(grid_config%ncol-3:grid_config%ncol)
    !
  else
    tmp = sqrt(grid_config%rmax*grid_config%rmin/r0/r0)
    !
    n1 = ceiling(real(grid_config%ncol) * 0.8 / (0.8 + tmp))
    n2 = ceiling(real(grid_config%ncol) * tmp / (0.8 + tmp) * 0.2)
    n3 = grid_config%ncol+1 - n1 - n2
    if (n1*n2*n3 .eq. 0) then
      write(*, '(A)') 'In get_column_locations:'
      write(*, '(A, 3I8)') 'n1*n2*n3 .eq. 0', n1, n2, n3
      stop
    end if
    !
    delr  = r0 * 8e-2
    delr1 = r0 * 1e-3
    !
    call logspace(locs(1 : n1), log10(grid_config%rmin), log10(r0-delr1), n1, 10D0)
    call logspace(locs(n1 : (n1+n2)), log10(r0-delr1), log10(r0+delr), n2+1, 10D0)
    call logspace(locs((n1+n2) : (n1+n2+n3)), log10(r0+delr), &
        log10(grid_config%rmax), n3+1, 10D0)
  end if

end subroutine get_column_locations



function get_ymax_here(x, y0, y1, fracmin)
  double precision get_ymax_here
  double precision, intent(in) :: x, y0, y1
  double precision, intent(in), optional :: fracmin
  integer :: i, n = 100
  double precision y, dy, del_ratio, val, frac
  get_ymax_here = 0D0
  if (present(fracmin)) then
    frac = fracmin
  else
    frac = 1D-4
  end if
  dy = (y1 - y0) * frac
  del_ratio = get_ratio_of_interval_log(y0, y1, dy, n)
  dy = dy * del_ratio**(n-1)
  y = y1
  do i=1, n
    if (grid_config%use_data_file_input) then
      val = get_RADMC_n(x, y)
    else
      val = get_density_analytic(x, y)
    end if
    if (val .ge. grid_config%min_val_considered) then
      get_ymax_here = y
      return
    end if
    y = y - dy
    dy = dy / del_ratio
  end do
end function get_ymax_here


recursive subroutine reset_cell_bounds(c)
  integer i
  type(type_cell), intent(inout) :: c
  if (c%nChildren .gt. 0) then
    c%ymax = 0D0
    c%ymin = 1D99
    do i=1, c%nChildren
      call reset_cell_bounds(c%children(i)%p)
      c%ymax = max(c%ymax, c%children(i)%p%ymax)
      c%ymin = min(c%ymin, c%children(i)%p%ymin)
    end do
  end if
end subroutine reset_cell_bounds



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
      c%nOffspring = 0
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


subroutine make_neighbors(id)
  integer, intent(in) :: id
  type(type_cell), pointer :: c
  integer i
  integer n, pos
  integer, parameter :: nnei_max = 512
  integer, dimension(:), allocatable :: idx_inner, idx_outer, idx_below, idx_above
  !double precision, dimension(:), allocatable :: fra_inner, fra_outer, fra_below, fra_above
  double precision frac, tol
  !
  allocate(idx_inner(nnei_max), &
           idx_outer(nnei_max), &
           idx_below(nnei_max), &
           idx_above(nnei_max)  &
           !fra_inner(nnei_max), &
           !fra_outer(nnei_max), &
           !fra_below(nnei_max), &
           !fra_above(nnei_max) &
          )
  !
  c => leaves%list(id)%p
  if (.not. associated(c%inner)) then
    allocate(c%inner, c%outer, c%below, c%above, c%around)
  else
    if (allocated(c%inner%idx)) then
      deallocate(c%inner%idx)
      !deallocate(c%inner%idx, c%inner%fra)
    end if
    if (allocated(c%outer%idx)) then
      deallocate(c%outer%idx)
      !deallocate(c%outer%idx, c%outer%fra)
    end if
    if (allocated(c%above%idx)) then
      deallocate(c%above%idx)
      !deallocate(c%above%idx, c%above%fra)
    end if
    if (allocated(c%below%idx)) then
      deallocate(c%below%idx)
      !deallocate(c%below%idx, c%below%fra)
    end if
    if (allocated(c%around%idx)) then
      deallocate(c%around%idx)
      !deallocate(c%around%idx, c%around%fra)
    end if
  end if
  c%inner%n = 0
  c%outer%n = 0
  c%above%n = 0
  c%below%n = 0
  c%around%n = 0
  do i=1, leaves%nlen
    if (i .eq. id) then
      cycle
    end if
    tol = min(1D-3 * min(c%xmax-c%xmin, c%ymax-c%ymin, &
        leaves%list(i)%p%xmax - leaves%list(i)%p%xmin, &
        leaves%list(i)%p%ymax - leaves%list(i)%p%ymin), &
        grid_config%very_small_len)
    if (is_neighbor(c, leaves%list(i)%p, tol, pos, frac)) then
      select case(pos)
        case(1)
          c%inner%n = c%inner%n + 1
          if (c%inner%n .gt. nnei_max) then
            write(*, '(A)') 'In make_neighbors:'
            write(*, '(A, 2I8)') 'c%inner%n too large:', c%inner%n, nnei_max
            stop
          end if
          idx_inner(c%inner%n) = i
          !fra_inner(c%inner%n) = frac
        case(2)
          c%outer%n = c%outer%n + 1
          if (c%outer%n .gt. nnei_max) then
            write(*, '(A)') 'In make_neighbors:'
            write(*, '(A, 2I8)') 'c%outer%n too large:', c%outer%n, nnei_max
            stop
          end if
          idx_outer(c%outer%n) = i
          !fra_outer(c%outer%n) = frac
        case(3)
          c%below%n = c%below%n + 1
          if (c%below%n .gt. nnei_max) then
            write(*, '(A)') 'In make_neighbors:'
            write(*, '(A, 2I8)') 'c%below%n too large:', c%below%n, nnei_max
            stop
          end if
          idx_below(c%below%n) = i
          !fra_below(c%below%n) = frac
        case(4)
          c%above%n = c%above%n + 1
          if (c%above%n .gt. nnei_max) then
            write(*, '(A)') 'In make_neighbors:'
            write(*, '(A, 2I8)') 'c%above%n too large:', c%above%n, nnei_max
            stop
          end if
          idx_above(c%above%n) = i
          !fra_above(c%above%n) = frac
      end select
    end if
  end do
  c%around%n = c%inner%n + c%outer%n + c%below%n + c%above%n
  if (c%around%n .gt. 0) then
    !allocate(c%around%idx(c%around%n), c%around%fra(c%around%n))
    allocate(c%around%idx(c%around%n))
  end if
  i = 0
  if (c%above%n .gt. 0) then
    !allocate(c%above%idx(c%above%n), c%above%fra(c%above%n))
    allocate(c%above%idx(c%above%n))
    c%above%idx = idx_above(1:c%above%n)
    !c%above%fra = fra_above(1:c%above%n)
    !c%above%fra_tot = sum(fra_above(1:c%above%n))
    c%around%idx(i+1:i+c%above%n) = c%above%idx
    !c%around%fra(i+1:i+c%above%n) = c%above%fra
    i = i + c%above%n
  end if
  if (c%inner%n .gt. 0) then
    !allocate(c%inner%idx(c%inner%n), c%inner%fra(c%inner%n))
    allocate(c%inner%idx(c%inner%n))
    c%inner%idx = idx_inner(1:c%inner%n)
    !c%inner%fra = fra_inner(1:c%inner%n)
    !c%inner%fra_tot = sum(fra_inner(1:c%inner%n))
    c%around%idx(i+1:i+c%inner%n) = c%inner%idx
    !c%around%fra(i+1:i+c%inner%n) = c%inner%fra
    i = i + c%inner%n
  end if
  if (c%outer%n .gt. 0) then
    !allocate(c%outer%idx(c%outer%n), c%outer%fra(c%outer%n))
    allocate(c%outer%idx(c%outer%n))
    c%outer%idx = idx_outer(1:c%outer%n)
    !c%outer%fra = fra_outer(1:c%outer%n)
    !c%outer%fra_tot = sum(fra_outer(1:c%outer%n))
    c%around%idx(i+1:i+c%outer%n) = c%outer%idx
    !c%around%fra(i+1:i+c%outer%n) = c%outer%fra
    i = i + c%outer%n
  end if
  if (c%below%n .gt. 0) then
    !allocate(c%below%idx(c%below%n), c%below%fra(c%below%n))
    allocate(c%below%idx(c%below%n))
    c%below%idx = idx_below(1:c%below%n)
    !c%below%fra = fra_below(1:c%below%n)
    !c%below%fra_tot = sum(fra_below(1:c%below%n))
    c%around%idx(i+1:i+c%below%n) = c%below%idx
    !c%around%fra(i+1:i+c%below%n) = c%below%fra
    i = i + c%below%n
  end if
  deallocate(idx_inner, &
             idx_outer, &
             idx_below, &
             idx_above  &
             !fra_inner, &
             !fra_outer, &
             !fra_below, &
             !fra_above &
            )
end subroutine make_neighbors


function is_neighbor(c1, c2, tol, pos, frac)
  logical is_neighbor
  type(type_cell), target, intent(in) :: c1, c2
  double precision, intent(in) :: tol
  integer, intent(out) :: pos
  double precision, intent(out) :: frac
  double precision, parameter :: frac_too_small = 1D-6
  is_neighbor = .true.
  if ((abs(c2%xmax - c1%xmin) .le. tol) .and. &
      Two_lines_has_common_section(c1%ymin, c1%ymax, c2%ymin, c2%ymax, frac)) then
    pos = 1
    frac = frac / (c1%ymax - c1%ymin)
  else if ((abs(c2%xmin - c1%xmax) .le. tol) .and. &
      Two_lines_has_common_section(c1%ymin, c1%ymax, c2%ymin, c2%ymax, frac)) then
    pos = 2
    frac = frac / (c1%ymax - c1%ymin)
  else if ((abs(c2%ymax - c1%ymin) .le. tol) .and. &
      Two_lines_has_common_section(c1%xmin, c1%xmax, c2%xmin, c2%xmax, frac)) then
    pos = 3
    frac = frac / (c1%xmax - c1%xmin)
  else if ((abs(c2%ymin - c1%ymax) .le. tol) .and. &
      Two_lines_has_common_section(c1%xmin, c1%xmax, c2%xmin, c2%xmax, frac)) then
    pos = 4
    frac = frac / (c1%xmax - c1%xmin)
  else
    is_neighbor = .false.
    pos = -1
    frac = 0D0
  end if
  if (frac .le. frac_too_small) then
    is_neighbor = .false.
    pos = -1
    frac = 0D0
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
  ! Such conditional statements are always dangerous!  Todo
  if (grid_config%columnwise) then
    if ((c%ymax-c%ymin) .gt. grid_config%largest_cell_size) then
      is_uniform = .false.
      return
    end if
    is_uniform = test_uniformity_columnwise(c%xmin, c%xmax, c%ymin, c%ymax)
  else
    if (min(c%xmax-c%xmin, c%ymax-c%ymin) .gt. &
        grid_config%largest_cell_size) then
      is_uniform = .false.
      return
    end if
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
                      grid_config%smallest_cell_size)
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
                      grid_config%smallest_cell_size)
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
  double precision minv, maxv, max_ratio_to_be_uniform_here!, Av_vert_0
  integer i
  integer, parameter :: n = 4
  double precision, dimension(n) :: vals
  double precision, parameter :: const_small_num = 1D-100
  vals(1) = get_density_analytic(xmin, ymin)
  vals(2) = get_density_analytic(xmax, ymin)
  vals(3) = get_density_analytic(xmax, ymax)
  vals(4) = get_density_analytic(xmin, ymax)
  maxv = maxval(vals)
  minv = minval(vals)
  max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
    ((log10(maxv) - grid_config%density_scale)/grid_config%density_log_range)**2
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
  double precision minv, maxv, max_ratio_to_be_uniform_here!, Av_vert_0
  integer i
  integer, parameter :: n = 3
  double precision, dimension(n) :: vals
  double precision, parameter :: const_small_num = 1D-100
  vals(1) = get_density_analytic(0.5D0*(xmin+xmax), ymin)
  vals(2) = get_density_analytic(0.5D0*(xmin+xmax), 0.5D0*(ymin+ymax))
  vals(3) = get_density_analytic(0.5D0*(xmin+xmax), ymax)
  maxv = maxval(vals)
  minv = minval(vals)
  max_ratio_to_be_uniform_here = grid_config%max_ratio_to_be_uniform + &
    ((log10(maxv) - grid_config%density_scale)/grid_config%density_log_range)**2
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
        ((log10(maxv) - grid_config%density_scale)/grid_config%density_log_range)**2
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
  integer i
  double precision dely, y
  integer, parameter :: n = 10
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
        ((log10(maxv) - grid_config%density_scale)/grid_config%density_log_range)**2
      if (maxv / (abs(minv) + const_small_num) .le. max_ratio_to_be_uniform_here) then
        test_uniformity_based_on_data_columnwise = .true.
      else
        test_uniformity_based_on_data_columnwise = .false.
      end if
    end if
    !write(*,*) xmin, xmax, ymin, ymax, maxv, minv, test_uniformity_based_on_data_columnwise
  end associate
end function test_uniformity_based_on_data_columnwise


function get_int_val_along_y(x, lmin, lmax, del_span, andrews)
  double precision get_int_val_along_y
  double precision, intent(in) :: x, lmin, lmax
  type(type_Andrews_disk), intent(in), optional :: andrews
  double precision, intent(out), optional :: del_span
  double precision del0, del, y, yy
  double precision, parameter :: del0_frac = 1D-4, del_ratio = 1.3D0
  integer i, n
  del0 = min((lmax-lmin)*0.1D0, max((lmax - lmin) * del0_frac, grid_config%very_small_len*0.01D0))
  n = ceiling(log( &
         (lmax-lmin)/del0 * (del_ratio - 1D0) + 1D0) / log(del_ratio)) + 5
  get_int_val_along_y = 0D0
  y = lmin
  del = del0
  if (present(del_span)) then
    del_span = 0D0
  end if
  do i=1, n
    yy = y + 0.5D0 * del
    if (yy .gt. lmax) then
      exit
    end if
    if (present(andrews)) then
      get_int_val_along_y = get_int_val_along_y + &
        get_density_analytic(x, yy, andrews) * del
    else
      get_int_val_along_y = get_int_val_along_y + &
        get_density_analytic(x, yy) * del
    end if
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
  type(type_cell), target, intent(inout) :: c
  type(type_cell), target, intent(in), optional :: parent
  integer, intent(in), optional :: nChildren
  integer i, stat
  if (present(parent)) then
    c%order = parent%order + 1
    c%parent => parent
  else
    c%order = 0
    nullify(c%parent)
  end if
  if (present(nChildren)) then
    c%nChildren = nChildren
    if (nChildren .gt. 0) then
      allocate(c%children(nChildren))
      do i=1, nChildren
        allocate(c%children(i)%p)
        call cell_init(c%children(i)%p, c, 0)
      end do
    end if
  else
    c%nChildren = 0
  end if
  c%using = .false.
end subroutine cell_init


subroutine init_children(c, nChildren)
  type(type_cell), target :: c
  integer i, nChildren
  if (nChildren .gt. 0) then
    if (allocated(c%children)) then
      if (size(c%children) .eq. nChildren) then
        return
      else
        deallocate(c%children)
      end if
    end if
    allocate(c%children(nChildren))
    do i=1, nChildren
      nullify(c%children(i)%p)
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
  double precision max_cell_size_frac, max_cell_size_0, aspect_ratio_max
  !
  max_cell_size_frac = 0.1D0
  max_cell_size_0 = 0.1D0
  aspect_ratio_max = 20.0D0
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
  double precision :: dx_0_frac, dy_0_frac, dx_0_min_frac, dy_0_min_frac, &
    dx_ratio, dy_ratio
  !
  dx_0_frac = 1D-4
  dy_0_frac = 1D-4
  dx_0_min_frac = 1D-1
  dy_0_min_frac = 1D-1
  dx_ratio = 1.05D0
  dy_ratio = 1.05D0
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


function get_density_analytic(x, y, andrews)
  double precision get_density_analytic
  double precision, intent(in) :: x, y
  type(type_Andrews_disk), intent(in), optional :: andrews
  select case(grid_config%analytical_to_use)
    case ('Hayashi')
      get_density_analytic = density_analytic_Hayashi(x, y)
    case ('Andrews')
      if (present(andrews)) then
        get_density_analytic = density_analytic_Andrews(x, y, andrews)
      else
        get_density_analytic = density_analytic_Andrews(x, y)
      end if
    case default
      write(*, '(/A)') 'In get_density_analytic:'
      write(*, '(A/)') 'Should not have this case!'
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


function density_analytic_Andrews(r, z, andrews)
  double precision, intent(in) :: r, z ! in AU
  type(type_Andrews_disk), intent(in), optional :: andrews
  double precision density_analytic_Andrews ! in number cm-3
  if (present(andrews)) then
    density_analytic_Andrews = Andrews_dens(r, z, andrews)
  else
    density_analytic_Andrews = Andrews_dens(r, z, a_andrews_4ini)
  end if
end function density_analytic_Andrews



pure function Andrews_dens(r, z, andrews)
  ! Andrews 2009, equations 1, 2, ...
  ! The value of these parameters are chosen arbitrarily based on his table.
  double precision Andrews_dens ! in number cm-3
  double precision, intent(in) :: r, z ! in AU
  type(type_Andrews_disk), intent(in) :: andrews
  !
  double precision sigma_c, sigma, rrc, h
  double precision tmp1, tmp2, rlog
  double precision Md, rc, hc, gam, psi
  double precision ftaper_in, ftaper_out
  double precision tmp3, tmp4, tmp5, tmp6
  !
  if ((r .lt. andrews%rin) .or. (r .gt. andrews%rout)) then
    Andrews_dens = 0D0
    return
  end if
  !
  Md  = andrews%Md
  rc  = andrews%rc
  hc  = andrews%hc
  gam = andrews%gam
  psi = andrews%psi
  !
  tmp3 = exp(-(andrews%rin/ andrews%rc)**(2D0-gam))
  tmp4 = exp(-(andrews%rout/andrews%rc)**(2D0-gam))
  !
  sigma_c = (2D0 - gam) * Md / (phy_2Pi * rc**2) / (tmp3 - tmp4)
  !
  rrc = r / rc
  rlog = log(rrc)
  tmp1 = exp(-gam * rlog) ! = rrc**(-gam)
  tmp2 = rrc * rrc * tmp1 ! = RRc**(2D0-gam)
  !
  ! Calculte the exponential taper
  if (r .lt. andrews%r0_in_exp) then
    ftaper_in = exp( ((r-andrews%r0_in_exp)/andrews%rs_in_exp)**andrews%p_in_exp ) &
        * andrews%f_in_exp
  else
    ftaper_in = 1D0
  end if
  if (r .gt. andrews%r0_out_exp) then
    ftaper_out = exp(((andrews%r0_out_exp-r)/andrews%rs_out_exp)**andrews%p_out_exp) &
        * andrews%f_out_exp
  else
    ftaper_out = 1D0
  end if
  ! Surf mass density in Msun/AU2
  sigma = sigma_c * tmp1 * exp(-tmp2) * (ftaper_in * ftaper_out)
  !
  h = hc * exp(psi * rlog) ! rrc**psi
  !
  ! Simulate an inner or outer bump (or valley)
  if (r .lt. andrews%r0_in_change) then
    h = h * andrews%f_in_change
  else if (r .gt. andrews%r0_out_change) then
    h = h * andrews%f_out_change
  end if
  !
  tmp1 = z / h
  tmp2 = tmp1 * tmp1 * 0.5D0
  !
  if (andrews%useNumDens) then
    Andrews_dens = sigma / (phy_sqrt2Pi * h) * exp(-tmp2) * &
      phy_Msun_CGS / ((phy_AU2cm)**3 * andrews%particlemass)
  else
    Andrews_dens = sigma / (phy_sqrt2Pi * h) * exp(-tmp2) * &
      phy_Msun_CGS / ((phy_AU2cm)**3)
  end if
end function Andrews_dens


subroutine load_data_from_RADMC
  integer nx, ny, i, j
  character(len=256) pathname
  character(len=32) :: commentstr
  commentstr = ''
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


recursive subroutine delete_tree(c)
  type(type_cell), pointer, intent(inout) :: c
  integer i
  if (c%nChildren .eq. 0) then
    return
  end if
  if (.not. allocated(c%children)) then
    return
  end if
  !
  do i=1, c%nChildren
    !
    call deallocate_when_not_using(c%children(i)%p)
    !
    call delete_tree(c%children(i)%p)
    !
    nullify(c%children(i)%p%parent)
    !deallocate(c%children(i)%p)
    nullify(c%children(i)%p)
  end do
  !
  deallocate(c%children)
  !
end subroutine delete_tree



subroutine deallocate_when_not_using(c)
  type(type_cell), pointer, intent(inout) :: c
  integer stat
  if (.not. associated(c)) then
    return
  end if
  if (associated(c%par)) then
    deallocate(c%par, c%h_c_rates, c%abundances)
    deallocate(c%col_den_toISM, c%col_den_toStar)
  end if
  !
  ! Ignore any deallocation error
  if (associated(c%inner)) then
    deallocate(c%inner%idx,  stat=stat)
    deallocate(c%inner, stat=stat)
  end if
  if (associated(c%outer)) then
    deallocate(c%outer%idx,  stat=stat)
    deallocate(c%outer, stat=stat)
  end if
  if (associated(c%below)) then
    deallocate(c%below%idx,  stat=stat)
    deallocate(c%below, stat=stat)
  end if
  if (associated(c%above)) then
    deallocate(c%above%idx,  stat=stat)
    deallocate(c%above, stat=stat)
  end if
  if (associated(c%around)) then
    deallocate(c%around%idx,  stat=stat)
    deallocate(c%around, stat=stat)
  end if
  !
  if (allocated(c%optical)) then
    if (allocated(c%optical%X)) then
      deallocate(c%optical%X, &
        c%optical%ext_tot, c%optical%flux, c%optical%phc, c%optical%dir_wei, stat=stat)
    end if
    deallocate(c%optical, stat=stat)
  end if
  !
  if (allocated(c%focc)) then
    if (allocated(c%focc%vals)) then
      deallocate(c%focc%vals, stat=stat)
    end if
    deallocate(c%focc)
  end if
  if (allocated(c%cont_lut)) then
    if (allocated(c%cont_lut%J)) then
      deallocate(c%cont_lut%J, stat=stat)
    end if
    deallocate(c%cont_lut)
  end if
end subroutine deallocate_when_not_using


end module grid
