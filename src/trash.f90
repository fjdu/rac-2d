subroutine sub_divide_another(c)
  type(type_cell), target :: c
  double precision :: xslope, yslope, delx, dely, xmid, ymid
  double precision :: diffx = 1D-6, diffy = 1D-6
  double precision v0, v1, v2
  v0 = get_density_analytic(c%xmin, c%ymin)
  v1 = get_density_analytic(c%xmin + diffx, c%ymin)
  v2 = get_density_analytic(c%xmin, c%ymin + diffy)
  xslope = abs((v1 - v0) / diffx)
  yslope = abs((v2 - v0) / diffy)
  delx = v0 * 0.5D0 / xslope
  dely = v0 * 0.5D0 / yslope
  if ((delx .ge. (c%xmax - c%xmin)) .and. &
      (dely .ge. (c%ymax - c%ymin))) then
    c%nChildren = 4
    call init_children(c, c%nChildren)
    xmid = 0.5D0 * (c%xmin + c%xmax)
    ymid = 0.5D0 * (c%ymin + c%ymax)
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
  end if
  if ((delx .lt. (c%xmax - c%xmin)) .and. &
      (dely .lt. (c%ymax - c%ymin))) then
    c%nChildren = 4
    call init_children(c, c%nChildren)
    xmid = c%xmin + delx
    ymid = c%ymin + dely
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
  end if
  if ((delx .lt. (c%xmax - c%xmin)) .and. &
      (dely .ge. (c%ymax - c%ymin))) then
    c%nChildren = 2
    call init_children(c, c%nChildren)
    xmid = c%xmin + delx
    c%children(1)%p%xmin = c%xmin
    c%children(1)%p%xmax = xmid
    c%children(1)%p%ymin = c%ymin
    c%children(1)%p%ymax = c%ymax
    c%children(2)%p%xmin = xmid
    c%children(2)%p%xmax = c%xmax
    c%children(2)%p%ymin = c%ymin
    c%children(2)%p%ymax = c%ymax
  end if
  if ((delx .ge. (c%xmax - c%xmin)) .and. &
      (dely .lt. (c%ymax - c%ymin))) then
    c%nChildren = 2
    call init_children(c, c%nChildren)
    ymid = c%ymin + dely
    c%children(1)%p%xmin = c%xmin
    c%children(1)%p%xmax = c%xmax
    c%children(1)%p%ymin = c%ymin
    c%children(1)%p%ymax = ymid
    c%children(2)%p%xmin = c%xmin
    c%children(2)%p%xmax = c%xmax
    c%children(2)%p%ymin = ymid
    c%children(2)%p%ymax = c%ymax
  end if
  !write(*,*) 'XX', c%order, c%xmin, c%xmax, xmid, c%ymin, c%ymax, ymid, delx, dely
end subroutine sub_divide_another



function is_uniform_another(c)
  logical is_uniform_another
  type(type_cell), target :: c
  integer, parameter :: nx=5, ny=5
  double precision, dimension(nx, ny) :: vals
  integer i, j
  double precision dx, dy, x, y, maxv, minv
  dx = (c%xmax - c%xmin) / dble(nx-1)
  dy = (c%ymax - c%ymin) / dble(ny-1)
  x = c%xmin
  do i=1, nx
    y = c%ymin
    do j=1, ny
      vals(i, j) = get_density_analytic(x, y)
      y = y + dy
    end do
    x = x + dx
  end do
  maxv = maxval(vals)
  minv = minval(vals)
  if (maxv .gt. min_val_considered) then
    if (maxv / (minv + tiny(0D0)) .gt. max_ratio_to_be_uniform) then
      is_uniform_another = .false.
    else
      is_uniform_another = .true.
    end if
  else
    is_uniform_another = .true.
  end if
end function is_uniform_another


function is_uniform_analytic(c)
  logical is_uniform_analytic
  type(type_cell), target :: c
  integer nx, ny
  integer i, j
  double precision val
  double precision dx_0, dy_0, dx, dy, x, y
  double precision :: dx_ratio=2D0, dy_ratio=2D0
  double precision xmid, ymid, max_ratio_to_be_uniform_here
  associate(d    => refinement_data, &
            n    => refinement_data%n_idx_incell, &
            idx  => refinement_data%idx_incell, &
            maxv => refinement_data%max_val, &
            minv => refinement_data%min_val, &
            avev => refinement_data%ave_val, &
            xyw  => refinement_data%xy_weighted)
    n = -1
    idx = -1
    !
    dx_0 = max(c%xmin * dx_0_min_frac, (c%xmax - c%xmin) * dx_0_frac)
    dy_0 = (c%ymax - c%ymin) * dy_0_frac
    nx = ceiling(log( &
           (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
    ny = ceiling(log( &
           (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
    !
    is_uniform_analytic = .true.
    !
    ymid = c%ymax
    x = c%xmin
    dx = dx_0
    do i=1, nx
      y = c%ymin
      dy = dy_0
      maxv = 0D0
      minv = huge(0D0)
      do j=1, ny
        val = get_density_analytic(x, y)
        maxv = max(maxv, val)
        minv = min(minv, val)
        max_ratio_to_be_uniform_here = max_ratio_to_be_uniform * log10(3D16/maxv)
        if ((maxv .gt. min_val_considered) .and. &
            (maxv / (minv + tiny(0D0)) .gt. max_ratio_to_be_uniform)) then
          ymid = min(ymid, y)
          is_uniform_analytic = .false.
        end if
        y = y + dy
        dy = dy * dy_ratio
      end do
      x = x + dx
      dx = dx * dx_ratio
    end do
    !
    xmid = c%xmax
    y = c%ymin
    dy = dy_0
    do j=1, ny
      x = c%xmin
      dx = dx_0
      maxv = 0D0
      minv = huge(0D0)
      do i=1, nx
        val = get_density_analytic(x, y)
        maxv = max(maxv, val)
        minv = min(minv, val)
        max_ratio_to_be_uniform_here = max_ratio_to_be_uniform * log10(3D15/maxv)
        if ((maxv .gt. min_val_considered) .and. &
            (maxv / (minv + tiny(0D0)) .gt. max_ratio_to_be_uniform)) then
          xmid = min(xmid, x)
          is_uniform_analytic = .false.
        end if
        x = x + dx
        dx = dx * dx_ratio
      end do
      y = y + dy
      dy = dy * dy_ratio
    end do
    !
    if (xmid .ge. c%xmax) then
      xmid = 0.5D0 * (c%xmin + c%xmax)
    end if
    if (ymid .ge. c%ymax) then
      ymid = 0.5D0 * (c%ymin + c%ymax)
    end if
    xyw(1) = xmid
    xyw(2) = ymid
    !write(*,*) c%order, c%xmin, c%xmax, xyw(1), c%ymin, c%ymax, xyw(2), maxv, minv, avev, x, y, dy
  end associate
end function is_uniform_analytic


subroutine find_mid_4(c)
  type(type_cell), target :: c
  integer, parameter :: n = 4
  double precision, dimension(n) :: score_s, xmid_s, ymid_s
  integer i, iwin
  associate(d    => refinement_data, &
            score => refinement_data%score, &
            xmid => refinement_data%xy_weighted(1), &
            ymid => refinement_data%xy_weighted(2))
    score = 0D0
    call find_mid_LL(c, score_s(1), xmid_s(1), ymid_s(1))
    call find_mid_LR(c, score_s(2), xmid_s(2), ymid_s(2))
    call find_mid_UR(c, score_s(3), xmid_s(3), ymid_s(3))
    call find_mid_UL(c, score_s(4), xmid_s(4), ymid_s(4))
    do i=1, n
      if (score_s(i) .gt. score) then
        score = score_s(i)
        xmid = xmid_s(i)
        ymid = ymid_s(i)
        iwin = i
      end if
    end do
    !write(*,*) c%xmin, c%xmax, c%ymin, c%ymax, xmid, ymid, get_density_analytic(xmid, ymid), &
    !    get_density_analytic(c%xmin, c%ymin)/get_density_analytic(xmid, ymid)
  end associate
end subroutine find_mid_4


subroutine find_mid_LL(c, score, xmid, ymid)
  type(type_cell), target :: c
  double precision, intent(out) :: score, xmid, ymid
  integer nx, ny
  integer i, j
  double precision dx_0, dy_0, dx, dy, del_x, del_y, x, y
  double precision :: dx_ratio=1.5D0, dy_ratio=1.5D0
  double precision area_this, aspect_this, score_this
  logical is_uniform_this
  dx_0 = max((c%xmax - c%xmin) * 1D-4, very_small_len)
  dy_0 = max((c%ymax - c%ymin) * 1D-4, very_small_len)
  nx = ceiling(log( &
         (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  !
  score = 0D0
  xmid = c%xmin
  ymid = c%ymin
  if ((nx .le. 1) .and. (ny .le. 1)) then
    return
  end if
  !
  if ((c%xmin + dx_0) .gt. c%xmax) then
    write(*,*) 'Error: dx might be too large!'
    write(*,*) dx_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  if ((c%ymin + dy_0) .gt. c%ymax) then
    write(*,*) 'Error: dy might be too large!'
    write(*,*) dy_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  x = c%xmin
  dx = dx_0
  do i=1, nx
    x = x + dx
    dx = dx * dx_ratio
    if (x .gt. c%xmax) then
      exit
    end if
    y = c%ymin
    dy = dy_0
    do j=1, ny
      y = y + dy
      dy = dy * dy_ratio
      if (y .gt. c%ymax) then
        exit
      end if
      del_x = x - c%xmin
      del_y = y - c%ymin
      score_this  = min(del_x, 5D0*del_y)
      if (score_this .gt. score) then
        is_uniform_this = test_uniformity_based_Av(c%xmin, x, c%ymin, y)
        if (is_uniform_this) then
          score = score_this
          xmid = x
          ymid = y
        end if
      end if
    end do
  end do
  if ((c%xmax - xmid) .le. (xmid - c%xmin)*0.1D0) then
    xmid = c%xmax
  end if
  if ((c%ymax - ymid) .le. (ymid - c%ymin)*0.1D0) then
    ymid = c%ymax
  end if
end subroutine find_mid_LL


subroutine find_mid_LR(c, score, xmid, ymid)
  type(type_cell), target :: c
  double precision, intent(out) :: score, xmid, ymid
  integer nx, ny
  integer i, j
  double precision dx_0, dy_0, dx, dy, del_x, del_y, x, y
  double precision :: dx_ratio=1.5D0, dy_ratio=1.5D0
  double precision area_this, aspect_this, score_this
  logical is_uniform_this
  dx_0 = max((c%xmax - c%xmin) * 1D-4, very_small_len)
  dy_0 = max((c%ymax - c%ymin) * 1D-4, very_small_len)
  nx = ceiling(log( &
         (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  !
  score = 0D0
  xmid = c%xmax
  ymid = c%ymin
  if ((nx .le. 1) .and. (ny .le. 1)) then
    return
  end if
  !
  if ((c%xmax - dx_0) .lt. c%xmin) then
    write(*,*) 'Error: dx might be too large!'
    write(*,*) dx_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  if ((c%ymin + dy_0) .gt. c%ymax) then
    write(*,*) 'Error: dy might be too large!'
    write(*,*) dy_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  x = c%xmax
  dx = dx_0
  do i=1, nx
    x = x - dx
    dx = dx * dx_ratio
    if (x .lt. c%xmin) then
      exit
    end if
    y = c%ymin
    dy = dy_0
    do j=1, ny
      y = y + dy
      dy = dy * dy_ratio
      if (y .gt. c%ymax) then
        exit
      end if
      del_x = c%xmax - x
      del_y = y - c%ymin
      score_this  = min(del_x, 5D0*del_y)
      if (score_this .gt. score) then
        is_uniform_this = test_uniformity_based_Av(x, c%xmax, c%ymin, y)
        if (is_uniform_this) then
          score = score_this
          xmid = x
          ymid = y
        end if
      end if
    end do
  end do
  if ((xmid - c%xmin) .le. (c%xmax - xmid)*0.1D0) then
    xmid = c%xmin
  end if
  if ((c%ymax - ymid) .le. (ymid - c%ymin)*0.1D0) then
    ymid = c%ymax
  end if
end subroutine find_mid_LR


subroutine find_mid_UR(c, score, xmid, ymid)
  type(type_cell), target :: c
  double precision, intent(out) :: score, xmid, ymid
  integer nx, ny
  integer i, j
  double precision dx_0, dy_0, dx, dy, del_x, del_y, x, y
  double precision :: dx_ratio=1.5D0, dy_ratio=1.5D0
  double precision area_this, aspect_this, score_this
  logical is_uniform_this
  dx_0 = max((c%xmax - c%xmin) * 1D-4, very_small_len)
  dy_0 = max((c%ymax - c%ymin) * 1D-4, very_small_len)
  nx = ceiling(log( &
         (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  !
  score = 0D0
  xmid = c%xmax
  ymid = c%ymax
  if ((nx .le. 1) .and. (ny .le. 1)) then
    return
  end if
  !
  if ((c%xmax - dx_0) .lt. c%xmin) then
    write(*,*) 'Error: dx might be too large!'
    write(*,*) dx_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  if ((c%ymax - dy_0) .lt. c%ymin) then
    write(*,*) 'Error: dy might be too large!'
    write(*,*) dy_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  x = c%xmax
  dx = dx_0
  do i=1, nx
    x = x - dx
    dx = dx * dx_ratio
    if (x .lt. c%xmin) then
      exit
    end if
    y = c%ymax
    dy = dy_0
    do j=1, ny
      y = y - dy
      dy = dy * dy_ratio
      if (y .lt. c%ymin) then
        exit
      end if
      del_x = c%xmax - x
      del_y = c%ymax - y
      score_this  = min(del_x, 5D0*del_y)
      if (score_this .gt. score) then
        is_uniform_this = test_uniformity_based_Av(x, c%xmax, y, c%ymax)
        if (is_uniform_this) then
          score = score_this
          xmid = x
          ymid = y
        end if
      end if
    end do
  end do
  if ((xmid - c%xmin) .le. (c%xmax - xmid) * 0.1D0) then
    xmid = c%xmin
  end if
  if ((ymid - c%ymin) .le. (c%ymax - ymid) * 0.1D0) then
    ymid = c%ymin
  end if
end subroutine find_mid_UR


subroutine find_mid_UL(c, score, xmid, ymid)
  type(type_cell), target :: c
  double precision, intent(out) :: score, xmid, ymid
  integer nx, ny
  integer i, j
  double precision dx_0, dy_0, dx, dy, del_x, del_y, x, y
  double precision :: dx_ratio=1.5D0, dy_ratio=1.5D0
  double precision area_this, aspect_this, score_this
  logical is_uniform_this
  dx_0 = max((c%xmax - c%xmin) * 1D-4, very_small_len)
  dy_0 = max((c%ymax - c%ymin) * 1D-4, very_small_len)
  nx = ceiling(log( &
         (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
  ny = ceiling(log( &
         (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
  !
  score = 0D0
  xmid = c%xmin
  ymid = c%ymax
  if ((nx .le. 1) .and. (ny .le. 1)) then
    return
  end if
  !
  if ((c%xmin + dx_0) .gt. c%xmax) then
    write(*,*) 'Error: dx might be too large!'
    write(*,*) dx_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  if ((c%ymax - dy_0) .lt. c%ymin) then
    write(*,*) 'Error: dy might be too large!'
    write(*,*) dy_0, c%xmin, c%xmax, c%ymin, c%ymax
    return
  end if
  x = c%xmin
  dx = dx_0
  do i=1, nx
    x = x + dx
    dx = dx * dx_ratio
    if (x .gt. c%xmax) then
      exit
    end if
    y = c%ymax
    dy = dy_0
    do j=1, ny
      y = y - dy
      dy = dy * dy_ratio
      if (y .lt. c%ymin) then
        exit
      end if
      del_x = x - c%xmin
      del_y = c%ymax - y
      score_this  = min(del_x, 5D0*del_y)
      if (score_this .gt. score) then
        is_uniform_this = test_uniformity_based_Av(c%xmin, x, y, c%ymax)
        if (is_uniform_this) then
          score = score_this
          xmid = x
          ymid = y
        end if
      end if
    end do
  end do
  if ((c%xmax - xmid) .le. (xmid - c%xmin) * 0.1D0) then
    xmid = c%xmax
  end if
  if ((ymid - c%ymin) .le. (c%ymax - ymid) * 0.1D0) then
    ymid = c%ymin
  end if
end subroutine find_mid_UL


subroutine find_mid_analytic(c)
  type(type_cell), target :: c
  integer nx, ny
  integer i, j
  double precision dx_0, dy_0, dx, dy, del_x, del_y, x, y
  double precision :: dx_ratio=1.5D0, dy_ratio=1.5D0
  !double precision Av_from_top
  double precision area_this, aspect_this, score_this
  logical is_uniform_this
  associate(d    => refinement_data, &
            area => refinement_data%area, &
            aspect => refinement_data%aspect, &
            score => refinement_data%score, &
            xmid => refinement_data%xy_weighted(1), &
            ymid => refinement_data%xy_weighted(2))
    dx_0 = max((c%xmax - c%xmin) * 1D-4, very_small_len)
    dy_0 = max((c%ymax - c%ymin) * 1D-4, very_small_len)
    nx = ceiling(log( &
           (c%xmax-c%xmin)/dx_0 * (dx_ratio - 1D0) + 1D0) / log(dx_ratio))
    ny = ceiling(log( &
           (c%ymax-c%ymin)/dy_0 * (dy_ratio - 1D0) + 1D0) / log(dy_ratio))
    !
    area = 0D0
    aspect = huge(0D0)
    score = 0D0
    xmid = c%xmax
    ymid = c%ymax
    if ((nx .le. 1) .and. (ny .le. 1)) then
      return
    end if
    !
    if ((c%xmin + dx_0) .gt. c%xmax) then
      write(*,*) 'Error: dx might be too large!'
      write(*,*) d%iChild, dx_0, c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    if ((c%ymin + dy_0) .gt. c%ymax) then
      write(*,*) 'Error: dy might be too large!'
      write(*,*) d%iChild, dy_0, c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    !
    !Av_from_top = get_int_val_analytic_along_y((c%xmin+c%xmax)*0.5D0, c%ymax, root%ymax) &
    !  * phy_AU2cm * colDen2Av_coeff
    !
    x = c%xmin
    dx = dx_0
    do i=1, nx
      x = x + dx
      dx = dx * dx_ratio
      if (x .gt. c%xmax) then
        exit
      end if
      y = c%ymin
      dy = dy_0
      do j=1, ny
        y = y + dy
        dy = dy * dy_ratio
        if (y .gt. c%ymax) then
          exit
        end if
        del_x = x - c%xmin
        del_y = y - c%ymin
        area_this   = del_x * del_y
        aspect_this = max(del_x, del_y) / (min(del_x, del_y) + tiny(0D0))
        score_this  = min(del_x, 5D0*del_y)
        if (score_this .gt. score) then
          is_uniform_this = test_uniformity_based_Av(c%xmin, x, c%ymin, y)
          if (is_uniform_this) then
            !area = area_this
            !aspect = aspect_this
            score = score_this
            xmid = x
            ymid = y
          end if
        end if
      end do
    end do
  end associate
end subroutine find_mid_analytic


function test_uniformity_based_Av(xmin, xmax, ymin, ymax)
  logical test_uniformity_based_Av
  double precision, intent(in) :: xmin, xmax, ymin, ymax
  integer :: nx = 20, ny = 20
  integer n_violate
  double precision x, y, dx, dy
  double precision ave_val, ave_this, int_val, del_span
  double precision Av_from_top, Av_abs, Av_diff
  double precision, parameter :: max_Av_error_allowed = 0.1D0
  real, parameter :: err_tol = 0.1
  integer i
  test_uniformity_based_Av = .true.
  dx = (xmax - xmin) / dble(nx-1)
  dy = (ymax - ymin) / dble(ny-1)
  if (((xmax-xmin) .le. very_small_len) .and. ((ymax-ymin) .le. very_small_len)) then
    return
  end if
  ave_val = get_ave_val_analytic(xmin, xmax, ymin, ymax)
  if (ave_val .le. min_val_considered) then
    return
  end if
  !
  n_violate = 0
  x = xmin
  do i=1, nx
    int_val = get_int_val_analytic_along_y(x, ymin, ymax, del_span)
    Av_abs = int_val * phy_AU2cm * colDen2Av_coeff
    Av_diff = abs(int_val - ave_val * del_span) * phy_AU2cm * colDen2Av_coeff
    Av_from_top = get_int_val_analytic_along_y(x, ymin, root%ymax) * phy_AU2cm * colDen2Av_coeff
    if (Av_from_top .gt. Av_opaque) then
      ave_this = int_val / del_span
      if (max(ave_val, ave_this) / (min(ave_val, ave_this) + tiny(0D0)) .gt. max_ratio_to_be_uniform) then
        n_violate = n_violate + 1
      end if
    else if (Av_diff .gt. max_Av_error_allowed) then
      n_violate = n_violate + 1
    end if
    x = x + dx
  end do
  if (n_violate .gt. floor(real(nx) * err_tol)) then
    test_uniformity_based_Av = .false.
    return
  end if
  !
  n_violate = 0
  y = ymin
  Av_from_top = get_int_val_analytic_along_y(0.5D0*(xmin+xmax), ymin, root%ymax) * &
      phy_AU2cm * colDen2Av_coeff
  do i=1, ny
    int_val = get_int_val_analytic_along_x(y, xmin, xmax, del_span)
    Av_abs = int_val * phy_AU2cm * colDen2Av_coeff
    Av_diff = abs(int_val - ave_val * del_span) * phy_AU2cm * colDen2Av_coeff
    if (Av_from_top .gt. Av_opaque) then
      ave_this = int_val / del_span
      if (max(ave_val, ave_this) / (min(ave_val, ave_this) + tiny(0D0)) .gt. max_ratio_to_be_uniform) then
        n_violate = n_violate + 1
      end if
    else if (Av_diff .gt. max_Av_error_allowed) then
      n_violate = n_violate + 1
    end if
    y = y + dy
  end do
  if (n_violate .gt. floor(real(ny) * err_tol)) then
    test_uniformity_based_Av = .false.
    return
  end if
end function test_uniformity_based_Av


subroutine merge_child(c)
  type(type_cell), target :: c
  double precision xmid, ymid
  logical :: flag_merged = .false.
  if (c%nChildren .eq. 2) then
    if (test_uniformity_simple(c%xmin, c%xmax, c%ymin, c%ymax)) then
      deallocate(c%children(1)%p, c%children(2)%p)
      c%nChildren = 0
      return
    end if
  else if (c%nChildren .eq. 4) then
    xmid = c%children(1)%p%xmax
    ymid = c%children(1)%p%ymax
    if (test_uniformity_simple(c%xmin, c%xmax, c%ymin, c%ymax)) then
      deallocate(c%children(1)%p, c%children(2)%p, c%children(3)%p, c%children(4)%p)
      c%nChildren = 0
      return
    end if
    if (test_uniformity_simple(c%xmin, c%xmax, c%ymin, ymid)) then
      deallocate(c%children(1)%p, c%children(2)%p)
      c%nChildren = 3
      !!!!!!!!!!!!
      flag_merged = .true.
    end if
    if (test_uniformity_simple(c%xmin, c%xmax, ymid, c%ymax)) then
      flag_merged = .true.
    end if
    if (.not. flag_merged) then
      if (test_uniformity_simple(c%xmin, xmid, c%ymin, c%ymax)) then
        flag_merged = .true.
      end if
      if (test_uniformity_simple(xmid, c%xmax, c%ymin, c%ymax)) then
        flag_merged = .true.
      end if
    end if
  else
    return
  end if
end subroutine merge_child


subroutine sub_divide_4(c)
  type(type_cell), target :: c
  double precision xmid, ymid, del_x_1, del_x_2, del_y_1, del_y_2
  double precision min_del_x, min_del_y
  integer i
  !
  xmid = refinement_data%xy_weighted(1)
  ymid = refinement_data%xy_weighted(2)
  del_x_1 = xmid - c%xmin
  del_x_2 = c%xmax - xmid
  del_y_1 = ymid - c%ymin
  del_y_2 = c%ymax - ymid
  !
  min_del_x = min(del_x_1, del_x_2)
  min_del_y = min(del_y_1, del_y_2)
  if ((min_del_x .LE. 10D0*very_small_len) .and. &
      (min_del_y .LE. 10D0*very_small_len)) then
    return
  else if ((min_del_x .gt. very_small_len) .and. &
           (min_del_y .le. very_small_len)) then
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
  else if ((min_del_x .le. very_small_len) .and. &
           (min_del_y .gt. very_small_len)) then
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
  else
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
  end if
end subroutine sub_divide_4


subroutine sub_divide_based_on_data(c)
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
    write(*, '(A, 6ES15.7, 2I5)') 'Find ill-defined cell! (1)', &
      c%xmin, c%xmax, xmid, c%ymin, c%ymax, ymid, c%order, c%nChildren
    return
  end if
  !
  c%nChildren = 2
  call init_children(c, c%nChildren)
  !
  aspect_1 = max(max((c%ymax - c%ymin), del_x_1) / min((c%ymax - c%ymin), del_x_1), &
                 max((c%ymax - c%ymin), del_x_2) / min((c%ymax - c%ymin), del_x_2))
  aspect_2 = max(max((c%xmax - c%xmin), del_y_1) / min((c%xmax - c%xmin), del_y_1), &
                 max((c%xmax - c%xmin), del_y_2) / min((c%xmax - c%xmin), del_y_2))
  if ((aspect_1 .LE. aspect_2) .AND. (min(del_x_1, del_x_2) .GT. very_small_len)) then
    c%children(1)%p%xmin = c%xmin
    c%children(1)%p%xmax = xmid
    c%children(1)%p%ymin = c%ymin
    c%children(1)%p%ymax = c%ymax
    c%children(2)%p%xmin = xmid
    c%children(2)%p%xmax = c%xmax
    c%children(2)%p%ymin = c%ymin
    c%children(2)%p%ymax = c%ymax
  else
    c%children(1)%p%xmin = c%xmin
    c%children(1)%p%xmax = c%xmax
    c%children(1)%p%ymin = c%ymin
    c%children(1)%p%ymax = ymid
    c%children(2)%p%xmin = c%xmin
    c%children(2)%p%xmax = c%xmax
    c%children(2)%p%ymin = ymid
    c%children(2)%p%ymax = c%ymax
  end if
end subroutine sub_divide_based_on_data


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


function get_density_from_data(x, y)
  double precision get_density_from_data
  double precision, intent(in) :: x, y
  double precision local_scale
  integer :: idx_min, n
  idx_min = minloc(abs(refinement_data%xyv(1, :) - x) + abs(refinement_data%xyv(2, :) - y), 1)
  if (y .le. x * 0.3D0) then
    get_density_from_data = refinement_data%xyv(3, idx_min)
  else
    get_density_from_data = 0D0
  end if
  !if (idx_min .eq. 1) then
  !  local_scale = min( &
  !    abs(refinement_data%xyv(1, 1) - refinement_data%xyv(1, 2)), &
  !    abs(refinement_data%xyv(1, 2) - refinement_data%xyv(1, 3))) &
  !    + min( &
  !    abs(refinement_data%xyv(2, 1) - refinement_data%xyv(2, 2)), &
  !    abs(refinement_data%xyv(2, 2) - refinement_data%xyv(2, 3)))
  !else if (idx_min .eq. refinement_data%nlen) then
  !  n = refinement_data%nlen
  !  local_scale = min( &
  !    abs(refinement_data%xyv(1, n-1) - refinement_data%xyv(1, n)), &
  !    abs(refinement_data%xyv(1, n-2) - refinement_data%xyv(1, n-1))) &
  !    + min( &
  !    abs(refinement_data%xyv(2, n-1) - refinement_data%xyv(2, n)), &
  !    abs(refinement_data%xyv(2, n-2) - refinement_data%xyv(2, n-1)))
  !else
  !  local_scale = min( &
  !    abs(refinement_data%xyv(1, idx_min + 1) - refinement_data%xyv(1, idx_min)), &
  !    abs(refinement_data%xyv(1, idx_min - 1) - refinement_data%xyv(1, idx_min))) &
  !    + min( &
  !    abs(refinement_data%xyv(2, idx_min + 1) - refinement_data%xyv(2, idx_min)), &
  !    abs(refinement_data%xyv(2, idx_min - 1) - refinement_data%xyv(2, idx_min)))
  !end if
  !local_scale = local_scale * 2D0
  !if (((refinement_data%xyv(1, idx_min) - x)**2 + (refinement_data%xyv(2, idx_min) - y)**2) &
  !    .lt. local_scale**2) then
  !  get_density_from_data = refinement_data%xyv(3, idx_min)
  !else
  !  get_density_from_data = 0D0
  !end if
end function get_density_from_data


subroutine disk_set_cell_params
  integer i, j
  do i=1, a_disk%n_columns
    a_disk%columns(i)%cells(1)%params%dz   = a_disk%columns(i)%params%dz_0
    a_disk%columns(i)%cells(1)%params%zmin = a_disk%columns(i)%params%zmin
    a_disk%columns(i)%cells(1)%params%zmax = a_disk%columns(i)%cells(1)%params%zmin + &
                                             a_disk%columns(i)%cells(1)%params%dz
    do j=2, a_disk%columns(i)%n_cells
      a_disk%columns(i)%cells(j)%params%dz   = a_disk%columns(i)%cells(j-1)%params%dz * &
                                               a_disk%columns(i)%params%dz_ratio
      a_disk%columns(i)%cells(j)%params%zmin = a_disk%columns(i)%cells(j-1)%params%zmax
      a_disk%columns(i)%cells(j)%params%zmax = &
        a_disk%columns(i)%cells(j)%params%zmin + a_disk%columns(i)%cells(j)%params%dz
    end do
    do j=1, a_disk%columns(i)%n_cells
      associate(aijp => a_disk%columns(i)%cells(j)%params)
        aijp%rmin = a_disk%columns(i)%params%rmin
        aijp%rmax = a_disk%columns(i)%params%rmax
        aijp%rcen = a_disk%columns(i)%params%rcen
        aijp%dr   = a_disk%columns(i)%params%dr
        aijp%zcen = 0.5D0 * (aijp%zmin + aijp%zmax)
        aijp%daz         = a_disk%columns(i)%params%daz
        ! Provide a simple temperature and density profile.
        ! Will not be used if a tabulated input is provided.
        aijp%Tgas = 100D0 / (1D0 + aijp%rcen) * (1D0 + aijp%zcen)
        aijp%Tdust = aijp%Tgas
        aijp%n_gas = a_disk%params%n_gas_scale * &
          (a_disk%params%r_charac/aijp%rcen)**(a_disk%params%r_gamma) * &
          exp(-(aijp%rcen/a_disk%params%r_charac)**(2D0-a_disk%params%r_gamma)) * &
          exp(-(aijp%zcen / a_disk%columns(i)%params%scale_height)**2)
        aijp%ratioDust2HnucNum = & ! n_Grain/n_H
          aijp%ratioDust2GasMass * (phy_mProton_CGS * aijp%MeanMolWeight) &
          / (4.0D0*phy_Pi/3.0D0 * (aijp%GrainRadius_CGS)**3 * &
             aijp%GrainMaterialDensity_CGS)
        aijp%dust_depletion = aijp%ratioDust2GasMass / ratioDust2GasMass_ISM
        aijp%velo_width_turb = 1D5 ! Todo
      end associate
      associate( &
        G     => phy_GravitationConst_CGS, &
        M     => a_disk%params%star_mass_in_Msun * phy_Msun_CGS, &
        r     => a_disk%columns(i)%cells(j)%params%rcen * phy_AU2cm, &
        v     => a_disk%columns(i)%cells(j)%params%velo_Kepler, &
        w     => a_disk%columns(i)%cells(j)%params%omega_Kepler, &
        dv_dr => a_disk%columns(i)%cells(j)%params%velo_gradient)
        v = sqrt(G * M / r)
        w = v / r
        dv_dr = 0.5D0 * v / r
      end associate
      !a_disk%columns(i)%cells(j)%params%omega_albedo = 0.5D0
      !a_disk%columns(i)%cells(j)%params%zeta_cosmicray_H2 = 1.36D-17
      !a_disk%columns(i)%cells(j)%params%R_H2_form_rate = 0D0 ! Will be calculated.
      !a_disk%columns(i)%cells(j)%params%GrainMaterialDensity_CGS = 2D0
      !a_disk%columns(i)%cells(j)%params%GrainRadius_CGS = 0.1D-4
      !a_disk%columns(i)%cells(j)%params%aGrainMin_CGS = 0.1D-4
      !a_disk%columns(i)%cells(j)%params%aGrainMax_CGS = 0.1D-4
      !a_disk%columns(i)%cells(j)%params%ratioDust2GasMass = 0.01D0
      !a_disk%columns(i)%cells(j)%params%MeanMolWeight = 1.4D0
      !a_disk%columns(i)%cells(j)%params%stickCoeffH = 1.0D0
    end do
  end do
end subroutine disk_set_cell_params


subroutine disk_set_column_params
  integer, parameter :: const_n_columns_guess = 256
  double precision, dimension(const_n_columns_guess) :: r_s, dr_s
  integer i, j, stat
  r_s(1) = a_disk%params%rmin
  do i=1, const_n_columns_guess
    if (r_s(i) .GT. a_disk%params%rmax) then
      exit
    end if
    dr_s(i) = min(a_disk%params%dr_factor * r_s(i), &
      get_local_doppler_kepler_scale( &
        a_disk%params%star_mass_in_Msun, &
        r_s(i), get_local_dv_microturb(), 2D0))
    r_s(i+1) = r_s(i) + dr_s(i)
  end do
  a_disk%n_columns = i - 1
  allocate(a_disk%columns(a_disk%n_columns))
  a_disk%n_cell_total = 0
  do i=1, a_disk%n_columns
    associate(aip => a_disk%columns(i)%params)
      aip%rmin = r_s(i)
      aip%rmax = r_s(i+1)
      aip%dr   = dr_s(i)
      aip%rcen = 0.5D0 * (aip%rmin + aip%rmax)
      aip%zmin = 0D0
      aip%scale_height = (aip%rcen + 0.25D0)**0.5D0 * a_disk%params%height_factor
      aip%zmax = aip%zmin + aip%scale_height
      aip%dz_0 = max(aip%zmax * a_disk%params%dz_factor, dr_s(i) * 0.1D0)
      aip%dz_ratio = a_disk%params%dz_ratio
      aip%daz = aip%dr
      !
      ! A consistent angle correction should be added.
      ! The continuum UV flux is relative to the ISRF, while the Lyman alpha flux is absolte.
      aip%UV_flux_top = 1D0 + &
        a_disk%params%UV_cont_phlumi_star_surface &
           / (4D0*phy_Pi * (aip%rcen * phy_AU2cm)**2) &
           / phy_Habing_photon_flux_CGS &
           !/ phy_UV_flux_ISM & Using this value makes the code very slow.
           * const_geometric_factor_UV
      ! No geometric dilution for Lyman alpha.
      aip%LymanAlpha_flux_top = &
        a_disk%params%Lyman_phlumi_star_surface &
           / (4D0*phy_Pi * (aip%rcen * phy_AU2cm)**2)
      aip%Xray_flux_top = &
        a_disk%params%Xray_phlumi_star_surface &
           / (4D0*phy_Pi * (aip%rcen * phy_AU2cm)**2) &
           * const_geometric_factor_Xray
      aip%cosmicray_flux_top = 1D0
      !
      a_disk%columns(i)%n_cells = ceiling( &
        log(aip%zmax / aip%dz_0 * (aip%dz_ratio - 1D0) + 1D0) &
        / log(aip%dz_ratio))
    end associate
    a_disk%columns(i)%n_cells = min(a_disk%columns(i)%n_cells, max(1,a_disk%params%ncellpercol_max))
    allocate(a_disk%columns(i)%cells(a_disk%columns(i)%n_cells))
    do j=1, a_disk%columns(i)%n_cells
      allocate(a_disk%columns(i)%cells(j)%abundances(chem_species%nSpecies), &
               a_disk%columns(i)%cells(j)%col_den(chem_idx_some_spe%nItem), &
               a_disk%columns(i)%cells(j)%col_den_acc(chem_idx_some_spe%nItem), STAT=stat)
      if (stat .NE. 0) then
        write(*,*) 'Fail to allocate memory!'
        stop
      end if
      ! Set the initial parameters for each cell.
      ! These initial values are from the config file, and are
      ! unlikely to vary with position or be be changed during the iteration.
      a_disk%columns(i)%cells(j)%params%type_cell_rz_phy_basic = cell_params_ini
    end do
    a_disk%n_cell_total = a_disk%n_cell_total + a_disk%columns(i)%n_cells
  end do
end subroutine disk_set_column_params


!type :: type_chemical_evol_parameters
!  double precision :: Tgas, Tdust, n_gas, UV_G0_factor, Xray_flux_0, Av, Ncol, &
!    LymanAlpha_flux_0, &
!    omega_albedo, zeta_cosmicray_H2, R_H2_form_rate, &
!    f_selfshielding_H2, f_selfshielding_CO, f_selfshielding_H2O, f_selfshielding_OH, &
!    GrainMaterialDensity_CGS, GrainRadius_CGS, aGrainMin_CGS, aGrainMax_CGS, &
!    ratioDust2GasMass, ratioDust2HnucNum, dust_depletion, MeanMolWeight, &
!    stickCoeffH
!end type type_chemical_evol_parameters

type :: a__cell
  type(phy_chem_rad_cell_params) :: params
  double precision, dimension(:), allocatable :: abundances
  double precision, dimension(:), allocatable :: col_den, col_den_acc
  logical :: converged = .false.
end type a__cell


type :: a__column
  integer n_cells
  type(phy_chem_rad_column_params) :: params
  type(a__cell), dimension(:), allocatable :: cells
end type a__column


subroutine load_dust_dens_T_from_RADMC
  ! This subroutine is supposed to be temporary only.
  double precision, dimension(:), allocatable :: r, z, rho, T
  integer fU, ios, i, j, k, l, i0
  integer, dimension(1) :: idx
  character(len=128) str
  if (.NOT. getFileUnit(fU)) then
    write(*,*) 'No free file unit!'
    stop
  end if
  l = GetFileLen_comment_blank(combine_dir_filename(a_disk%params%input_struct_dir, &
                               a_disk%params%filename_RADMC), '!')
  allocate(r(l), z(l), rho(l), T(l))
  call openFileSequentialRead(fU, combine_dir_filename(a_disk%params%input_struct_dir, &
                                  a_disk%params%filename_RADMC), 99999)
  i = 0
  do
    read(fU, '(A)', IOSTAT=ios) str
    if (ios .LT. 0) then
      exit
    end if
    if ((str(1:1) .NE. '!') .AND. (len_trim(str) .GT. 0)) then
       i = i + 1
      read(str, '(2X, 3(ES11.2, 2X), 4X, F7.1)') r(i), z(i), rho(i), T(i)
    end if
  end do
  close(fU)
  r = r / phy_AU2cm
  z = z / phy_AU2cm
  !write(*,*) 'minval(r), maxval(r), minval(z), maxval(z)', &
  !            minval(r), maxval(r), minval(z), maxval(z)
  if (a_disk%params%rescale_with_input) then
    ! Scale the vertical height to the RADMC file.
    do i=1, a_disk%n_columns
      a_disk%params%height_rescale_factor = 0D0
      do k=1, l
        if ((r(k) .ge. a_disk%columns(i)%params%rmin) .and. &
            (r(k) .le. a_disk%columns(i)%params%rmax)) then
          a_disk%params%height_rescale_factor = &
            max(a_disk%params%height_rescale_factor, &
                z(k) / a_disk%columns(i)%params%zmax)
        end if
      end do
      call rescale_column_heights_this_do(i)
    end do
  end if
  do i=1, a_disk%n_columns
    do j=1, a_disk%columns(i)%n_cells
      i0 = 0
      associate(Tdust  => a_disk%columns(i)%cells(j)%params%Tdust, &
                Tgas   => a_disk%columns(i)%cells(j)%params%Tgas, &
                n_gas  => a_disk%columns(i)%cells(j)%params%n_gas)
        Tdust = 0D0
        n_gas = 0D0
        do k=1, l
          if ( &
              (r(k) .GE. a_disk%columns(i)%cells(j)%params%rmin) .AND. &
              (r(k) .LE. a_disk%columns(i)%cells(j)%params%rmax) .AND. &
              (z(k) .GE. a_disk%columns(i)%cells(j)%params%zmin) .AND. &
              (z(k) .LE. a_disk%columns(i)%cells(j)%params%zmax)) then
            i0 = i0 + 1
            Tdust = Tdust + T(k)
            n_gas = n_gas + rho(k)
          end if
        end do
        if (i0 .GT. 0) then
          Tdust = Tdust / dble(i0)
          n_gas = n_gas / dble(i0)
        else
          idx = minloc((r - a_disk%columns(i)%cells(j)%params%rcen)**2 + &
                       (z - a_disk%columns(i)%cells(j)%params%zcen)**2)
          Tdust = T(idx(1))
          n_gas = rho(idx(1))
        end if
        Tgas = Tdust
        associate(d2g_mass_ratio => a_disk%columns(i)%cells(j)%params%ratioDust2GasMass, &
                  mean_gas_mass => a_disk%columns(i)%cells(j)%params%MeanMolWeight * phy_mProton_CGS)
          !n_gas = n_gas / d2g_mass_ratio / mean_gas_mass
          n_gas = n_gas / 0.01D0 / mean_gas_mass
        end associate
      end associate
    end do
  end do
end subroutine load_dust_dens_T_from_RADMC


subroutine rescale_column_heights_do
  integer i, j
  associate(f => a_disk%params%height_rescale_factor)
    do i=1, a_disk%n_columns
      associate(aip => a_disk%columns(i)%params)
        aip%zmin = aip%zmin * f
        aip%zmax = aip%zmax * f
        aip%dz_0 = aip%dz_0 * f
      end associate
      do j=1, a_disk%columns(i)%n_cells
        associate(aijp => a_disk%columns(i)%cells(j)%params)
          aijp%zmin = aijp%zmin * f
          aijp%zmax = aijp%zmax * f
          aijp%zcen = aijp%zcen * f
          aijp%dz   = aijp%dz   * f
        end associate
      end do
    end do
  end associate
end subroutine rescale_column_heights_do


subroutine rescale_column_heights_this_do(icol)
  integer j, icol
  associate(f => a_disk%params%height_rescale_factor)
    associate(aip => a_disk%columns(icol)%params)
      aip%zmin = aip%zmin * f
      aip%zmax = aip%zmax * f
      aip%dz_0 = aip%dz_0 * f
    end associate
    do j=1, a_disk%columns(icol)%n_cells
      associate(aijp => a_disk%columns(icol)%cells(j)%params)
        aijp%zmin = aijp%zmin * f
        aijp%zmax = aijp%zmax * f
        aijp%zcen = aijp%zcen * f
        aijp%dz   = aijp%dz   * f
      end associate
    end do
    if (FileUnitOpened(a_book_keeping%fU)) then
      write(a_book_keeping%fU, '(A, 2X, I4, 2X, A, F9.4)') &
        '! Rescale column', icol, 'by a factor of ', f
      flush(a_book_keeping%fU)
    end if
  end associate
end subroutine rescale_column_heights_this_do


type :: phy_chem_rad_column_params
  double precision :: &
    UV_flux_top, &
    LymanAlpha_flux_top, &
    Xray_flux_top, &
    cosmicray_flux_top, &
    scale_height
  double precision dv_microturb
  double precision rmin, rmax, rcen, dr, zmin, zmax, dz_0, dz_ratio, daz
end type phy_chem_rad_column_params


type, extends(type_cell_rz_phy_basic) :: phy_chem_rad_cell_params
  type(type_heating_cooling_rates_list) h_c_rates
end type phy_chem_rad_cell_params


!subroutine disk_calc_column_densities
!  integer i, j, i0
!  do i=1, cell_leaves%nlen
!    associate(p => cell_leaves%list(i)%p, &
!              dz => cell_leaves%list(i)%p%par%dz * phy_AU2cm)
!      p%col_den   = p%abundances(chem_idx_some_spe%idx) * p%par%n_gas * dz
!      p%par%dNcol = p%par%n_gas * dz
!      p%col_den_acc = 0D0
!      p%par%Ncol = 0D0
!    end associate
!  end do
!  cell_leaves%stat = 0
!  do i=1, cell_leaves%nlen
!    associate(p => cell_leaves%list(i)%p)
!      if (cell_leaves%stat(i) .eq. 1) then
!        cycle
!      end if
!      do j=1, p%above%n
!        i0 = p%above%idx(j)
!        call calc_col_den_recursive(i0)
!        p%col_den_acc = p%col_den_acc + &
!          (cell_leaves%list(i0)%p%col_den_acc + &
!           cell_leaves%list(i0)%p%col_den) * p%above%fra(j)
!        p%par%Ncol = p%par%Ncol + &
!          (cell_leaves%list(i0)%p%par%Ncol + &
!           cell_leaves%list(i0)%p%par%dNcol) * p%above%fra(j)
!      end do
!      cell_leaves%stat(i) = 1
!    end associate
!  end do
!end subroutine disk_calc_column_densities
!
!
!recursive subroutine calc_col_den_recursive(ic)
!  integer, intent(in) :: ic
!  integer i, i0
!  if (cell_leaves%stat(ic) .eq. 1) then
!    return
!  else
!    cell_leaves%stat(ic) = 1
!    do i=1, cell_leaves%list(ic)%p%above%n
!      i0 = cell_leaves%list(ic)%p%above%idx(i)
!      call calc_col_den_recursive(i0)
!      cell_leaves%list(ic)%p%col_den_acc  = cell_leaves%list(ic)%p%col_den_acc + &
!        (cell_leaves%list(i0)%p%col_den_acc + cell_leaves%list(i0)%p%col_den) * &
!        cell_leaves%list(ic)%p%above%fra(i)
!      cell_leaves%list(ic)%p%par%Ncol  = cell_leaves%list(ic)%p%par%Ncol + &
!        (cell_leaves%list(i0)%p%par%Ncol + cell_leaves%list(i0)%p%par%dNcol) * &
!        cell_leaves%list(ic)%p%above%fra(i)
!    end do
!  end if
!end subroutine calc_col_den_recursive


subroutine disk_set_cell_coupled_params
  use load_Visser_CO_selfshielding
  integer i
  !
  !call disk_calc_column_densities
  !
  do i=1, cell_leaves%nlen
    associate( &
      p        => cell_leaves%list(i)%p%par, &
      Ncol     => cell_leaves%list(i)%p%par%Ncol, &
      dNcol    => cell_leaves%list(i)%p%par%dNcol, &
      Ncol_H2  => cell_leaves%list(i)%p%col_den_acc(chem_idx_some_spe%iiH2), &
      dcol_H2  => cell_leaves%list(i)%p%col_den(chem_idx_some_spe%iiH2), &
      Ncol_H   => cell_leaves%list(i)%p%col_den_acc(chem_idx_some_spe%iiHI), &
      dcol_H   => cell_leaves%list(i)%p%col_den(chem_idx_some_spe%iiHI), &
      Ncol_H2O => cell_leaves%list(i)%p%col_den_acc(chem_idx_some_spe%iiH2O), &
      dcol_H2O => cell_leaves%list(i)%p%col_den(chem_idx_some_spe%iiH2O), &
      Ncol_OH  => cell_leaves%list(i)%p%col_den_acc(chem_idx_some_spe%iiOH), &
      dcol_OH  => cell_leaves%list(i)%p%col_den(chem_idx_some_spe%iiOH), &
      Ncol_CO  => cell_leaves%list(i)%p%col_den_acc(chem_idx_some_spe%iiCO), &
      dcol_CO  => cell_leaves%list(i)%p%col_den(chem_idx_some_spe%iiCO) &
      )
      ! Kwok eq 10.20
      p%Av = 1.086D0 * p%ratioDust2HnucNum * &
        (phy_Pi * p%GrainRadius_CGS**2) * 2D0 * &
        (Ncol + dNcol * 0.5D0)
      p%f_selfshielding_H2  = &
        min(1D0, ((Ncol_H2 + dcol_H2*0.5D0)/1D14)**(-0.75D0)) ! Tielens 2005, equation 8.39
      p%f_selfshielding_H2O = &
        min(1D0, exp(-(Ncol_H2O * const_LyAlpha_cross_H2O))) * &
        tau2beta(dcol_H2O * const_LyAlpha_cross_H2O)
      p%f_selfshielding_OH  = &
        min(1D0, exp(-(Ncol_OH * const_LyAlpha_cross_OH))) * &
        tau2beta(dcol_OH * const_LyAlpha_cross_OH)
      p%f_selfshielding_CO = get_12CO_shielding(Ncol_H2, Ncol_CO)
      write(*,'(/2ES16.6/)') Ncol_H2, Ncol_CO
    end associate
  end do
end subroutine disk_set_cell_coupled_params




!!double precision abs_sig, sca_sig, sca_g
!!allocate(a_dust_set%list(2))
!!call load_dust_optical_property('silicate.Kappa', 1)
!!do i=1, a_dust_set%list(1)%p%nlen
!!  write(*,'(I5, 4ES14.4)') i, a_dust_set%list(1)%p%wavelen(i), a_dust_set%list(1)%p%freq(i), &
!!    a_dust_set%list(1)%p%sig_abs(i), a_dust_set%list(1)%p%sig_sca(i)
!!end do
!!call get_dust_abs_sca_crosssec(1, 2D15, abs_sig, sca_sig, sca_g)
!!write(*, '(3ES14.4)') abs_sig, sca_sig, sca_g
!!
!!allocate(a_cell_dust%list(1))
!!a_cell_dust%list(1)%material = 'silicate'
!!a_cell_dust%list(1)%radius_CGS = 0.1D-4
!!a_cell_dust%list(1)%material_density_CGS = 2D0
!!call calc_dust_particle_params_basic(a_cell_dust%list(1))
!!write(*, '(ES14.4)') a_cell_dust%list(1)%particle_mass_CGS * abs_sig
!!stop

!!! Testing
type(type_spline_2D) :: sp2d
type(type_spline_1D) :: sp1d
double precision, dimension(201) ::  x, y
integer s1, s2
!
sp2d%nx = 11
sp2d%ny = 11
allocate( &
  sp2d%xi(sp2d%nx), &
  sp2d%yi(sp2d%ny), &
  sp2d%vi(sp2d%nx, sp2d%ny))
sp2d%xi = (/(dble(i)*0.5D0*(1D0+dble(i)/1D3), i=0,10)/)
sp2d%yi = sp2d%xi
call openFileSequentialWrite(10, 'tmp1', 99999)
do i=1, sp2d%nx
do j=1, sp2d%ny
  sp2d%vi(i, j) = sin(sp2d%xi(i)**2)*exp(-sp2d%yi(j)/2D0) + cos(sp2d%yi(j))
  write(10,'(3ES14.4)') sp2d%xi(i), sp2d%yi(j), sp2d%vi(i, j)
end do
end do
call spline2d_prepare(sp2d)
close(10)
!
x = (/(dble(i)*0.1D0-5D0, i=0,200)/)
y = x
call openFileSequentialWrite(10, 'tmp2', 99999)
do i=1, 201
do j=1, 201
write(10,'(3ES14.4)') x(i), y(j), spline2d_interpol(x(i), y(j), sp2d, s1, s2)
end do
end do
close(10)
!
!call spline1d_prepare(sp2d%sp1d_using)
!write(*,*) sp2d%sp1d_using%ddyi
!write(*,*) spline1d_interpol(1.5D0, sp2d%sp1d_using)
!!
!sp1d%n = 5
!allocate(sp1d%xi(5), sp1d%yi(5), sp1d%ddyi(5))
!sp1d%xi = sp2d%xi
!sp1d%yi = sp2d%vi(:, 1)
!call spline1d_prepare(sp1d)
!write(*,*) spline1d_interpol(1.5D0, sp1d)
stop
return
!!! END Testing



! Test
! Tested against matlab.
! x=[1D0,2D0,3D0,4D0,5D0]
! y=[1D0,1D1,1D2,4D1,5D0]
! xx and yy are taken from the output of this code.
! plot(x,y,'b.',xx,yy,'b',xx,spline(x,y,xx),'mo');
type(type_spline_1D) :: sp1d
type(type_barycentric_1d) :: bary1d
double precision, dimension(505) :: x, y
sp1d%n = 5
allocate(sp1d%xi(5), sp1d%yi(5), sp1d%ddyi(5))
sp1d%xi = (/1D0,2D0,3D0,4D0,5D0/)
sp1d%yi = (/1D0,1D1,1D2,4D1,5D0/)
call spline1d_prepare(sp1d)
do i=1,505
  x(i) = dble(i) * 0.01D0
  y(i) = spline1d_interpol(x(i), sp1d)
end do
write(*,'(505ES14.4)') x
write(*,'(505ES14.4)') y
!
bary1d%n = 5
bary1d%d = -1
allocate(bary1d%xi(5), bary1d%yi(5), bary1d%wi(5))
bary1d%xi = (/1D0,2D0,3D0,4D0,5D0/)
bary1d%yi = (/1D0,1D1,1D2,4D1,5D0/)
call barycentric1d_prepare(bary1d)
do i=1,505
  y(i) = barycentric1d_interpol(x(i), bary1d)
end do
write(*,'(505ES14.4)') y
stop
!!! END Test


character(len=64) :: str = ' 1.23 aaqw12Jhdsbb ksf >2.1Ad?? sfsd .23abdbsad     '
character(len=8), dimension(8) :: str_split
integer nout
call split_str_by_space(str, str_split, 8, nout)
do i=1, nout
  write(*,*) i, str_split(i)
end do
