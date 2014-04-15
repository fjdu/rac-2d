module ray_propagating

use grid
use data_struct

implicit none

integer, parameter, private :: Undef_dirtype = -100

contains

subroutine enter_the_domain(ray, cstart, cnext, found)
  ! The photon location will be updated to the entry point if it does enter
  ! the cell.
  ! The entry point is shifted a little bit along the photon propagation
  ! direction.
  ! cnext will point to the leaf cell containing the entry point.
  type(type_ray), intent(inout) :: ray
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer, intent(out) :: cnext
  logical, intent(out) :: found
  double precision length, r, z, eps
  integer dirtype
  !
  r = ray%x**2 + ray%y**2
  z = ray%z
  if (is_inside_cell_sq(r, z, cstart)) then
    found = .true.
    length = 0D0
    eps = 0D0
  else
    call calc_intersection_ray_cell(ray, cstart, &
      length, r, z, eps, found, dirtype)
  end if
  if (found) then
    call locate_photon_cell_by_tree(r, z, cstart, cnext, found)
    if (found) then
      !eps = min(eps, 1D-2*(cnext%xmax-cnext%xmin))
      ray%x = ray%x + ray%vx * (length + eps)
      ray%y = ray%y + ray%vy * (length + eps)
      ray%z = ray%z + ray%vz * (length + eps)
      ! A further check to make sure the photon is really inside the cell.
      call calc_intersection_ray_cell(ray, cnext, &
        length, r, z, eps, found, dirtype)
    end if
  end if
end subroutine enter_the_domain



subroutine enter_the_domain_mirror(ray, cstart, cnext, found)
  ! The photon location will be updated to the entry point if it does enter
  ! the cell.
  ! The entry point is shifted a little bit along the photon propagation
  ! direction.
  ! cnext will point to the leaf cell containing the entry point.
  type(type_ray), intent(inout) :: ray
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer, intent(out) :: cnext
  logical, intent(out) :: found
  double precision length, r, z, eps
  integer dirtype
  !
  r = ray%x**2 + ray%y**2
  z = ray%z
  if (is_inside_cell_mirror_sq(r, z, cstart)) then
    found = .true.
    length = 0D0
    eps = 0D0
  else
    call calc_intersection_ray_cell_mirror(ray, cstart, &
      length, r, z, eps, found, dirtype)
  end if
  if (found) then
    call locate_photon_cell_mirror(r, z, cstart, cnext, found)
    if (found) then
      ray%x = ray%x + ray%vx * (length + eps)
      ray%y = ray%y + ray%vy * (length + eps)
      ray%z = ray%z + ray%vz * (length + eps)
      ! A further check to make sure the photon is really inside the cell.
      call calc_intersection_ray_cell_mirror(ray, cnext, &
        length, r, z, eps, found, dirtype)
    end if
  end if
end subroutine enter_the_domain_mirror



subroutine locate_photon_cell_alt(r, z, c, dirtype,  cout, found)
  ! Given r and z and start from c, find out a cell containing (r,z).
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  integer, intent(in) :: dirtype
  type(type_cell), pointer, intent(out) :: cout
  logical, intent(out) :: found
  integer i
  type(type_neighbor), pointer :: neib
  if (c%using) then
    select case (dirtype)
      case (1)
        neib => c%above
      case (2)
        neib => c%below
      case (3,4)
        neib => c%inner
      case (5,6)
        neib => c%outer
      case (Undef_dirtype)
        found = .false.
        write(*, '(A)') 'In locate_photon_cell_alt:'
        write(*, '(A)') 'dirtype undefined!'
        write(*, '(A, 6ES16.9/)') 'r,z,cxy = ', r, z, c%xmin**2, c%xmax**2, c%ymin, c%ymax
        return
      case default
        write(*, '(A)') 'In locate_photon_cell_alt:'
        write(*, '(A, I4)') 'dirtype = ', dirtype
        write(*, '(A)') 'Photon still in, probably due to numerical err.'
        write(*, '(A, 6ES16.9/)') 'r,z,cxy = ', r, z, c%xmin**2, c%xmax**2, c%ymin, c%ymax
        cout => c
        found = .true.
        return
    end select
    do i=1, neib%n
      if (is_inside_cell_sq(r, z, leaves%list(neib%idx(i))%p)) then
        cout => leaves%list(neib%idx(i))%p
        found = .true.
        return
      end if
    end do
  end if
  call locate_photon_cell_by_tree(r, z, c, cout, found)
end subroutine locate_photon_cell_alt



subroutine locate_photon_cell_by_tree(r, z, c, cout, found)
  ! Given r and z and start from c, find out a cell containing (r,z).
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  type(type_cell), pointer, intent(out) :: cout
  logical, intent(out) :: found
  integer i, j
  logical flag
  integer, parameter :: NMAX = 1000000
  !
  found = .false.
  cout => c
  !
  do j=1, NMAX
    if (is_inside_cell_sq(r, z, cout)) then
      if (cout%nChildren .eq. 0) then
        found = .true.
        return
      end if
      flag = .true.
      do i=1, cout%nChildren
        if (is_inside_cell_sq(r, z, cout%children(i)%p)) then
          cout => cout%children(i)%p
          flag = .false.
          exit
        end if
      end do
      if (flag) then
        cout => null()
        return
      end if
    else
      if (associated(cout%parent)) then
        cout => cout%parent
      else
        cout => null()
        return
      end if
    end if
  end do
  cout => null()
  return
end subroutine locate_photon_cell_by_tree



subroutine locate_photon_cell_mirror(r, z, c, cout, found)
  ! Given r and z and start from c, find out a cell containing (r,z).
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  type(type_cell), pointer, intent(out) :: cout
  logical, intent(out) :: found
  integer i, j
  logical flag
  integer, parameter :: NMAX = 1000000
  !
  found = .false.
  cout => c
  !
  do j=1, NMAX
    if (is_inside_cell_mirror_sq(r, z, cout)) then
      if (cout%nChildren .eq. 0) then
        found = .true.
        return
      end if
      flag = .true.
      do i=1, cout%nChildren
        if (is_inside_cell_mirror_sq(r, z, cout%children(i)%p)) then
          cout => cout%children(i)%p
          flag = .false.
          exit
        end if
      end do
      if (flag) then
        cout => null()
        return
      end if
    else
      if (associated(cout%parent)) then
        cout => cout%parent
      else
        cout => null()
        return
      end if
    end if
  end do
  cout => null()
  return
end subroutine locate_photon_cell_mirror



pure subroutine calc_intersection_ray_cell_mirror(ray, c, length, r, z, eps, found, dirtype)
  type(type_ray), intent(in) :: ray
  type(type_cell), pointer, intent(in) :: c
  double precision, intent(out) :: length, r, z, eps
  logical, intent(out) :: found
  integer, intent(out) :: dirtype
  type(type_ray) ray2
  !
  double precision :: length1, r1, z1, eps1
  double precision :: length2, r2, z2, eps2
  logical :: found1
  logical :: found2
  integer :: dirtype1
  integer :: dirtype2
  !
  call calc_intersection_ray_cell(ray,  c, length1, r1, z1, eps1, found1, dirtype1)
  !
  ray2%x  =  ray%x
  ray2%y  =  ray%y
  ray2%z  = -ray%z
  ray2%vx =  ray%vx
  ray2%vy =  ray%vy
  ray2%vz = -ray%vz
  !
  ! Here may be made more efficient.
  call calc_intersection_ray_cell(ray2, c, length2, r2, z2, eps2, found2, dirtype2)
  !
  found = found1 .or. found2
  !
  if (found1 .and. found2) then
    if (length1 .le. length2) then
      length = length1
      r = r1
      z = z1
      eps = eps1
      dirtype = dirtype1
    else
      length = length2
      r = r2
      z = -z2
      eps = eps2
      dirtype = dirtype2
      if (dirtype .eq. 1) then
        dirtype = 2
      else if (dirtype .eq. 2) then
        dirtype = 1
      end if
    end if
  else if (found1) then
    length = length1
    r = r1
    z = z1
    eps = eps1
    dirtype = dirtype1
  else if (found2) then
    length = length2
    r = r2
    z = -z2
    eps = eps2
    dirtype = dirtype2
    if (dirtype .eq. 1) then
      dirtype = 2
    else if (dirtype .eq. 2) then
      dirtype = 1
    end if
  else
    dirtype = Undef_dirtype
  end if
  !
end subroutine calc_intersection_ray_cell_mirror


pure subroutine calc_intersection_ray_cell(ray, c, length, rsq, z, eps, found, dirtype)
  type(type_ray), intent(in) :: ray
  type(type_cell), intent(in) :: c
  double precision, intent(out) :: length, rsq, z, eps
  logical, intent(out) :: found
  integer, intent(out) :: dirtype
  double precision rr, zz, A, B, C1, C2, D1, D2
  double precision t1, t2
  double precision, dimension(6) :: L
  double precision, parameter :: eps_ratio = 1D-6
  integer idx, i
  logical flag_inside_cell
  double precision, parameter :: FL = -1D0 ! False length
  double precision, parameter :: MinLen = 1D-30
  double precision, parameter :: MinVz  = 1D-20
  double precision, parameter :: MinVxy = 1D-40
  double precision, parameter :: MinLenFrac = 1D-6
  !
  ! For intesection with top and bottom surfaces
  if (abs(ray%vz) .ge. MinVz) then
    L(1) = (c%ymax - ray%z) / ray%vz
    L(2) = (c%ymin - ray%z) / ray%vz
  else
    L(1) = FL
    L(2) = FL
  end if
  if (L(1) .ge. 0D0) then
    t1 = ray%x + L(1)*ray%vx
    t2 = ray%y + L(1)*ray%vy
    rr = t1 * t1 + t2 * t2
    if ((rr .lt. c%xmin*c%xmin) .or. (rr .gt. c%xmax*c%xmax)) then
      L(1) = FL
    end if
  end if
  if (L(2) .ge. 0D0) then
    t1 = ray%x + L(2)*ray%vx
    t2 = ray%y + L(2)*ray%vy
    rr = t1 * t1 + t2 * t2
    if ((rr .lt. c%xmin*c%xmin) .or. (rr .gt. c%xmax*c%xmax)) then
      L(2) = FL
    end if
  end if
  !
  ! For intersection with the inner and outer cylindrical boundaries
  A = ray%vx * ray%vx + ray%vy * ray%vy
  B = 2D0 * (ray%x*ray%vx + ray%y*ray%vy)
  t1 = ray%x * ray%x + ray%y * ray%y
  C1 = t1 - c%xmin * c%xmin
  C2 = t1 - c%xmax * c%xmax
  t1 = B * B
  t2 = 4D0 * A
  D1 = t1 - t2 * C1
  D2 = t1 - t2 * C2
  ! Inner cylinder
  if (D1 .gt. 0D0) then
    if (abs(A) .gt. MinVxy) then
      t1 = sqrt(D1)
      L(3) = (-B + t1) / (2D0*A)
      zz = ray%z + ray%vz * L(3)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(3) = FL
      end if
      L(4) = (-B - t1) / (2D0*A)
      zz = ray%z + ray%vz * L(4)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(4) = FL
      end if
    else
      L(3) = FL
      L(4) = FL
    end if
  else
    L(3) = FL
    L(4) = FL
  end if
  ! Outer cylinder
  if (D2 .gt. 0D0) then
    if (abs(A) .gt. MinVxy) then
      t1 = sqrt(D2)
      L(5) = (-B + t1) / (2D0*A)
      zz = ray%z + ray%vz * L(5)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(5) = FL
      end if
      L(6) = (-B - t1) / (2D0*A)
      zz = ray%z + ray%vz * L(6)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(6) = FL
      end if
    else
      L(5) = FL
      L(6) = FL
    end if
  else
    L(5) = FL
    L(6) = FL
  end if
  ! The closest one is what we want.
  rr   = 1D100
  idx = 0
  do i=1, 6
    if ((L(i) .gt. MinLen) .and. (L(i) .lt. rr)) then
      rr = L(i)
      idx = i
    end if
  end do
  if (idx .eq. 0) then
    found = .false.
    dirtype = Undef_dirtype
  else
    found = .true.
    length = L(idx)
    !eps = eps_ratio * max(min(c%xmax-c%xmin, c%ymax-c%ymin), L(idx))
    !eps = eps_ratio * max(min(c%xmax-c%xmin, c%ymax-c%ymin), 1D-2*L(idx))
    !eps = eps_ratio * min(c%xmax-c%xmin, c%ymax-c%ymin)
    eps = min(c%xmax-c%xmin, c%ymax-c%ymin) * MinLenFrac
    L(idx) = L(idx) + eps
    t1 = ray%x + ray%vx * L(idx)
    t2 = ray%y + ray%vy * L(idx)
    rsq = t1 * t1 + t2 * t2  ! r squared
    z = ray%z + ray%vz * L(idx)
    flag_inside_cell = &
      (c%xmin**2 .le. rsq) .and. &
      (c%xmax**2 .ge. rsq) .and. &
      (c%ymin    .le. z) .and. &
      (c%ymax    .ge. z)
    if (.not. flag_inside_cell) then
      ! Exit the cell c instead of entering c
      ! 1: Exit through the top
      ! 2: Exit through the bottom
      ! 3: Exit through the inner edge
      ! 4: Exit through the inner edge also
      ! 5: Exit through the outer edge
      ! 6: Exit through the outer edge also
      dirtype = idx
    else
      !!! found = .false.
      dirtype = -idx
    end if
  end if
end subroutine calc_intersection_ray_cell



pure function is_inside_cell_sq(rsq, z, c)
  logical is_inside_cell_sq
  double precision, intent(in) :: rsq, z
  type(type_cell), pointer, intent(in) :: c
  is_inside_cell_sq = &
    (c%xmin**2 .le. rsq) .and. (c%xmax**2 .ge. rsq) .and. &
    (c%ymin .le. z) .and. (c%ymax .ge. z)
end function is_inside_cell_sq



pure function is_inside_cell_mirror_sq(rsq, z, c)
  logical is_inside_cell_mirror_sq
  double precision, intent(in) :: rsq, z
  type(type_cell), pointer, intent(in) :: c
  is_inside_cell_mirror_sq = &
    (c%xmin**2 .le. rsq) .and. (c%xmax**2 .ge. rsq) .and. &
    (c%ymin .le. abs(z)) .and. (c%ymax .ge. abs(z))
end function is_inside_cell_mirror_sq




end module ray_propagating
