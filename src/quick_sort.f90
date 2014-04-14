module quick_sort
implicit none

public :: quick_sort_array
private :: partition, LE_vec, GE_vec

contains


subroutine unique_vector_idx(a, n, idx_unique, n_unique, rtol, atol, idx_reverse)
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: a
  integer, intent(out) :: n_unique
  integer, dimension(n), intent(out) :: idx_unique 
  double precision, intent(in), optional :: rtol, atol
  integer, dimension(n), intent(out), optional :: idx_reverse 
  integer, dimension(n) :: idx_sorted 
  integer i, i0, i1
  double precision rt, at
  !
  if (present(rtol)) then
    rt = rtol
  else
    rt = 0D0
  end if
  !
  if (present(atol)) then
    at = atol
  else
    at = 0D0
  end if
  !
  call  quick_sort_vector_idx(a, n, idx_sorted)
  idx_unique(1) = idx_sorted(1)
  n_unique = 1
  if (present(idx_reverse)) then
    idx_reverse(idx_sorted(1)) = 1
  end if
  do i=2, n
    i0 = idx_unique(n_unique)
    i1 = idx_sorted(i)
    if (abs(a(i0) - a(i1)) .gt. (0.5D0*rt*(a(i0) + a(i1)) + at)) then
      n_unique = n_unique + 1
      idx_unique(n_unique) = i1
    end if
    !
    if (present(idx_reverse)) then
      idx_reverse(i1) = n_unique
    end if
  end do
end subroutine unique_vector_idx


subroutine unique_vector(a, n, n_unique, rtol, atol)
  integer, intent(in) :: n
  double precision, dimension(n), intent(inout) :: a
  integer, intent(out) :: n_unique
  double precision, intent(in), optional :: rtol, atol
  double precision, dimension(1, n) :: atmp
  integer i
  double precision rt, at
  !
  if (present(rtol)) then
    rt = rtol
  else
    rt = 0D0
  end if
  !
  if (present(atol)) then
    at = atol
  else
    at = 0D0
  end if
  !
  atmp(1, :) = a
  call quick_sort_array(atmp, 1, n, 1, (/1/))
  a(1) = atmp(1, i)
  n_unique = 1
  do i=2,n
    if (abs(a(i-1) - atmp(1, i)) .gt. &
        (0.5D0*rt*(a(i-1) + atmp(1, i)) + at)) then
      n_unique = n_unique + 1
      a(n_unique) = a(i)
    end if
  end do
end subroutine unique_vector


pure subroutine quick_sort_vector_idx(a, n, idx_sorted)
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: a
  integer, dimension(n), intent(out) :: idx_sorted
  double precision, dimension(2, n) :: atmp
  integer i
  do i=1, n
    atmp(1, i) = a(i)
    atmp(2, i) = dble(i)
  end do
  call quick_sort_array(atmp, 2, n, 1, (/1/))
  do i=1, n
    idx_sorted(i) = int(atmp(2, i))
  end do
end subroutine quick_sort_vector_idx


pure recursive subroutine quick_sort_array(a, n, m, ncmp, icmp)
  ! Array to be sorted: a
  ! Array a has n columns and m rows.
  ! The sorting is performed along the column direction, namely,
  ! each row is treated as a whole.
  ! icmp contains the indices of the columns to be compared.
  ! ncmp = len(icmp)
  integer, intent(in) :: n, m, ncmp
  double precision, dimension(n, m), intent(inout) :: a
  integer, dimension(ncmp), intent(in) :: icmp
  integer ip
  if (m > 1) then
    call partition(a, n, m, ncmp, icmp, ip)
    if (ip > 2) then
      call quick_sort_array(A(:, 1:ip-1), n, ip-1, ncmp, icmp)
    end if
    if (m-ip > 0) then
      call quick_sort_array(A(:, ip:m), n, m-ip+1, ncmp, icmp)
    end if
  end if
end subroutine quick_sort_array


pure subroutine partition(a, n, m, ncmp, icmp, marker)
  integer, intent(in) :: n, m, ncmp
  double precision, dimension(n, m), intent(inout) :: a
  integer, dimension(ncmp), intent(in) :: icmp
  integer, intent(out) :: marker
  integer i, j
  double precision, dimension(n) :: x, tmp
  x = a(:, 1)
  i = 0
  j = m + 1
  do
    i = i + 1
    do
      if (GE_vec(a(:, i), x, ncmp, icmp)) then
        exit
      else
        i = i + 1
      end if
    end do
    j = j - 1
    do
      if (LE_vec(a(:, j), x, ncmp, icmp)) then
        exit
      else
        j = j - 1
      end if
    end do
    if (i < j) then
      tmp = a(:, i)
      a(:, i) = a(:, j)
      a(:, j) = tmp
    else if (i == j) then
      marker = i + 1
      return
    else
      marker = i
      return
    end if
  end do
end subroutine partition


pure function LE_vec(x, y, ncmp, icmp)
  ! Return true if x <= y in the lexical sense
  logical LE_vec
  double precision, dimension(:), intent(in) :: x, y
  integer, intent(in) :: ncmp
  integer, dimension(:), intent(in) :: icmp
  integer i, j
  do i=1, ncmp
    j = icmp(i)
    if (x(j) .lt. y(j)) then
      LE_vec = .true.
      return
    else if (x(j) .gt. y(j)) then
      LE_vec = .false.
      return
    else
      LE_vec = .true.
    end if
  end do
end function LE_vec


pure function GE_vec(x, y, ncmp, icmp)
  logical GE_vec
  double precision, dimension(:), intent(in) :: x, y
  integer, intent(in) :: ncmp
  integer, dimension(:), intent(in) :: icmp
  integer i, j
  do i=1, ncmp
    j = icmp(i)
    if (x(j) .gt. y(j)) then
      GE_vec = .true.
      return
    else if (x(j) .lt. y(j)) then
      GE_vec = .false.
      return
    else
      GE_vec = .true.
    end if
  end do
end function GE_vec

end module quick_sort
