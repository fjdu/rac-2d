module binary_array_io

use trivials
use spline_1d_2d

implicit none


type :: type_table_2d
  integer :: nx, ny
  double precision, dimension(:), allocatable :: x, y
  double precision, dimension(:,:), allocatable :: val
end type type_table_2d


contains


subroutine read_binary_array(filename, res)
  character(len=*), intent(in) :: filename
  type(type_table_2d), intent(out) :: res
  integer i, fU, ios
  double precision ndim
  double precision, dimension(:), allocatable :: dims
  !
  if (.not. getFileUnit(fU)) then
     write(*, *) 'In read_binary_array:'
     write(*, *) 'Cannot get a file unit.'
     return
  end if
  open(UNIT=fU, FILE=filename, &
     IOSTAT=ios, &
     ACCESS='STREAM', &
     FORM='UNFORMATTED', ACTION='READ')
  read(fU) ndim
  !write(*,*) 'ndim=', ndim
  allocate(dims(int(ndim)))
  do i=1, int(ndim)
    read(fU) dims(i)
  end do
  res%nx = int(dims(1))
  res%ny = int(dims(2))
  !write(*,*) 'nx,ny=', res%nx, res%ny
  if (allocated(res%x)) then
    deallocate(res%x)
  end if
  if (allocated(res%y)) then
    deallocate(res%y)
  end if
  if (allocated(res%val)) then
    deallocate(res%val)
  end if
  allocate(res%x(res%nx), res%y(res%ny), &
           res%val(res%nx, res%ny))
  read(fU) res%x
  read(fU) res%y
  read(fU) res%val
  close(fU)
end subroutine read_binary_array



subroutine create_spline2d_from_table(tab, spl)
  type(type_table_2d), intent(in) :: tab
  type(type_spline_2D), intent(out) :: spl
  spl%nx = tab%nx
  spl%ny = tab%ny
  spl%itype = 0
  !
  !if (allocated(spl%xi)) then
  !  deallocate(spl%xi, spl%yi, spl%vi)
  !end if
  allocate(spl%xi(spl%nx), &
           spl%yi(spl%ny), &
           spl%vi(spl%nx, spl%ny))
  spl%xi = tab%x
  spl%yi = tab%y
  spl%vi = tab%val
  call spline2d_prepare(spl)
end subroutine create_spline2d_from_table


subroutine create_spline2d_from_file(filename, spl)
  character(len=*), intent(in) :: filename
  type(type_spline_2D), intent(out) :: spl
  type(type_table_2d) :: tab
  call read_binary_array(filename, tab)
  call create_spline2d_from_table(tab, spl)
  deallocate(tab%x, tab%y, tab%val)
end subroutine create_spline2d_from_file


subroutine get_idx_in_table(ix, iy, x0, y0, tab, is_regular)
  integer, intent(out) :: ix, iy
  type(type_table_2d), intent(in) :: tab
  double precision, intent(in) :: x0, y0
  logical, intent(in) :: is_regular
  double precision del
  if (is_regular) then
    del = tab%x(2) - tab%x(1)
    ix = floor((x0 - tab%x(1)) / del) + 1
    del = tab%y(2) - tab%y(1)
    iy = floor((y0 - tab%y(1)) / del) + 1
  else
    if (tab%x(1) .gt. x0) then
      ix = 0
    else if (tab%x(tab%nx) .lt. x0) then
      ix = tab%nx + 1
    else
      do ix=1, tab%nx
        if (tab%x(ix) .le. x0) then
          exit
        end if
      end do
    end if
    if (tab%y(1) .gt. y0) then
      iy = 0
    else if (tab%y(tab%ny) .lt. y0) then
      iy = tab%ny + 1
    else
      do iy=1, tab%ny
        if (tab%y(iy) .le. y0) then
          exit
        end if
      end do
    end if
  end if
end subroutine get_idx_in_table



end module binary_array_io
