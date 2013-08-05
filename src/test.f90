program test
use quick_sort
use grid
implicit none

!integer, parameter :: n=3, m=10, ncmp=3
!double precision, dimension(n, m) :: a = &
!  reshape((/1D0, 2D0, -5D0, &
!            0D0, 2D0, -4D0, &
!            4D0, 2D0, -2D0, &
!            9D0, 5D0, -3D0, &
!            3D0, 4D0, -1D0, &
!            0D0, 2D0, -2D0, &
!            3D0, 4D0, -6D0, &
!            1D0, 5D0, -2D0, &
!            0D0, 5D0, -3D0, &
!            2D0, 0D0, -1D0  /), (/n, m/))
!
!integer, dimension(ncmp) :: icmp = (/2, 1, 2/)

integer, dimension(4) :: xx

type :: test_memory_
 double precision, dimension(10) :: z = 1D0
end type test_memory_

type :: test_memory
  integer m
  !type(test_memory_), dimension(1024, 1024, 1024) :: x
  type(test_memory_), dimension(100, 100) :: x
end type test_memory

type :: test_memory_x
  integer n
  type(test_memory), dimension(1) :: x
end type test_memory_x

integer i, j

double precision, dimension(100,1000,1000) :: zz

type(test_memory_x), pointer :: yy

allocate(yy)
do i=1,100
do j=1,100
yy%x(1)%x(i,j)%z = dble(i*j)
end do
end do

zz = 0D0

!yy%n = 12
!yy%x%m = 23
!yy%x%x(1024,1024,1)%z = 2D2

i = 1
associate(qq => xx(i), ww => i)
  ww = 2
  write(*,*) ww
  i = 4
  write(*,*) ww
  qq = 23
end associate
write(*,*) xx
write(*,*) i

read(*,*)

!write(*, *) 'Before'
!do i=1, m
!  write(*, '(3F7.2)') a(:, i)
!end do
!call quick_sort_array(a, n, m, ncmp, icmp)
!write(*, *) 'After'
!do i=1, m
!  write(*, '(3F7.2)') a(:, i)
!end do

grid_config%rmin = 0.5D0
grid_config%rmax = 50D0
grid_config%zmin = 0.0D0
grid_config%zmax = 50D0


call make_grid
!do i=1, refinement_data%nlen
!  write(*, '(3ES12.4)') refinement_data%xyv(:, i)
!end do
write(*, '(A, ES14.4)') 'Max density: ', get_density_analytic(root%xmin, root%ymin)
write(*,*) root%nOffspring
write(*,*) cell_leaves%nlen

!do i=1, cell_leaves%nlen
!  if ((cell_leaves%list(i)%p%above%n .ge. 1) .and. &
!      (cell_leaves%list(i)%p%below%n .ge. 1) .and. &
!      (cell_leaves%list(i)%p%inner%n .ge. 1) .and. &
!      (cell_leaves%list(i)%p%outer%n .ge. 1)) cycle
!  write(*, '(I6, 5ES14.5, 4I4)') &
!       cell_leaves%list(i)%p%order, &
!       cell_leaves%list(i)%p%par%n_gas, &
!       cell_leaves%list(i)%p%xmin, &
!       cell_leaves%list(i)%p%ymin, &
!       cell_leaves%list(i)%p%xmax, &
!       cell_leaves%list(i)%p%ymax, &
!       cell_leaves%list(i)%p%inner%n, &
!       cell_leaves%list(i)%p%outer%n, &
!       cell_leaves%list(i)%p%above%n, &
!       cell_leaves%list(i)%p%below%n
!end do

do i=1, surf_cells%nlen
  j = surf_cells%idx(i)
  write(*, '(I6, 5ES14.5, 4I4)') &
       cell_leaves%list(j)%p%order, &
       cell_leaves%list(j)%p%val, &
       cell_leaves%list(j)%p%xmin, &
       cell_leaves%list(j)%p%ymin, &
       cell_leaves%list(j)%p%xmax, &
       cell_leaves%list(j)%p%ymax, &
       cell_leaves%list(j)%p%inner%n, &
       cell_leaves%list(j)%p%outer%n, &
       cell_leaves%list(j)%p%above%n, &
       cell_leaves%list(j)%p%below%n
end do
do i=1, bott_cells%nlen
  j = bott_cells%idx(i)
  write(*, '(I6, 5ES14.5, 4I4)') &
       cell_leaves%list(j)%p%order, &
       cell_leaves%list(j)%p%val, &
       cell_leaves%list(j)%p%xmin, &
       cell_leaves%list(j)%p%ymin, &
       cell_leaves%list(j)%p%xmax, &
       cell_leaves%list(j)%p%ymax, &
       cell_leaves%list(j)%p%inner%n, &
       cell_leaves%list(j)%p%outer%n, &
       cell_leaves%list(j)%p%above%n, &
       cell_leaves%list(j)%p%below%n
end do


end program test

