module cell_class

use trivials
use data_struct
use phy_const

implicit none

integer, parameter :: n_cell_split = 2
integer, parameter :: nMaxOrder = 64
integer, parameter :: nMaxNumOfNeighbor = 32

double precision, parameter :: &
  vel_sample_range_in_doppler = 4D0, &
  vel_sample_intvl_in_doppler = 0.2D0

! To be set by the user
integer ::            n_cell_par = 11
double precision ::   smallest_size = 1D-4
double precision ::   max_ratio_within_cell = 1.5D0
double precision ::   min_den_considered = 1D0
double precision ::   aspect_ratio_max = 1D0
double precision ::   max_cell_size_frac = 0.5D0
double precision ::   max_cell_size_0 =    0.5D0
integer ::            n_closest_points_use = 4
double precision ::   factor_for_local_scale = 2D0
namelist /cell_config/ &
  n_cell_par, &
  smallest_size, &
  max_ratio_within_cell, &
  min_den_considered, &
  aspect_ratio_max, &
  max_cell_size_frac, &
  max_cell_size_0, &
  n_closest_points_use, &
  factor_for_local_scale

! To be set by the code
double precision very_small_len
double precision very_small_den

integer :: n_cell_leaves = 0

type :: cell_par
  double precision :: r=0D0, z=0D0, n=0D0, eps=0D0, sc=0D0, ab=0D0, &
    tot=0D0, alb=0D0, g=0D0, Tgas=0D0, Tdust=0D0, &
    coeff_len2tottau=0D0, coeff_len2abstau = 0D0, coeff_len2scatau=0D0, &
    surf_a=0D0, vol=0D0, dv_doppler=0D0, df_doppler=0D0, &
    fH = 0D0, velocity=0D0
  integer n_vel_intvl
  double precision, dimension(:), allocatable :: v_thermal, gaussVal, Voigt
  double precision sca_Lyman_alpha
end type cell_par

type :: cell_rad
  integer :: counts = 0
  double precision :: J=0D0, photon_den=0D0
  double precision, dimension(:), allocatable :: Jnu
  double precision, dimension(:), allocatable :: Fnu
end type cell_rad

type :: cell_neighbors
  integer :: nUp=0, nDown=0, nIn=0, nOut=0
  integer, dimension(nMaxNumOfNeighbor) :: iUp, iDown, iIn, iOut
end type cell_neighbors

type :: cell
  double precision x1, x2, y1, y2
  integer order
  integer nOffspr
  integer :: nChildren
  type(cell), pointer :: parent
  type(cell), pointer, dimension(:) :: children
  double precision, dimension(:), allocatable :: d
  ! d(1): r; d(2): z; d(3): nH; d(4): eps;
  ! d(5): sca; d(6): abs; d(7): g; d(8): fH
  ! d(9): Tgas; d(10): fH2O; d(11): fOH
  type(cell_par), pointer :: par
  type(cell_rad), pointer :: rad
  type(cell_neighbors), pointer :: nei
  !type(chem_single_point), pointer :: chem
end type cell

type :: cell_ptr
  type(cell), pointer :: p
end type cell_ptr

type :: polygon_generic
  integer n_vert
  double precision, dimension(:,:), allocatable :: xy
end type polygon_generic

type(cell_ptr), dimension(:), allocatable :: cell_leaves
type(polygon_generic) data_inputs_bd

contains

  subroutine initialize_cell_class
    implicit none
    double precision, parameter :: very_small_factor = 1D-3
    very_small_len = very_small_factor * smallest_size
    very_small_den = very_small_factor * min_den_considered
  end subroutine initialize_cell_class


  subroutine set_data_inputs_bd()
    implicit none
    double precision minr, minz, maxr, maxz
    minr = different_r(1)
    maxr = different_r(n_different_r)
    minz = minval(data_inputs(3, :))
    maxz = maxval(data_inputs(3, :))
    data_inputs_bd%n_vert = 5
    allocate(data_inputs_bd%xy(2,data_inputs_bd%n_vert))
    data_inputs_bd%xy(:,1) = (/minr, minz/)
    data_inputs_bd%xy(:,2) = (/maxr, minz/)
    data_inputs_bd%xy(:,3) = (/maxr, maxz/)
    data_inputs_bd%xy(:,4) = (/minr, 0.15D0/)
    data_inputs_bd%xy(:,5) = (/minr, minz/)
  end subroutine set_data_inputs_bd


  logical function is_inside_polygon(p, x, y)
    implicit none
    double precision x, y
    type(polygon_generic) p
    integer i, j
    is_inside_polygon = .FALSE.
    i = 1
    j = p%n_vert
    do
      if (i .GT. p%n_vert) exit
      if (&
        (((p%xy(2,i) .LT. y) .AND. (p%xy(2,j) .GE. y)) .OR. &
         ((p%xy(2,j) .LT. y) .AND. (p%xy(2,i) .GE. y))) .AND. &
        ((p%xy(1,i) .LE. x) .OR. (p%xy(1,j) .LE. x))) then
        is_inside_polygon = xor(is_inside_polygon, &
          ((p%xy(1,i) + (y-p%xy(2,i)) * (p%xy(1,j)-p%xy(1,i)) / &
                   (p%xy(2,j)-p%xy(2,i))) .LT. x))
      end if
      j = i
      i = i + 1
    end do
  end function is_inside_polygon


  !!logical function is_inside_polygon(p, x, y)
  !!  implicit none
  !!  double precision x, y
  !!  type(polygon_generic) p
  !!  integer i
  !!  complex(kind=8) z1, z2, dz
  !!  double precision angle_sum
  !!  double precision, parameter :: angle_tol = 1D0
  !!  angle_sum = 0D0
  !!  z1 = cmplx(p%xy(1,1)-x, p%xy(2,1)-y)
  !!  do i=1, p%n_vert-1
  !!    z2 = cmplx(p%xy(1,i+1)-x, p%xy(2,i+1)-y)
  !!    dz = z2/z1
  !!    angle_sum = angle_sum + atan2(aimag(dz), real(dz))
  !!    z1 = z2
  !!  end do
  !!  if (abs(angle_sum) .LT. angle_tol) then
  !!    is_inside_polygon = .FALSE.
  !!  else
  !!    is_inside_polygon = .TRUE.
  !!  end if
  !!end function is_inside_polygon

  subroutine set_bounding_box(this)
    implicit none
    type(cell), intent(inout), target :: this
    this%x1  = 0D0
    this%x2  = 5D2!different_r(n_different_r)
    this%y1  = 0D0!minval(data_inputs(3, :))
    this%y2  = 6.5D2!maxval(data_inputs(3, :))
  end subroutine set_bounding_box


  subroutine cell_initialize(this)
    implicit none
    type(cell), intent(inout), target :: this
    allocate(this%d(n_cell_par))
    nullify(this%parent)
    nullify(this%children)
    !
    this%d = dbl_NaN
    this%x1  = dbl_NaN
    this%x2  = dbl_NaN
    this%y1  = dbl_NaN
    this%y2  = dbl_NaN
    this%order = 1
    this%nOffspr = 0
    this%nChildren = 0
  end subroutine cell_initialize


  recursive subroutine make_children(this)
    implicit none
    type(cell), intent(inout), target :: this
    integer i
    integer idxr1, idxr2, i_r, ir, idxz1, idxz2, iz, idx
    double precision min_den, max_den
    double precision r_sep_local, z_sep_local
    integer n_points_in_cell
    logical uniform_in_cell, return_now
    double precision xmid, ymid
    double precision del_x_1, del_x_2, del_y_1, del_y_2, aspect_1, aspect_2
    ! Firt find out the data points that fall into this cell.
    ! The data of this cell is taken to be the average of these points.
    min_den = huge(0D0)
    max_den = 0D0
    this%d = 0D0
    n_points_in_cell = 0
    return_now = .FALSE.
    uniform_in_cell = .TRUE.
    idxr1 = binary_search(different_r, n_different_r, this%x1, 1)
    idxr2 = binary_search(different_r, n_different_r, this%x2, 1)
    do i_r = idxr1, idxr2
      ir = (i_r - 1) * n_z_step + 1
      idxz1 = binary_search(&
        data_inputs(3, ir:ir+n_z_step-1), n_z_step, this%y1, 1)
      idxz2 = binary_search(&
        data_inputs(3, ir:ir+n_z_step-1), n_z_step, this%y2, 1)
      do iz = idxz2, idxz1 ! idxz1 .GT. idxz2
        idx = ir + iz - 1
        if ((data_inputs(2, idx) .LT. this%x1-very_small_len) .OR. &
            (data_inputs(2, idx) .GT. this%x2+very_small_len) .OR. &
            (data_inputs(3, idx) .LT. this%y1-very_small_len) .OR. &
            (data_inputs(3, idx) .GT. this%y2+very_small_len)) then
          cycle
        end if
        n_points_in_cell = n_points_in_cell + 1
        min_den = min(min_den, data_inputs(4, idx))
        max_den = max(max_den, data_inputs(4, idx))
        this%d = this%d + data_inputs(2:n_cell_par+1, idx)
      end do
      if ((max_den .GT. min_den_considered) .AND. &
          ((max_den / (min_den+very_small_den)) .GT. max_ratio_within_cell)) then
        uniform_in_cell = .FALSE.
      end if
    end do
    this%d = this%d / dble(n_points_in_cell)
    ! If the cell is uniform enough, this cell will not be further splitted
    ! unless its aspect ratio is too large.
    if (n_points_in_cell .EQ. 0) then
      write(*,*) 'This cell does not contain any input points!'
      return_now = .TRUE.
    end if
    if (uniform_in_cell .OR. (n_points_in_cell .EQ. 1)) then
      !if (.NOT. is_inside_polygon(data_inputs_bd, &
      !     (this%x1+this%x2)*0.5D0, (this%y1+this%y2)*0.5D0)) then
      !  this%d(3:n_cell_par) = 0D0
      !end if
      return_now = .TRUE.
    end if
    if (this%order .GE. nMaxOrder) then
      write(*,*) 'Max tree depth reached!'
    end if
    if (return_now) then
      call regularize_cell(this)
      if (this%nChildren .GT. 0) then
        this%nOffspr = this%nChildren
        n_cell_leaves = n_cell_leaves + this%nChildren
      else
        nullify(this%children)
        this%nChildren = 0
        this%nOffspr = 0
        n_cell_leaves = n_cell_leaves + 1
      end if
      return
    end if
    ! This cell will be splitted.
    this%nChildren = n_cell_split
    allocate(this%children(this%nChildren))
    do i=1, this%nChildren
      call cell_initialize(this%children(i))
    end do
    ! The split point is taken to be the geometric center of the points
    ! within this cell.  It is a binary split.  The direction of split is
    ! chosen such that the resulting aspect ratio is as small as possible.
    xmid = this%d(1)
    ymid = this%d(2)
    del_x_1 = xmid-this%x1
    del_x_2 = this%x2-xmid
    del_y_1 = ymid-this%y1
    del_y_2 = this%y2-ymid
    if ((min(del_x_1, del_x_2) .LE. very_small_len) .AND. &
        (min(del_y_1, del_y_2) .LE. very_small_len)) &
      write(*,'(A, 10ES10.3, I5)') 'Something is wrong! (1)', &
        del_x_1, del_x_2, del_y_1, del_y_2, &
        this%x1, this%x2, xmid, this%y1, this%y2, ymid, n_points_in_cell
    aspect_1 = max(max((this%y2-this%y1), del_x_1) / min((this%y2-this%y1), del_x_1), &
                   max((this%y2-this%y1), del_x_2) / min((this%y2-this%y1), del_x_2))
    aspect_2 = max(max((this%x2-this%x1), del_y_1) / min((this%x2-this%x1), del_y_1), &
                   max((this%x2-this%x1), del_y_2) / min((this%x2-this%x1), del_y_2))
    if ((aspect_1 .LE. aspect_2) .AND. (min(del_x_1, del_x_2) .GT. very_small_len)) then
      this%children(1)%x1 = this%x1
      this%children(1)%x2 = xmid
      this%children(1)%y1 = this%y1
      this%children(1)%y2 = this%y2
      this%children(2)%x1 = xmid
      this%children(2)%x2 = this%x2
      this%children(2)%y1 = this%y1
      this%children(2)%y2 = this%y2
    else
      this%children(1)%x1 = this%x1
      this%children(1)%x2 = this%x2
      this%children(1)%y1 = this%y1
      this%children(1)%y2 = ymid
      this%children(2)%x1 = this%x1
      this%children(2)%x2 = this%x2
      this%children(2)%y1 = ymid
      this%children(2)%y2 = this%y2
    end if
    do i=1, this%nChildren
      this%children(i)%parent => this
      this%children(i)%order = this%order + 1
      call make_children(this%children(i))
      this%nOffspr = this%nOffspr + this%children(i)%nOffspr + 1
    end do
  end subroutine make_children


  subroutine regularize_cell(this)
    ! To further split cells that are too thin or too fat into smaller ones.
    ! To be called when a new terminate cell is constructed.
    implicit none
    type(cell), intent(inout), target :: this
    double precision aspect_ratio, delta
    integer i, j, iSub, nsplit
    double precision del_max_x, del_max_y
    double precision del_x, del_y, del1_x, del1_y
    integer ndiv_x, ndiv_y
    ! Max size allowed at this position
    del_max_x = max_cell_size_frac * sqrt(this%x1*this%x1 + this%y1*this%y1) + max_cell_size_0
    del_max_y = del_max_x
    del_x = this%x2 - this%x1
    del_y = this%y2 - this%y1
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
    if (nsplit .EQ. 1) return
    del_x = del_x / dble(ndiv_x)
    del_y = del_y / dble(ndiv_y)
    this%nChildren = nsplit
    allocate(this%children(nsplit))
    iSub = 1
    do i=0, ndiv_x-1
      do j=0, ndiv_y-1
        call cell_initialize(this%children(iSub))
        this%children(iSub)%parent => this
        this%children(iSub)%order = this%order + 1
        this%children(iSub)%x1 = this%x1 + del_x * dble(i)
        this%children(iSub)%x2 = this%children(iSub)%x1 + del_x
        this%children(iSub)%y1 = this%y1 + del_y * dble(j)
        this%children(iSub)%y2 = this%children(iSub)%y1 + del_y
        ! Get rid of the rounding errors.
        if (i .EQ. 0)          this%children(iSub)%x1 = this%x1
        if (i .EQ. (ndiv_x-1)) this%children(iSub)%x2 = this%x2
        if (j .EQ. 0)          this%children(iSub)%y1 = this%y1
        if (j .EQ. (ndiv_y-1)) this%children(iSub)%y2 = this%y2
        iSub = iSub + 1
      end do
    end do
    do i=1, nsplit
      call set_cell_data(this%children(i)) ! This is slow!
    end do
  end subroutine regularize_cell


  subroutine set_cell_data(this)
    type(cell), intent(inout), target :: this
    double precision x, y, dtmp, local_scale
    integer i, j
    integer, dimension(n_closest_points_use) :: idx
    double precision, dimension(n_closest_points_use) :: dist
    x =  0.5D0*(this%x1 + this%x2)
    y =  0.5D0*(this%y1 + this%y2)
    ! First find out the nearest neighbors.
    dist = huge(0D0)
    idx = 0
    do i=1, len_inputs
      dtmp = (x-data_inputs(2,i))*(x-data_inputs(2,i)) + (y-data_inputs(3,i))*(y-data_inputs(3,i))
      do j=1, n_closest_points_use
        if (dtmp .LT. dist(j)) then
          dist(j) = dtmp
          idx(j) = i
          exit
        end if
      end do
    end do
    ! Now calculate the local length scale, which is defined to be the
    ! longest among the distances between points that are close to the given
    ! point.
    local_scale = 0D0
    do i=1, n_closest_points_use
      do j=i+1, n_closest_points_use
        if ((idx(i) .LE. 0) .OR. (idx(j) .LE. 0)) cycle
        local_scale = max(local_scale, &
               (data_inputs(2, idx(i)) - data_inputs(2, idx(j))) * &
               (data_inputs(2, idx(i)) - data_inputs(2, idx(j))) + &
               (data_inputs(3, idx(i)) - data_inputs(3, idx(j))) * &
               (data_inputs(3, idx(i)) - data_inputs(3, idx(j))))
      end do
    end do
    ! Now set the data values.  Note that the dist and local_scale are not square-rooted.
    if (dist(1) .LE. local_scale*(factor_for_local_scale*factor_for_local_scale)) then
      this%d(3:n_cell_par) = data_inputs(4:n_cell_par+1, idx(1))
    else
      this%d = 0D0
    end if
  end subroutine set_cell_data


  recursive subroutine make_cell_leaves(this, idx)
    implicit none
    type(cell), intent(inout), target :: this
    integer i, idx
    double precision dv
    if (.NOT. allocated(cell_leaves)) then
      idx = 1
      allocate(cell_leaves(n_cell_leaves))
    end if
    if (this%nOffspr .EQ. 0) then
      cell_leaves(idx)%p => this
      idx = idx + 1
      allocate(this%par, this%rad)
      this%par%r = this%d(1)
      this%par%z = this%d(2)
      this%par%surf_a = 2D0*phy_Pi * (this%x2*this%x2 - this%x1*this%x1) + &
        2D0*phy_Pi * (this%x2 - this%x1) * (this%y2 - this%y1)
      this%par%vol = (this%x2 - this%x1) * (this%x2 + this%x1) &
          * (this%y2 - this%y1) * phy_Pi * 2D0 ! The factor 2D0 is due to the reflection symmetry.
      !!! Physical parameters; they can be manipulated for testing purpose.
      this%par%n      = this%d(3)
      this%par%eps    = this%d(4)
      this%par%sc     = this%d(5)
      this%par%ab     = this%d(6)
      this%par%tot    = this%par%sc + this%par%ab
      if (this%par%tot .LE. tiny(0D0)) then
        this%par%alb  = 1D0 ! Albedo = 1 for vacuum.
      else
        this%par%alb  = this%par%sc / this%par%tot
      end if
      this%par%g      = this%d(7)
      this%par%Tgas   = this%d(9)
      this%par%Tdust  = this%d(9)
      this%par%fH     = this%d(10)
      this%par%dv_doppler = sqrt(phy_kBoltzmann*this%par%Tgas/phy_mProton) !
      this%par%df_doppler = this%par%dv_doppler / phy_SpeedOfLight &
          * Lyman_alpha_freq
      this%par%coeff_len2tottau = phy_AU2cm * this%par%n * this%par%eps * this%par%tot
      this%par%coeff_len2abstau = phy_AU2cm * this%par%n * this%par%eps * this%par%ab
      this%par%coeff_len2scatau = phy_AU2cm * this%par%n * this%par%eps * this%par%sc
      this%par%velocity = sqrt(phy_GravitationConst * (phy_Msun*star_mass_in_Msun) / &
          (this%par%r*phy_AU2m))
      !
      this%par%n_vel_intvl = &
          1 + 2*int(vel_sample_range_in_doppler / vel_sample_intvl_in_doppler)
      allocate(this%par%v_thermal(this%par%n_vel_intvl), &
               this%par%gaussVal(this%par%n_vel_intvl), &
               this%par%Voigt(this%par%n_vel_intvl))
      dv = this%par%dv_doppler * vel_sample_intvl_in_doppler
      this%par%v_thermal(1) = -vel_sample_range_in_doppler * this%par%dv_doppler
      do i=2, this%par%n_vel_intvl
        this%par%v_thermal(i) = this%par%v_thermal(i-1) + dv
      end do
      this%par%gaussVal(1) = 0D0
      do i=1, this%par%n_vel_intvl
        this%par%gaussVal(i) = 1D0 / (phy_sqrt2Pi * this%par%dv_doppler) &
            * exp(-this%par%v_thermal(i)*this%par%v_thermal(i) &
                  /(2D0*this%par%dv_doppler*this%par%dv_doppler))
      end do
      ! Vacuum cells won't scatter or absorb photons.
      if (this%par%n .LE. min_den_considered) then
        this%par%eps    = 0D0
        this%par%sc     = 0D0
        this%par%ab     = 0D0
        this%par%tot    = 0D0
        this%par%alb    = 1D0
        this%par%g      = 0D0
        this%par%Tgas   = 0D0
        this%par%Tdust  = 0D0
        this%par%fH     = 0D0
        this%par%dv_doppler = 0D0
        this%par%df_doppler = 0D0
        this%par%coeff_len2tottau = 0D0
        this%par%coeff_len2abstau = 0D0
        this%par%coeff_len2scatau = 0D0
        this%par%gaussVal         = 0D0
      end if
      ! Only for testing purpose.  Otherwise the above lines should be used!
      ! 0.
      ! this%par%n      = this%d(3)
      ! this%par%eps    = this%d(4)
      ! this%par%sc     = 0D0
      ! this%par%ab     = this%d(6)
      ! this%par%tot    = this%par%sc + this%par%ab
      ! this%par%alb    = 0D0
      ! this%par%g      = 1D0
      ! this%par%coeff_len2tottau = 0D0
      ! this%par%coeff_len2abstau = 0D0
      ! this%par%coeff_len2scatau = 0D0
      ! 1.
      ! this%par%n      = this%d(3)
      ! this%par%eps    = this%d(4)
      ! this%par%sc     = 0D0
      ! this%par%ab     = this%d(6)
      ! this%par%tot    = this%par%sc + this%par%ab
      ! this%par%alb    = 0D0
      ! this%par%g      = 1D0
      ! this%par%coeff_len2tottau = 0.01D0
      ! this%par%coeff_len2abstau = 0.01D0
      ! this%par%coeff_len2scatau = 0D0
      ! 2.
      ! this%par%n      = this%d(3)
      ! this%par%eps    = this%d(4)
      ! this%par%sc     = 0D0
      ! this%par%ab     = this%d(6)
      ! this%par%tot    = this%par%sc + this%par%ab
      ! this%par%alb    = 0D0
      ! this%par%g      = 1D0
      ! this%par%coeff_len2tottau = 3D-5 * sqrt(this%x1*this%x1 + this%y1*this%y1)
      ! this%par%coeff_len2abstau = this%par%coeff_len2tottau
      ! this%par%coeff_len2scatau = 0D0
      !!!!!!!!!!!!!!
      return
    else
      do i=1, this%nChildren
        call make_cell_leaves(this%children(i), idx)
      end do
    end if
  end subroutine make_cell_leaves


  subroutine make_cell_neighbors
    implicit none
    integer i, j
    type(cell), pointer :: this, maybe
    type(cell_neighbors), pointer :: neib
    integer :: nMaxIn=0, nMaxOut=0, nMaxUp=0, nMaxDown=0
    integer :: iMaxIn, iMaxOut, iMaxUp, iMaxDown
    do i=1, n_cell_leaves
      this => cell_leaves(i)%p
      allocate(this%nei)
      neib => this%nei
      if (((this%x2-this%x1) .LE. very_small_len) .OR. &
          ((this%y2-this%y1) .LE. very_small_len)) then
        write(*,*) 'Your code has a problem!'
      end if
      if (this%x1 .LE. very_small_len) then
        neib%nIn = neib%nIn + 1
        neib%iIn(neib%nIn) = i
      end if
      if (this%y1 .LE. very_small_len) then
        neib%nDown = neib%nDown + 1
        neib%iDown(neib%nDown) = i
      end if
      do j=1, n_cell_leaves
        maybe => cell_leaves(j)%p
        ! Neighborhood in different directions are not exclusive.
        if (abs(this%x1 - maybe%x2) .LE. very_small_len) then
          if ((this%y1 .LE. maybe%y2+very_small_len) .AND. &
              (this%y2 .GE. maybe%y1-very_small_len)) then
            neib%nIn = neib%nIn + 1
            neib%iIn(neib%nIn) = j
            if (maybe%y1 .LE. this%y1) then
              neib%nDown = neib%nDown + 1
              neib%iDown(neib%nDown) = j
            end if
            if (maybe%y2 .GE. this%y2) then
              neib%nUp = neib%nUp + 1
              neib%iUp(neib%nUp) = j
            end if
          end if
        end if
        if (abs(this%x2 - maybe%x1) .LE. very_small_len) then
          if ((this%y1 .LE. maybe%y2+very_small_len) .AND. &
              (this%y2 .GE. maybe%y1-very_small_len)) then
            neib%nOut = neib%nOut + 1
            neib%iOut(neib%nOut) = j
            if (maybe%y1 .LE. this%y1) then
              neib%nDown = neib%nDown + 1
              neib%iDown(neib%nDown) = j
            end if
            if (maybe%y2 .GE. this%y2) then
              neib%nUp = neib%nUp + 1
              neib%iUp(neib%nUp) = j
            end if
          end if
        end if
        if (abs(this%y1 - maybe%y2) .LE. very_small_len) then
          if ((this%x1 .LE. maybe%x2+very_small_len) .AND. &
              (this%x2 .GE. maybe%x1-very_small_len)) then
            neib%nDown = neib%nDown + 1
            neib%iDown(neib%nDown) = j
            if (maybe%x1 .LE. this%x1) then
              neib%nIn = neib%nIn + 1
              neib%iIn(neib%nIn) = j
            end if
            if (maybe%y2 .GE. this%y2) then
              neib%nOut = neib%nOut + 1
              neib%iOut(neib%nOut) = j
            end if
          end if
        end if
        if (abs(this%y2 - maybe%y1) .LE. very_small_len) then
          if ((this%x1 .LE. maybe%x2+very_small_len) .AND. &
              (this%x2 .GE. maybe%x1-very_small_len)) then
            neib%nUp = neib%nUp + 1
            neib%iUp(neib%nUp) = j
            if (maybe%x1 .LE. this%x1) then
              neib%nIn = neib%nIn + 1
              neib%iIn(neib%nIn) = j
            end if
            if (maybe%y2 .GE. this%y2) then
              neib%nOut = neib%nOut + 1
              neib%iOut(neib%nOut) = j
            end if
          end if
        end if
      end do
      if (nMaxIn .LT. neib%nIn) then
        nMaxIn = neib%nIn
        iMaxIn = i
      end if
      if (nMaxOut .LT. neib%nOut) then
        nMaxOut = neib%nOut
        iMaxOut = i
      end if
      if (nMaxUp .LT. neib%nUp) then
        nMaxUp = neib%nUp
        iMaxUp = i
      end if
      if (nMaxDown .LT. neib%nDown) then
        nMaxDown = neib%nDown
        iMaxDown = i
      end if
    end do
    write(*,'(A32, I5, "@", I8)') 'Neighbor number: nMaxIn   =',  nMaxIn,  iMaxIn
    write(*,'(A32, I5, "@", I8)') 'Neighbor number: nMaxOut  =',  nMaxOut, iMaxOut
    write(*,'(A32, I5, "@", I8)') 'Neighbor number: nMaxUp   =',  nMaxUp,  iMaxUp
    write(*,'(A32, I5, "@", I8)') 'Neighbor number: nMaxDown =',  nMaxDown,iMaxDown
  end subroutine make_cell_neighbors


  subroutine find_out_boundary_points
    implicit none
    logical, dimension(4) :: hasPointsInQuadrants
    logical isBoundary
    integer :: i, j, nBoundary_maybe=0
    integer, dimension(2048) :: idxBoundary_maybe
    do i=1, len_inputs
      isBoundary = .TRUE.
      hasPointsInQuadrants = .FALSE.
      do j=1, len_inputs
        if      ((data_inputs(2, j) .GT. data_inputs(2, i)) .AND. &
                 (data_inputs(3, j) .GT. data_inputs(3, i))) then
          hasPointsInQuadrants(1) = .TRUE.
        else if ((data_inputs(2, j) .LT. data_inputs(2, i)) .AND. &
                 (data_inputs(3, j) .GT. data_inputs(3, i))) then
          hasPointsInQuadrants(2) = .TRUE.
        else if ((data_inputs(2, j) .LT. data_inputs(2, i)) .AND. &
                 (data_inputs(3, j) .LT. data_inputs(3, i))) then
          hasPointsInQuadrants(3) = .TRUE.
        else if ((data_inputs(2, j) .GT. data_inputs(2, i)) .AND. &
                 (data_inputs(3, j) .LT. data_inputs(3, i))) then
          hasPointsInQuadrants(4) = .TRUE.
        end if
        if (hasPointsInQuadrants(1) .AND. hasPointsInQuadrants(2) .AND. &
            hasPointsInQuadrants(3) .AND. hasPointsInQuadrants(4)) then
          isBoundary = .FALSE.
          exit
        end if
      end do
      if (isBoundary) then
        nBoundary_maybe = nBoundary_maybe + 1
        idxBoundary_maybe(nBoundary_maybe) = i
      end if
    end do
    do i=1, nBoundary_maybe
      write(*,*) i, data_inputs(2, idxBoundary_maybe(i)), data_inputs(3, idxBoundary_maybe(i))
    end do
  end subroutine find_out_boundary_points

end module cell_class
