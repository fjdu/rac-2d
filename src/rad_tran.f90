
module rad_tran_mc

  use trivials
  use generic_type
  use cell_class
  use global_var
  use phy_const
  use physical_par

  double precision photon_dir_w_min, photon_dir_w_max

  integer MC_n_photons, MC_step_in_path_max
  double precision, dimension(2) :: photon_elevation_range
  double precision photon_weight_threshold
  integer source_SED_num
  character(len=128) source_SED_file
  double precision, dimension(:,:), allocatable :: source_SED
  double precision, dimension(:,:), allocatable :: source_SED_freq
  integer, dimension(:), allocatable :: source_SED_nPhoton
  double precision, dimension(:), allocatable :: out_freq_intervals
  integer out_freq_num
  double precision out_freq_low, out_freq_high
  double precision out_wavelength_low, out_wavelength_high
  namelist /MonteCarlo_RadTran_par/ &
    MC_n_photons, MC_step_in_path_max, photon_elevation_range, &
    photon_weight_threshold, &
    source_SED_file, &
    out_freq_num, out_freq_low, out_freq_high, &
    out_wavelength_low, out_wavelength_high

  double precision, dimension(:), allocatable :: photon_freq_list

  type :: ray_seg
    type(type_point) :: bg, ed
  end type ray_seg

  type :: photon_vec
    type(type_point) :: pt
    type(type_direction_cartesian) :: dir
    type(cell), pointer :: cell
    double precision rho, z
    double precision freq
    double precision weight
    logical free, absorbed
    integer In_Out_Up_Down
    double precision :: vel_projected = 0D0
  end type photon_vec

  type :: receptor
    type(type_point) pt
    type(type_direction_cartesian) dir
    double precision radius
  end type receptor

  type :: rad_grid
    double precision x1, x2, y1, y2, z1, z2
    double precision dx, dy, dz
    integer nx, ny, nz
    integer, dimension(:,:,:,:), allocatable :: ph_count_dir
    integer, dimension(:,:,:), allocatable :: ph_count
    double precision, dimension(:,:,:), allocatable :: J
  end type rad_grid

  contains

  subroutine initialize_rad_tran_mc
    implicit none
    photon_dir_w_min = sin(photon_elevation_range(1)*phy_Deg2Rad)
    photon_dir_w_max = sin(photon_elevation_range(2)*phy_Deg2Rad)
    call read_source_SED()
    call make_out_freq_intervals()
    call allocate_cell_J_F()
  end subroutine initialize_rad_tran_mc


  type(photon_vec) function get_photon_initial(n, amr)
    implicit none
    integer n
    double precision, dimension(2) :: rnd
    type(cell), pointer, intent(in) :: amr
    double precision tmp
    ! Isotropic stellar field
    get_photon_initial%pt%x = 0D0
    get_photon_initial%pt%y = 0D0
    get_photon_initial%pt%z = 0D0
    get_photon_initial%free = .FALSE.
    get_photon_initial%absorbed = .FALSE.
    !get_photon_initial%dir = rand_pt_sphere_uniform()
    get_photon_initial%dir = rand_pt_sphere_uniform_range(photon_dir_w_min, photon_dir_w_max)
    get_photon_initial%rho = sqrt( &
      get_photon_initial%pt%x*get_photon_initial%pt%x + &
      get_photon_initial%pt%y*get_photon_initial%pt%y)
    get_photon_initial%z = abs(get_photon_initial%pt%z)
    get_photon_initial%weight = 1D0
    call random_number(tmp)
    get_photon_initial%freq = &
      source_SED_freq(1, n) * tmp + source_SED_freq(2, n) * (1D0-tmp)
    !
    if (is_inside_cell_(get_photon_initial, amr)) then
      call locate_photon_cell(get_photon_initial, amr)
      get_photon_initial%cell%rad%counts = get_photon_initial%cell%rad%counts + 1
    else
      get_photon_initial%cell => null()
    end if
    !
  end function get_photon_initial


  logical function is_inside_cell_(ph, amr)
    implicit none
    type(photon_vec) ph
    type(cell), pointer, intent(in) :: amr
    if ((ph%rho .LT. amr%x1) .OR. (ph%rho .GT. amr%x2) .OR. &
        (ph%z .LT. amr%y1) .OR. (ph%z .GT. amr%y2)) then
      is_inside_cell_ = .FALSE.
    else
      is_inside_cell_ = .TRUE.
    end if
  end function is_inside_cell_


  recursive subroutine locate_photon_cell(ph, amr)
    implicit none
    type(photon_vec) ph
    type(cell), pointer, intent(in) :: amr
    type(cell), pointer :: tmp
    integer i
    ! When calling this function, you must make sure that the photon is inside amr.
    if (amr%nOffspr .EQ. 0) then
      ph%cell => amr
      return
    end if
    do i=1, amr%nChildren
      tmp => amr%children(i)
      ! Be cautious here!
      if ((ph%rho .GE. tmp%x1) .AND. (ph%rho .LT. tmp%x2) .AND. &
          (ph%z .GE. tmp%y1)   .AND. (ph%z .LT. tmp%y2)) then
        call locate_photon_cell(ph, tmp)
        return
      end if
    end do
  end subroutine locate_photon_cell


  subroutine update_photon_cellInfo(ph, found)
    implicit none
    type(photon_vec) ph
    type(cell), pointer :: this
    integer i
    logical found
    type(cell_neighbors), pointer :: nei
    found = .FALSE.
    nei => ph%cell%nei
    select case(ph%In_Out_Up_Down)
      case (1) ! Inward
        do i=1, nei%nIn
          this => cell_leaves(nei%iIn(i))%p
          if ((ph%rho .GE. this%x1) .AND. (ph%rho .LT. this%x2) .AND. &
              (ph%z .GE. this%y1)   .AND. (ph%z .LT. this%y2)) then
            ph%cell => this
            ph%cell%rad%counts = ph%cell%rad%counts + 1
            found = .TRUE.
            return
          end if
        end do
      case (2) ! Outward
        do i=1, nei%nOut
          this => cell_leaves(nei%iOut(i))%p
          if ((ph%rho .GE. this%x1) .AND. (ph%rho .LT. this%x2) .AND. &
              (ph%z .GE. this%y1)   .AND. (ph%z .LT. this%y2)) then
            ph%cell => this
            ph%cell%rad%counts = ph%cell%rad%counts + 1
            found = .TRUE.
            return
          end if
        end do
      case (3) ! Upward
        do i=1, nei%nUp
          this => cell_leaves(nei%iUp(i))%p
          if ((ph%rho .GE. this%x1) .AND. (ph%rho .LT. this%x2) .AND. &
              (ph%z .GE. this%y1)   .AND. (ph%z .LT. this%y2)) then
            ph%cell => this
            ph%cell%rad%counts = ph%cell%rad%counts + 1
            found = .TRUE.
            return
          end if
        end do
      case (4) ! Downward
        do i=1, nei%nDown
          this => cell_leaves(nei%iDown(i))%p
          if ((ph%rho .GE. this%x1) .AND. (ph%rho .LT. this%x2) .AND. &
              (ph%z .GE. this%y1)   .AND. (ph%z .LT. this%y2)) then
            ph%cell => this
            ph%cell%rad%counts = ph%cell%rad%counts + 1
            found = .TRUE.
            return
          end if
        end do
    end select
  end subroutine update_photon_cellInfo


  function get_next_photon(this_photon) result(next_ph)
    use global_var
    implicit none
    type(photon_vec) this_photon
    type(photon_vec) next_ph
    double precision next_tau, tau_sca_tot, tau_cell_abs
    doube precision sca_Lyman_alpha
    double precision length, length_weighted
    logical found_end, found_cell
    double precision tmp
    integer, parameter :: n_max_cell_traverse = 16384
    integer i
    type(type_direction_cartesian) d_rel
    !
    call random_number(tmp)
    next_tau = -log(tmp)
    next_ph = this_photon ! Copy the basic information.
    found_end = .FALSE.
    !
    do i=1, n_max_cell_traverse
      call get_tau_in_this_cell(next_ph, length)
      call calcVoigt(next_ph, sca_Lyman_alpha)
      tau_sca_dust = length * ph%cell%par%coeff_len2scatau
      tau_sca_Lyman_alpha = sca_Lyman_alpha * length
      tau_sca_tot = tau_sca_dust + tau_sca_Lyman_alpha
      tau_abs = length * ph%cell%par%coeff_len2abstau
      if (tau_sca_tot .GE. next_tau) then ! Scatter within the current cell.
        length       = length *       (next_tau / tau_sca_tot)
        tau_cell_abs = tau_cell_abs * (next_tau / tau_sca_tot)
        d_rel = rand_pt_HenyeyGreenstein(next_ph%cell%par%g)
        next_ph%dir = get_anisotropic_dir(next_ph%dir, d_rel)
        found_end = .TRUE.
      else ! Will go to the next cell
        next_tau = next_tau - tau_sca_tot
        length = length + very_small_len
      end if
      next_ph%pt%x = next_ph%pt%x + this_photon%dir%u * length
      next_ph%pt%y = next_ph%pt%y + this_photon%dir%v * length
      next_ph%pt%z = next_ph%pt%z + this_photon%dir%w * length
      next_ph%rho = sqrt(next_ph%pt%x*next_ph%pt%x + next_ph%pt%y*next_ph%pt%y)
      next_ph%z = abs(next_ph%pt%z)
      call get_cell_vel_where_photon(next_ph)
      tmp = exp(-tau_cell_abs)
      if (tau_cell_abs .LE. 1D-2) then
        length_weighted = length * next_ph%weight * &
            (1D0 - 0.5D0*tau_cell_abs*(1D0 - tau_cell_abs/3D0))
      else
        length_weighted = length * next_ph%weight * ((1D0 - tmp) / tau_cell_abs)
      end if
      next_ph%cell%rad%J = next_ph%cell%rad%J + length_weighted
      next_ph%weight = next_ph%weight * tmp
      if (next_ph%weight .LT. photon_weight_threshold) then
        next_ph%absorbed = .TRUE.
        return
      end if
      if (found_end) then
        if (.NOT. is_inside_cell_(next_ph, next_ph%cell)) then
          ! This is to take into account the rounding error.
          !call update_photon_cellInfo(next_ph, next_ph%cell%parent, found_cell)
          call update_photon_cellInfo(next_ph, found_cell)
          if (.NOT. found_cell) next_ph%free = .TRUE.
        end if
        return
      else
        if (.NOT. associated(next_ph%cell%parent)) then
          ! This should never happen.
          next_ph%free = .TRUE.
          write(*,*) 'Already at the top-level cell!'
          return
        end if
        call update_photon_cellInfo(next_ph, found_cell)
        if (.NOT. found_cell) then
          next_ph%free = .TRUE.
          return
        end if
      end if
    end do
    if (.NOT. found_end) write(*, '(A, I6, 7ES12.2/)') &
      '! Warning: end point not found!', &
      i, next_ph%pt%x, next_ph%pt%y, next_ph%pt%z, &
      next_ph%cell%x1, next_ph%cell%x2, next_ph%cell%y1, next_ph%cell%y2
  end function get_next_photon


  subroutine get_tau_in_this_cell(ph, length)
    ! Calculate the path length and tau from the location of a photon to the
    !   boundary of the cell that this photon resides in.
    use phy_const
    implicit none
    type(photon_vec), intent(inout) :: ph
    double precision lam1, lam2, length
    double precision A, B, C1, C2, D1, D2
    double precision, dimension(4) :: x
    double precision, parameter :: very_big_num  = huge(0D0)
    double precision, parameter :: very_tiny_num = tiny(0D0)
    integer i, itmp
    A = ph%pt%x*ph%dir%u + ph%pt%y*ph%dir%v
    B = ph%dir%u*ph%dir%u + ph%dir%v*ph%dir%v
    C1 = (ph%rho - ph%cell%x1) * (ph%rho + ph%cell%x1)
    C2 = (ph%rho - ph%cell%x2) * (ph%rho + ph%cell%x2)
    D1 = A*A - B * C1
    D2 = A*A - B * C2
    x = very_big_num
    if (ph%rho .GT. 0D0) then
      if (D1 .GE. 0D0) then
        x(1) = (-A - sqrt(D1)) / (B + very_tiny_num)
        x(2) = (-A + sqrt(D1)) / (B + very_tiny_num)
      end if
    end if
    if (D2 .GE. 0D0) then
      x(3) = (-A - sqrt(D2)) / (B + very_tiny_num)
      x(4) = (-A + sqrt(D2)) / (B + very_tiny_num)
    end if
    lam1 = very_big_num
    do i=1,4
      if (x(i) .GE. 0D0) then
        if (x(i) .LT. lam1) then 
          lam1 = x(i)
          ph%In_Out_Up_Down = floor(dble(i+1)/2D0) ! 1,2,3,4 -> 1,1,2,2; 1 = in; 2 = out
        end if
      end if
    end do
    if (ph%dir%w .GT. 0D0) then ! Upward
      if (ph%pt%z .GE. 0D0) then
        lam2 = (ph%cell%y2 - ph%z) / abs(ph%dir%w)
        itmp = 3
      else
        lam2 = (ph%z - ph%cell%y1) / abs(ph%dir%w)
        itmp = 4
      end if
    else if (ph%dir%w .LT. 0D0) then! Downward
      if (ph%pt%z .GT. 0D0) then
        lam2 = (ph%z - ph%cell%y1) / abs(ph%dir%w)
        itmp = 4
      else
        lam2 = (ph%cell%y2 - ph%z) / abs(ph%dir%w)
        itmp = 3
      end if
    else
      lam2 = very_big_num
    end if
    if (lam1 .LE. lam2) then
      length = lam1
    else
      length = lam2
      ph%In_Out_Up_Down = itmp
    end if
  end subroutine get_tau_in_this_cell


  subroutine calcVoigt(this, sca_Lyman_alpha)
    implicit none
    type(cell), pointer :: this
    double precision sca_Lyman_alpha
    double precision dv, nu_1
    this%par%Voigt(1) = 0D0
    sca_Lyman_alpha = 0D0
    dv = this%par%v_thermal(2) - this%par%v_thermal(1)
    do i=2, this%par%n_vel_intvl
      nu_1 = ph%freq - Lyman_alpha_freq * &
        (1D0 + (ph%vel_projected + this%par%v_thermal(i)) * / phy_SpeedOfLight)
      this%par%Voigt(i) = this%par%Voigt(i-1) + &
        dv * this%par%gaussVal(i) * (this%n*1D6) * &
          Lyman_alpha_crosssec_peak * &
          Lyman_alpha_Lorentz_width * Lyman_alpha_Lorentz_width / &
          (nu_1*nu_1 + Lyman_alpha_Lorentz_width * Lyman_alpha_Lorentz_width)
    end do
    sca_Lyman_alpha = this%par%Voigt(this%par%n_vel_intvl)
  end subroutine calcVoigt


  subroutine get_cell_vel_where_photon(ph)
    implicit none
    type(photon_vec) ph
    double precision tmp
    double precision, parameter :: const_small_num = 1D-10
    tmp = ph%pt%x * ph%pt%x + ph%pt%y * ph%pt%y
    if (tmp .GT. const_small_num) then
      ph%vel_projected = ph%cell%velocity * &
        (-ph%pt%y * ph%dir%u + ph%pt%x * ph%dir%v) / sqrt(tmp)
    else
      ph%vel_projected = 0D0
    end if
  end subroutine get_cell_vel_where_photon


  subroutine print_traverse(fU)
    implicit none
    type(cell), pointer :: this
    integer, intent(in) :: fU
    integer i
    double precision coeff_counts2flux, star_photon_num_freq, solid_angle, freq
    freq = phy_SpeedOfLight / (1216E-10)
    solid_angle = 2D0*phy_Pi * (photon_dir_w_max - photon_dir_w_min)
    star_photon_num_freq = &
      (star_luminosity_in_Lsun * phy_Lsun * phy_erg2joule) &
        * (ratio_uv2total * ratio_lyman2uv) & !!! Something is still lacking here!!!
        / (phy_hPlanck*freq)
    coeff_counts2flux = &
      star_photon_num_freq &
        * (solid_angle/(4D0*phy_Pi)) &
        / MC_n_photons
    write(*, '(A, ES16.5)') 'Photon frequency (Hz) = ', freq
    write(*, '(A, ES16.5)') 'Lyman alpha photons emitted by the star = ', star_photon_num_freq
    write(fU, '("! Be careful with the convention of J.")')
    write(fU, '("!", A5, A4, 5A12, A9, 8A12)') &
        '0', '  1',  ' 2',    ' 3',    ' 4',    ' 5',    &
        '6',        '     7',        '8',        '   9',              &
        ' 10',       ' 11',       ' 12',        ' 13',        &
        '14',       '   15'
    write(fU, '("!", A5, A4, 5A12, A9, 8A12)') &
        'i', 'ord',  'x1',    'y1',    'x2',    'y2',    &
        'n',        'counts',        'J',        'n_ph',              &
        'sca',       'abs',       'alb',        'eps',        &
        'g',        'coeff'
    do i=1, n_cell_leaves
      this => cell_leaves(i)%p
      this%rad%J = this%rad%J * coeff_counts2flux &
        / (4D0*phy_Pi * this%par%vol * (phy_AU2cm*phy_AU2cm))
      this%rad%photon_den = this%rad%J * (4D0*phy_Pi / (phy_SpeedOfLight * phy_m2cm))
      if (isnan(this%d(3))) write(*,*) 'Cell value is not set!', &
        i, this%order, this%x1, this%y1, this%x2, this%y2
      write(fU, '(I6, I4, 5ES12.3E3, I9, 8ES12.3E3)') &
        i, this%order, this%x1, this%y1, this%x2, this%y2, &
        this%par%n, this%rad%counts, this%rad%J, this%rad%photon_den, &
        this%par%sc, this%par%ab, this%par%alb, this%par%eps, &
        this%par%g, this%par%coeff_len2tottau
    end do
  end subroutine print_traverse


  double precision function get_Euclidean_dist(pt1, pt2)
    implicit none
    type(type_point), intent(in) :: pt1, pt2
    get_Euclidean_dist = sqrt( &
      (pt1%x - pt2%x) * (pt1%x - pt2%x) + &
      (pt1%y - pt2%y) * (pt1%y - pt2%y) + &
      (pt1%z - pt2%z) * (pt1%z - pt2%z))
  end function get_Euclidean_dist


  type(type_direction_cartesian) function get_dir_pt(pt1, pt2)
    implicit none
    type(type_point) pt1, pt2
    double precision dist
    !dist = get_Euclidean_dist(pt1, pt2)
    dist = sqrt( &
      (pt1%x - pt2%x) * (pt1%x - pt2%x) + &
      (pt1%y - pt2%y) * (pt1%y - pt2%y) + &
      (pt1%z - pt2%z) * (pt1%z - pt2%z))
    get_dir_pt%u = (pt2%x - pt1%x) / dist
    get_dir_pt%v = (pt2%y - pt1%y) / dist
    get_dir_pt%w = (pt2%z - pt1%z) / dist
  end function get_dir_pt


  type(type_direction_cartesian) function rand_pt_sphere_uniform()
    implicit none
    double precision x1, x2, s
    double precision, dimension(2) :: rnd
    do
      call random_number(rnd)
      x1 = rnd(1)*2D0 - 1D0
      x2 = rnd(2)*2D0 - 1D0
      s = x1*x1 + x2*x2
      if (s .LT. 1D0) then
        rand_pt_sphere_uniform%u = 2D0*x1*sqrt(1D0 - s)
        rand_pt_sphere_uniform%v = 2D0*x2*sqrt(1D0 - s)
        rand_pt_sphere_uniform%w = 1D0 - 2D0 * s
        exit
      end if
    end do
  end function rand_pt_sphere_uniform


  type(type_direction_cartesian) function rand_pt_sphere_uniform_range(minw, maxw)
    implicit none
    double precision x1, x2, s, minw, maxw
    double precision, dimension(2) :: rnd
    do
      call random_number(rnd)
      x1 = rnd(1)*2D0 - 1D0
      x2 = rnd(2)*2D0 - 1D0
      s = x1*x1 + x2*x2
      if (s .LE. 1D0) then
        rand_pt_sphere_uniform_range%w = 1D0 - 2D0 * s
        if ((rand_pt_sphere_uniform_range%w .LE. maxw) .AND. &
            (rand_pt_sphere_uniform_range%w .GE. minw)) then
          rand_pt_sphere_uniform_range%u = 2D0*x1*sqrt(1D0 - s)
          rand_pt_sphere_uniform_range%v = 2D0*x2*sqrt(1D0 - s)
          exit
      end if
      end if
    end do
  end function rand_pt_sphere_uniform_range


  type(type_direction_cartesian) function rand_pt_HenyeyGreenstein(g)
    use phy_const
    implicit none
    double precision, intent(in) :: g
    double precision t, phi, costheta, sintheta
    double precision, dimension(2) :: p
    call random_number(p)
    t = (1D0 - g*g) / (1D0 + g * (2D0 * p(1) - 1D0))
    costheta = 0.5D0 / g * &
               (1D0 + g * g - t * t)
    sintheta = sqrt(1D0 - costheta * costheta)
    phi = (2D0 * phy_Pi) * p(2)
    rand_pt_HenyeyGreenstein%u = sintheta * cos(phi)
    rand_pt_HenyeyGreenstein%v = sintheta * sin(phi)
    rand_pt_HenyeyGreenstein%w = costheta
  end function rand_pt_HenyeyGreenstein


  type(type_direction_cartesian) function get_anisotropic_dir(dir0, dir_rel)
    type(type_direction_cartesian) dir0, dir_rel
    type(type_sphere_coor_quat) d
    d = cartesian_to_spherical_quat(dir0)
    get_anisotropic_dir = rot_around_Y(dir_rel, d%costheta, d%sintheta)
    get_anisotropic_dir = rot_around_Z(get_anisotropic_dir, &
      d%cosphi, d%sinphi)
  end function get_anisotropic_dir


  type(type_sphere_coor_quat) function cartesian_to_spherical_quat(d)
  ! Convert (u, v, w) into (costheta, sintheta, cosphi, sinphi)
    implicit none
    type(type_direction_cartesian) d
    cartesian_to_spherical_quat%costheta = d%w
    cartesian_to_spherical_quat%sintheta = sqrt(1D0 - d%w * d%w)
    if (cartesian_to_spherical_quat%sintheta .GT. 0D0) then
      cartesian_to_spherical_quat%cosphi = d%u / cartesian_to_spherical_quat%sintheta
      cartesian_to_spherical_quat%sinphi = d%v / cartesian_to_spherical_quat%sintheta
    else
      cartesian_to_spherical_quat%cosphi = 0D0
      cartesian_to_spherical_quat%sinphi = 1D0
    end if
  end function cartesian_to_spherical_quat


  type(type_direction_cartesian) function rot_around_X(dir, cosa, sina)
  ! Follow the right rand convention
    implicit none
    type(type_direction_cartesian) dir
    double precision cosa, sina
    rot_around_X%u = dir%u
    rot_around_X%v = dir%v * cosa - dir%w * sina
    rot_around_X%w = dir%w * cosa + dir%v * sina
  end function rot_around_X


  type(type_direction_cartesian) function rot_around_Y(dir, cosa, sina)
  ! Follow the right rand convention
    implicit none
    type(type_direction_cartesian) dir
    double precision cosa, sina
    rot_around_Y%u = dir%u * cosa + dir%w * sina
    rot_around_Y%v = dir%v
    rot_around_Y%w = dir%w * cosa - dir%u * sina
  end function rot_around_Y


  type(type_direction_cartesian) function rot_around_Z(dir, cosa, sina)
  ! Follow the right rand convention
    implicit none
    type(type_direction_cartesian) dir
    double precision cosa, sina
    rot_around_Z%u = dir%u * cosa - dir%v * sina
    rot_around_Z%v = dir%v * cosa + dir%u * sina
    rot_around_Z%w = dir%w
  end function rot_around_Z


  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
  end subroutine


  subroutine read_source_SED()
    implicit none
    character(len=64) strtmp, strtmp1
    character :: comment_char='!'
    integer fU, ios, i
    double precision source_SED_tot
    if (.NOT. getFileUnit(fU)) then
      write(*,*) 'Cannot get a free file unit!'
      stop
    end if
    source_SED_num = GetFileLen_comment(source_SED_file, comment_char)
    allocate(source_SED(3, source_SED_num), &
             source_SED_freq(2, source_SED_num), &
             source_SED_nPhoton(source_SED_num))
    call openFileSequentialRead(fU, source_SED_file, 999)
    source_SED = 0D0
    source_SED_freq = 0D0
    source_SED_num = 0
    source_SED_tot = 0D0
    do
      read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
      if (ios .LT. 0) exit
      strtmp1 = adjustl(strtmp)
      if (IsDigitChar(strtmp1(1:1))) then
        source_SED_num = source_SED_num + 1
        read(strtmp, '(3F8.2)') source_SED(1,source_SED_num), &
          source_SED(2,source_SED_num), source_SED(3,source_SED_num)
        source_SED_freq(1:2, source_SED_num) = &
          phy_SpeedOfLight/source_SED(1:2, source_SED_num) * 1D10
        source_SED_tot = source_SED_tot + source_SED(3,source_SED_num)
      end if
    end do
    close(fU)
    do i=1, source_SED_num
      source_SED_nPhoton(i) = int(dble(MC_n_photons) &
        * source_SED(3, i)/source_SED_tot)
    end do
  end subroutine read_source_SED


  subroutine make_out_freq_intervals()
    implicit none
    double precision del_freq
    integer i
    if ((out_freq_low .LE. 1D-1) .OR. (out_freq_high .LE. 1D-1)) then
      if ((out_wavelength_low  .LE. 1D-1) .OR. &
          (out_wavelength_high .LE. 1D-1)) then
        out_wavelength_low  = source_SED(1,1)
        out_wavelength_high = source_SED(2, source_SED_num)
        stop
      end if
      out_freq_low  = phy_SpeedOfLight/out_wavelength_high*1D10
      out_freq_high = phy_SpeedOfLight/out_wavelength_low*1D10
    end if
    if (out_freq_num .EQ. 0) then
      out_freq_num = source_SED_num
    end if
    allocate(out_freq_intervals(out_freq_num+2))
    del_freq = (out_freq_high - out_freq_low) / dble(out_freq_num-1)
    out_freq_intervals(1) = 0D0
    out_freq_intervals(out_freq_num+2) = huge(0D0)
    out_freq_intervals(2) = out_freq_low
    do i=3, out_freq_num+1
      out_freq_intervals(i) = out_freq_intervals(i-1) + del_freq
    end do
  end subroutine make_out_freq_intervals


  subroutine allocate_cell_J_F()
    implicit none
    integer i
    do i=1, n_cell_leaves
      allocate(cell_leaves(i)%p%rad%Jnu(out_freq_num+2), &
               cell_leaves(i)%p%rad%Fnu(out_freq_num+2))
      cell_leaves(i)%p%rad%Jnu = 0D0
      cell_leaves(i)%p%rad%Fnu = 0D0
    end do
  end subroutine allocate_cell_J_F


  subroutine getGaussRnd2(rnd, cen, sigma)
    implicit none
    double precision, dimension(2) :: rnd
    double precision, dimension(2) :: uv
    double precision cen, sigma
    double precision tmp
    call random_number(uv)
    tmp = sqrt(-2D0*log(uv(1))) * sigma
    rnd(1) = tmp * cos(phy_2Pi*uv(2)) + cen
    rnd(2) = tmp * sin(phy_2Pi*uv(2)) + cen
  end subroutine getGaussRnd2


  !subroutine pseudoVoigtProfile(nu, nu0, v0, a, vth)
  !  implicit none
  !end subroutine pseudoVoigtProfile


end module rad_tran_mc
