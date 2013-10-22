module montecarlo

use phy_const
use trivials
use grid
use data_struct

implicit none


type(type_optical_property) dust_0, HI_0, water_0
type(type_stellar_spectrum) star_0
type(type_LUT_Tdust) lut_0
type(type_montecarlo_config) mc_conf
type(type_distribution_table) p4lam
type(type_global_material_collection) gl_coll_0

double precision, dimension(2), parameter :: lam_range_UV = (/9D2, 3D3/)
double precision, dimension(2), parameter :: lam_range_LyA = (/1210D0, 1220D0/)
double precision, dimension(2), parameter :: lam_range_NIR = (/8D3, 5D4/)
double precision, dimension(2), parameter :: lam_range_MIR = (/5D4, 3D5/)
double precision, dimension(2), parameter :: lam_range_FIR = (/3D5, 2D6/)

double precision, parameter :: refine_UV = 0.01D0
double precision, parameter :: refine_LyA = 0.001D0


namelist /montecarlo_configure/ mc_conf


contains


subroutine save_photon(ph, mc)
  type(type_photon_packet), intent(in) :: ph
  type(type_montecarlo_config), intent(in) :: mc
  if (mc%savephoton) then
    write(mc%fU, '(9ES19.10, I10)') & 
        ph%ray%x, ph%ray%y, ph%ray%z, &
        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
        ph%lam, ph%en, ph%wei, ph%e_count
  end if
end subroutine save_photon



subroutine montecarlo_prep
  integer n, n0
  double precision lam_max, lam_start
  double precision, dimension(:), allocatable :: ltmp, vtmp
  type(type_stellar_spectrum) stmp
  double precision, parameter :: T_Lya = 300D0
  !
  mc_conf%minw = sin(mc_conf%min_ang*phy_Deg2Rad)
  mc_conf%maxw = sin(mc_conf%max_ang*phy_Deg2Rad)
  mc_conf%maxw = get_surf_max_angle()
  write(*,'(/A, 2ES12.4)') 'minw,maxw = ', mc_conf%minw, mc_conf%maxw
  !
  call load_dust_data( &
    combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_dust), dust_0)
  !
  call load_H2O_ab_crosssection( &
    combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_water), water_0)
  !
  call make_H_Lya(T_Lya, HI_0)
  !
  call align_optical_data
  !
  call make_LUT_Tdust(dust_0, lut_0)
  !
  star_0%mass = mc_conf%star_mass
  star_0%radius = mc_conf%star_radius
  star_0%T = mc_conf%star_temperature
  !
  lam_max = min(1D6, dust_0%lam(dust_0%n)) ! in angstrom
  if (mc_conf%use_blackbody_star) then
    call make_stellar_spectrum(dust_0%lam(1), &
      lam_max, 10000, star_0)
  else
    !
    call load_stellar_spectrum( &
      combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_star), star_0)
    if (star_0%lam(star_0%n) .lt. lam_max) then
      ! Fill the remainder with blackbody radiation
      stmp%mass = star_0%mass
      stmp%radius = star_0%radius
      stmp%T = star_0%T
      lam_start = 2D0 * star_0%lam(star_0%n) - star_0%lam(star_0%n - 1)
      n = max(10, int(lam_max / lam_start * 10))
      call make_stellar_spectrum(lam_start, &
        lam_max, n, stmp)
      !
      n0 = star_0%n
      allocate(ltmp(n0), vtmp(n0))
      ltmp = star_0%lam
      vtmp = star_0%vals
      !
      star_0%n = star_0%n + n
      deallocate(star_0%lam, star_0%vals)
      allocate(star_0%lam(star_0%n), star_0%vals(star_0%n))
      !
      star_0%lam(1:n0) = ltmp
      star_0%vals(1:n0) = vtmp
      star_0%lam((n0+1) : star_0%n) = stmp%lam
      star_0%vals((n0+1): star_0%n) = stmp%vals
      !
      deallocate(ltmp, vtmp, stmp%lam, stmp%vals)
    end if
    !
  end if
  !
  star_0%lumi = get_stellar_luminosity(star_0)
  star_0%lumi_UV = get_stellar_luminosity(star_0, lam_range_UV(1), lam_range_UV(2))
  write(*,'(A, ES16.6, A)') 'Stellar total luminosity: ', star_0%lumi, ' erg s-1.'
  write(*,'(A, ES16.6, A)') 'Stellar UV luminosity: ', star_0%lumi_UV, ' erg s-1.'
  !
  star_0%minw = mc_conf%minw
  star_0%maxw = mc_conf%maxw
  !
  star_0%vals = star_0%vals * (star_0%maxw - star_0%minw) / 2D0
  star_0%lumi = star_0%lumi * (star_0%maxw - star_0%minw) / 2D0
  star_0%lumi_UV = star_0%lumi_UV * (star_0%maxw - star_0%minw) / 2D0
  !
  mc_conf%eph = star_0%lumi / dble(mc_conf%nph)
  !
  write(*,'(A, ES16.6, A)') 'Stellar luminosity within (minw,maxw): ', star_0%lumi, ' erg s-1.'
  write(*,'(A, ES16.6, A//)') 'Lumi per photon: ', mc_conf%eph, ' erg s-1.'
  !
  call make_global_coll
  !
  p4lam%n = lut_0%m
  allocate(p4lam%pvals(0:p4lam%n))
  !
end subroutine montecarlo_prep



subroutine align_optical_data
  ! The imported optical data for dust, HI, and water are not aligned.
  ! So here I will align them to the same lambda vector.
  double precision, dimension(:), allocatable :: v, v1
  integer n, n1, n_using
  integer i
  !
  n1 = HI_0%n + water_0%n
  n = dust_0%n + HI_0%n + water_0%n
  allocate(v1(n1), v(n))
  call merge_vec(water_0%n, water_0%lam, HI_0%n, HI_0%lam, n1, v1, n_using)
  call merge_vec(dust_0%n, dust_0%lam, n1, v1, n, v, n_using)
  !
  call reasign_optical(dust_0, n, v)
  call reasign_optical(HI_0, n, v)
  call reasign_optical(water_0, n, v)
  !
  !do i=1, n
  !  !write(*,*) i, dust_0%lam(i), dust_0%ab(i), HI_0%ab(i), water_0%ab(i)
  !  !write(*,*) i, dust_0%lam(i), dust_0%sc(i), HI_0%sc(i), water_0%sc(i)
  !  !write(*,*) i, dust_0%lam(i), dust_0%g(i), HI_0%g(i), water_0%g(i)
  !end do
  deallocate(v1, v)
end subroutine align_optical_data



subroutine reasign_optical(op, n, v)
  type(type_optical_property), intent(inout) :: op
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: v
  double precision, allocatable, dimension(:) :: y
  !
  allocate(y(n))
  !
  call transfer_value(op%n, op%lam, op%ab, n, v, y)
  deallocate(op%ab)
  allocate(op%ab(n))
  op%ab = y
  !
  call transfer_value(op%n, op%lam, op%sc, n, v, y)
  deallocate(op%sc)
  allocate(op%sc(n))
  op%sc = y
  !
  call transfer_value(op%n, op%lam, op%g, n, v, y)
  deallocate(op%g)
  allocate(op%g(n))
  op%g = y
  !
  deallocate(op%lam)
  allocate(op%lam(n))
  op%lam = v
  !
  op%n = n
  !
  deallocate(y)
end subroutine reasign_optical



subroutine transfer_value(n1, x1, y1, n2, x2, y2)
  integer, intent(in) :: n1, n2
  double precision, dimension(n1), intent(in) :: x1, y1
  double precision, dimension(n2), intent(in) :: x2
  double precision, dimension(n2), intent(out) :: y2
  integer i, i0, j
  i0 = 1
  do i=1, n2
    if (x2(i) .lt. x1(1)) then
      y2(i) = 0D0
    else if (x2(i) .gt. x1(n1)) then
      y2(i) = 0D0
    else
      do j=i0, n1-1
        if ((x1(j) .le. x2(i)) .and. (x1(j+1) .ge. x2(i))) then
          y2(i) = y1(j) + (x2(i) - x1(j)) * (y1(j+1) - y1(j)) / (x1(j+1) - x1(j))
          i0 = j
          exit
        end if
      end do
    end if
  end do
end subroutine transfer_value


subroutine merge_vec(n1, v1, n2, v2, n, v, n_using)
  use qsort_c_module
  integer, intent(in) :: n1, n2, n
  double precision, dimension(n1), intent(in) :: v1
  double precision, dimension(n2), intent(in) :: v2
  double precision, dimension(n), intent(out) :: v
  integer, intent(out) :: n_using
  !
  v(1:n1) = v1
  v(n1+1:n1+n2) = v2
  v(n1+n2+1:n) = huge(0D0)
  call QsortC(v)
  n_using = n1 + n2
end subroutine merge_vec



subroutine make_global_coll
  gl_coll_0%ntype = 3
  allocate(gl_coll_0%list(gl_coll_0%ntype))
  gl_coll_0%list(1) = dust_0
  gl_coll_0%list(2) = HI_0
  gl_coll_0%list(3) = water_0
end subroutine make_global_coll




subroutine prep_local_optics(c, gl, dust)
  type(type_cell), intent(inout), pointer :: c
  type(type_global_material_collection), intent(in) :: gl
  type(type_optical_property), intent(in) :: dust
  if (.not. c%using) then
    return
  end if
  if (allocated(c%optical%X)) then
    deallocate(c%optical%X, &
               c%optical%flux, &
               c%optical%acc, &
               c%optical%summed, &
               c%optical%summed_ab, &
               c%optical%summed_sc, &
               c%optical%dir_wei)
  end if
  c%optical%ntype = gl%ntype*2
  c%optical%nlam = dust%n
  allocate(c%optical%X(c%optical%ntype), &
           c%optical%flux(c%optical%nlam), &
           c%optical%acc(c%optical%nlam, c%optical%ntype), &
           c%optical%summed(c%optical%nlam), &
           c%optical%summed_ab(c%optical%nlam), &
           c%optical%summed_sc(c%optical%nlam), &
           c%optical%dir_wei(c%optical%nlam))
end subroutine prep_local_optics



subroutine reset_local_optics(c)
  type(type_cell), intent(inout), pointer :: c
  integer i
  if (.not. c%using) then
    return
  end if
  c%optical%en_gain = 0D0
  c%optical%en_gain_abso = 0D0
  c%optical%en_prev = 0D0
  c%optical%ab_count = 0
  c%optical%cr_count = 0
  c%optical%kph = 0D0
  do i=1, c%optical%nlam
    c%optical%dir_wei(i)%u = 0D0
    c%optical%dir_wei(i)%v = 0D0
    c%optical%dir_wei(i)%w = 0D0
    c%optical%flux(i) = 0D0
  end do
  call update_local_coll(c)
  call make_local_material_collection(c%optical, gl_coll_0)
end subroutine reset_local_optics



subroutine update_local_coll(c)
  type(type_cell), intent(inout), pointer :: c
  if (c%using) then
    c%optical%X(1) = c%par%n_dust * c%par%mdust
    c%optical%X(2) = c%par%n_gas * c%par%X_HI
    c%optical%X(3) = c%par%n_gas * c%par%X_H2O
  else
    c%optical%X = 0D0
  end if
end subroutine update_local_coll



function Voigt(a, x)
  ! 2007MNRAS_375_1043Zaghloul
  double precision Voigt
  double precision, intent(in) :: x, a
  double precision v1, v2, u, du
  !
  if (a .le. 26.6D0) then
    v1 = exp(a*a - x*x) * erfc(a) * cos(2D0 * a * x)
  else
    v1 = (1D0 - 1D0 / (2D0 * a**2) &
      + 3D0 / (4D0 * a**4) - 15D0 / (8D0 * a**6) &
      + 105D0 / (16D0 * a**8) - 945D0 / (32D0 * a**10) &
      + 10395D0 / (64D0 * a**12)) &
      / (a * sqrt(phy_Pi)) &
      * exp(-x*x) * cos(2D0 * a * x)
  end if
  v2 = 0D0
  if (x .gt. 0D0) then
    u = 0D0
    du = min(1D0/(2D0*a)/100D0, 0.001D0, x*1D-3)
    do
      if (u .gt. x) then
        exit
      end if
      v2 = v2 + du/6D0 * ( &
        exp((u-x)*(u+x)) * sin(2D0*a*(x-u)) + &
        exp((u+du*0.5D0-x)*(u+du*0.5D0+x)) * sin(2D0*a*(x-u-du*0.5D0)) * 4D0 + &
        exp((u+du-x)*(u+du+x)) * sin(2D0*a*(x-u-du)) &
        )
      u = u + du
    end do
  end if
  Voigt = v1 + v2 * (2D0 / sqrt(phy_Pi))
end function Voigt



subroutine montecarlo_do(mc, cstart)
  type(type_montecarlo_config), intent(inout) :: mc
  type(type_photon_packet) ph, ph0
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer :: cthis, cnext
  integer(kind=LongInt) i
  logical found, escaped, destructed
  double precision eph_acc
  !
  call init_random_seed
  !
  eph_acc = 0D0
  !
  if (mc%savephoton) then
    call openFileSequentialWrite(mc%fU, &
      combine_dir_filename(mc%mc_dir_out, mc%fname_photons), 512)
  end if
  !
  ph0%lam = star_0%lam(1)
  ph0%iSpec = 1
  !
  i = 0
  do ! i=1, mc%nph*3
    i = i + 1
    ! The exact number of photons will not be exactly equal to mc%nph, because
    ! each photon may have a different energy due to refinement at certain
    ! wavelength (Lya most likely).
    !
    call emit_a_photon(mc, ph0) ! ph%lam is in angstrom
    ph = ph0
    !
    eph_acc = eph_acc + ph%en
    if (eph_acc .gt. star_0%lumi) then
      exit
    end if
    !
    call save_photon(ph, mc)
    !
    if (mod(i, 10) .eq. 0) then
      write (*, &
        '(A, 4X, A, I10, 2X, "(", F0.4, "%)", 4X, A, ES14.6, I6, " of ", I6)') &
        CHAR(27)//'[A', "Monte Carlo...  Photon ", i, &
        eph_acc*1E2/star_0%lumi, "lam = ", ph%lam, ph%iSpec, star_0%n
    end if
    !
    ! Get the index for accessing the optical data
    ph%iKap = get_idx_for_kappa(ph%lam, dust_0)
    if (ph%iKap .eq. 0) then ! No optical data for this lambda.
      call save_photon(ph, mc)
      cycle
    end if
    !
    ! Enter the disk domain
    call enter_the_domain(ph, cstart, cnext, found)
    if (.not. found) then
      call save_photon(ph, mc)
      cycle
    end if
    cthis => cnext
    !
    ! Increase the crossing count
    cthis%optical%cr_count = cthis%optical%cr_count + 1
    !
    call walk_scatter_absorb_reemit(ph, cthis, cstart, mc%nmax_cross, escaped, destructed)
    if (escaped) then ! Photon escaped from the disk domain
      call save_photon(ph, mc)
    else if (destructed) then
      ! Do nothing
    else
      write(*,'(A/)') 'Premature end of photon transport!'
    end if
  end do
  !
  if (mc%savephoton) then
    close(mc%fU)
  end if
  !
end subroutine montecarlo_do



subroutine enter_the_domain(ph, cstart, cnext, found)
  ! The photon location will be updated to the entry point if it does enter the
  ! cell.
  ! The entry point is shifted a little bit along the photon propagation
  ! direction.
  ! cnext will point to the leaf cell containing the entry point.
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer, intent(out) :: cnext
  logical, intent(out) :: found
  double precision length, r, z, eps
  integer dirtype
  !
  call calc_intersection_ray_cell(ph%ray, cstart, length, r, z, eps, found, dirtype)
  if (found) then
    call locate_photon_cell(r, z, cstart, cnext, found)
    if (found) then
      ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
      ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
      ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    end if
  end if
end subroutine enter_the_domain



function get_reemit_lam(T0, T1, kph, lut, dust, idx1)
  double precision get_reemit_lam
  double precision, intent(in) :: T0, T1, kph
  type(type_LUT_Tdust), intent(in) :: lut
  type(type_optical_property), intent(in) :: dust
  integer, intent(in) :: idx1
  integer i, ilam, idx0
  double precision a, b, c1, c0
  double precision r
  !
  call random_number(r)
  !
  idx0 = 0
  do i=idx1, 0, -1
    if (lut%Tds(i) .le. T0) then
      idx0 = i
      exit
    end if
  end do
  if (idx1 .lt. lut%n) then
    a = (T1 - lut%Tds(idx1)) / (lut%Tds(idx1+1) - lut%Tds(idx1))
  else
    a = 1D0
  end if
  c1 = (1D0-a) * lut%table(lut%m, idx1) + a * lut%table(lut%m, idx1+1)
  if (idx0 .eq. 0) then
    p4lam%pvals = &
      (kph + 1D0) * &
        ((1D0-a) * lut%table(:, idx1) + a * lut%table(:, idx1+1)) / c1
  else
    b = (T0 - lut%Tds(idx0)) / (lut%Tds(idx0+1) - lut%Tds(idx0))
    c0 = (1D0-b) * lut%table(lut%m, idx0) + b * lut%table(lut%m, idx0+1)
    p4lam%pvals = &
      (kph + 1D0) * &
        ((1D0-a) * lut%table(:, idx1) + a * lut%table(:, idx1+1)) / c1 - &
      kph * &
        ((1D0-b) * lut%table(:, idx0) + b * lut%table(:, idx0+1)) / c0
  end if
  p4lam%pvals = p4lam%pvals / p4lam%pvals(lut%m)
  ilam = get_a_sample(p4lam)
  if (ilam .eq. 0) then
    write(*,'(A)') 'ilam = 0:'
    write(*,'(21ES12.4,/)') p4lam%pvals(0:20)
  end if
  get_reemit_lam = dust%lam(ilam)
end function get_reemit_lam



subroutine walk_scatter_absorb_reemit(ph, c, cstart, imax, escaped, destructed)
  ! ph must be guaranteed to be inside c.
  ! An intersection between ph and c must exist, unless there is a numerical
  ! error.
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  type(type_cell), intent(in), pointer :: cstart
  integer(kind=LongInt), intent(in) :: imax
  logical, intent(out) :: escaped, destructed
  logical found, encountered
  type(type_cell), pointer :: cnext
  double precision tau_this, frac_abso
  integer(kind=LongInt) i
  integer itype
  double precision length, r, z, eps
  double precision rnd, tau, albedo, t
  integer dirtype
  !
  escaped = .false.
  destructed = .false.
  !
  call random_number(rnd)
  tau = -log(rnd)
  !
  ! imax is usually set to a large number
  do i=1, imax
    call calc_intersection_ray_cell(ph%ray, c, length, r, z, eps, found, dirtype)
    if (.not. found) then
      write(*,'(A, 6ES16.6/)') 'ph not in c: ', &
        sqrt(ph%ray%x**2+ph%ray%y**2), ph%ray%z, c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    if ((ph%iKap .gt. 0) .and. c%using) then
      tau_this = (c%optical%summed_sc(ph%iKap) + c%optical%summed_ab(ph%iKap)) * length * phy_AU2cm
    else
      tau_this = 0D0
    end if
    if (tau_this .ge. tau) then
      length = length * (tau/tau_this)
      encountered = .true.
      tau = 0D0
      !
      ph%e_count = ph%e_count + 1
      ph%ray%x = ph%ray%x + ph%ray%vx * length
      ph%ray%y = ph%ray%y + ph%ray%vy * length
      ph%ray%z = ph%ray%z + ph%ray%vz * length
    else
      encountered = .false.
      tau = tau - tau_this
      ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
      ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
      ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    end if
    !!!! Todo
    !if (ph%ray%z .lt. 0D0) then
    !  ph%ray%z  = -ph%ray%z
    !  ph%ray%vz = -ph%ray%vz
    !end if
    !
    if (c%using) then
      albedo = c%optical%summed_sc(ph%iKap) / &
        (c%optical%summed_sc(ph%iKap) + c%optical%summed_ab(ph%iKap) + 1D-100)
      frac_abso = tau2frac(tau_this) * (1D0 - albedo)
      !
      c%optical%en_prev = c%optical%en_gain_abso
      c%optical%en_gain_abso = c%optical%en_gain_abso + frac_abso * ph%en
      c%optical%flux(ph%iKap) = c%optical%flux(ph%iKap) + length * ph%en
      !
      ! Project the photon direction to local radial frame.
      t = sqrt(ph%ray%x**2 + ph%ray%y**2) + 1D-100
      c%optical%dir_wei(ph%iKap)%u = c%optical%dir_wei(ph%iKap)%u + &
        length * ph%en * (ph%ray%x * ph%ray%vx + ph%ray%y * ph%ray%vy) / t
      c%optical%dir_wei(ph%iKap)%v = c%optical%dir_wei(ph%iKap)%v + &
        length * ph%en * (ph%ray%x * ph%ray%vy - ph%ray%y * ph%ray%vx) / t
      c%optical%dir_wei(ph%iKap)%w = c%optical%dir_wei(ph%iKap)%w + length * ph%en * ph%ray%vz
    end if
    !write(*, '(I5, 6ES14.4/)') i, c%optical%en_gain_abso, c%par%mdust_cell, frac_abso, albedo, &
    !  c%optical%summed_sc(ph%iKap), c%optical%summed_ab(ph%iKap)
    !write(*,'(20X, 4ES15.5)') maxval(dust_0%ab), maxval(HI_0%ab), &
    !  maxval(water_0%ab), maxval(c%optical%summed_ab)
    !
    if (encountered) then
      call find_encounter_type(ph%iKap, c%optical, itype)
      select case (itype)
        case (1) ! Dust absorption
          call reemit_dust(ph, c, cstart, escaped)
          if (escaped) then
            return
          end if
        case (2) ! Dust scattering
          call get_reemit_dir_HenyeyGreenstein(ph%ray, dust_0%g(ph%iKap))
        case (3)
          write(*, '(A/)') 'H absorption: not possible!'
        case (4) ! H scattering
          call get_reemit_dir_uniform(ph%ray)
          ! write(*,'(A/)') 'Scattered by H atom!'
        case (5) ! Water absorption; no reemission
          destructed = .true.
          ! write(*,'(A/)') 'Absorbed by water!'
          return
        case (6)
          write(*, '(A/)') 'Water absorption: not possible!'
        case default
          write(*,'(A/)') 'Should not have this case.'
      end select
      call random_number(rnd)
      tau = -log(rnd)
    else
      !call locate_photon_cell(r, z, c, cnext, found)
      call locate_photon_cell_alt(r, z, c, dirtype, cnext, found)
      if (.not. found) then! Not entering a neighboring cell
        ! May be entering a non-neighboring cell
        call enter_the_domain(ph, cstart, cnext, found)
        if (.not. found) then ! Escape
          escaped = .true.
          return
        end if
      end if
      c => cnext
      c%optical%cr_count = c%optical%cr_count + 1
    end if
  end do
end subroutine walk_scatter_absorb_reemit



subroutine reemit_dust(ph, c, cstart, escaped)
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  type(type_cell), intent(in), pointer :: cstart
  logical, intent(out) :: escaped
  double precision Tdust_old
  integer idx
  double precision, parameter :: update_threshold = 1.01D0
  !if (c%optical%en_gain_abso .ge. update_threshold * c%optical%en_prev) then
    c%optical%ab_count = c%optical%ab_count + 1
    c%optical%kph = c%optical%en_prev / ph%en
    c%optical%en_gain = c%optical%en_gain + ph%en
    !
    Tdust_old = c%par%Tdust1
    !
    c%par%Tdust1 = get_Tdust_from_LUT( &
      c%optical%en_gain_abso / &
      (4D0*phy_Pi * c%par%mdust_cell), lut_0, idx)
    !
    ph%lam = get_reemit_lam(Tdust_old, c%par%Tdust1, c%optical%kph, lut_0, dust_0, idx)
    ph%iKap = get_idx_for_kappa(ph%lam, dust_0)
    !
    call get_reemit_dir_uniform(ph%ray)
  !else
  !  reemit = .false.
  !end if
end subroutine reemit_dust



subroutine emit_a_photon(mc, ph)
  type(type_montecarlo_config), intent(inout) :: mc
  type(type_photon_packet), intent(inout) :: ph
  !
  if ((ph%lam .lt. lam_range_UV(1)) .or. (ph%lam .gt. lam_range_UV(2))) then
    ! Not UV
    ph%en = mc%eph
    call get_next_lam(ph%lam, ph%iSpec, star_0, ph%en)
    ph%wei = 1D0
  else if ((ph%lam .lt. lam_range_LyA(1)) .or. (ph%lam .gt. lam_range_LyA(2))) then
    ! UV but not LyA
    ph%en = mc%eph * refine_UV
    call get_next_lam(ph%lam, ph%iSpec, star_0, ph%en)
    ph%wei = refine_UV
  else
    ! LyA
    ph%en = mc%eph * refine_LyA
    call get_next_lam(ph%lam, ph%iSpec, star_0, ph%en)
    ph%wei = refine_LyA
  end if
  ph%ray%x = 0D0
  ph%ray%y = 0D0
  ph%ray%z = 0D0
  ph%e_count = 0
  call get_emit_dir_uniform(ph%ray, mc%minw, mc%maxw)
end subroutine emit_a_photon



subroutine make_local_material_collection(loc, glo)
  type(type_local_encounter_collection), intent(inout) :: loc
  type(type_global_material_collection), intent(in) :: glo
  integer i
  loc%acc = 0D0
  loc%acc(:, 1) = glo%list(1)%ab * loc%X(1)
  loc%acc(:, 2) = glo%list(1)%sc * loc%X(1) + loc%acc(:, 1)
  loc%summed_ab = glo%list(1)%ab * loc%X(1)
  loc%summed_sc = glo%list(1)%sc * loc%X(1)
  do i=4, glo%ntype*2, 2
    loc%acc(:, i-1) = loc%acc(:, i-2) + glo%list(i/2)%ab * loc%X(i/2)
    loc%acc(:, i)   = loc%acc(:, i-1) + glo%list(i/2)%sc * loc%X(i/2)
    loc%summed_ab = loc%summed_ab + glo%list(i/2)%ab * loc%X(i/2)
    loc%summed_sc = loc%summed_sc + glo%list(i/2)%sc * loc%X(i/2)
  end do
  loc%summed = loc%acc(:, glo%ntype*2)
  do i=1, glo%ntype*2
    loc%acc(:, i) = loc%acc(:, i) / (loc%summed + 1D-100)
  end do
end subroutine make_local_material_collection



subroutine find_encounter_type(iKap, coll, itype)
  integer, intent(in) :: iKap
  type(type_local_encounter_collection), intent(in) :: coll
  integer, intent(out) :: itype
  double precision r
  integer i
  !
  call random_number(r)
  itype = coll%ntype
  do i=1, coll%ntype
    if (r .lt. coll%acc(iKap, i)) then
      itype = i
      return
    end if
  end do
end subroutine find_encounter_type



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
      case (-1)
        found = .false.
        return
      case default
        write(*, '(A/)') 'Should not have this case!'
    end select
    do i=1, neib%n
      if (is_inside_cell(r, z, cell_leaves%list(neib%idx(i))%p)) then
        cout => cell_leaves%list(neib%idx(i))%p
        found = .true.
        return
      end if
    end do
  end if
  call locate_photon_cell(r, z, c, cout, found)
end subroutine locate_photon_cell_alt



recursive subroutine locate_photon_cell(r, z, c, cout, found)
  ! Given r and z and start from c, find out a cell containing (r,z).
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  type(type_cell), pointer, intent(out) :: cout
  logical, intent(out) :: found
  integer i
  if (is_inside_cell(r, z, c)) then
    if (c%nChildren .eq. 0) then
      cout => c
      found = .true.
      return
    else
      do i=1, c%nChildren
        if (is_inside_cell(r, z, c%children(i)%p)) then
          call locate_photon_cell(r, z, c%children(i)%p, cout, found)
          if (found) then
            return
          end if
        end if
      end do
      cout => null()
      found = .false.
      return
    end if
  else
    if (associated(c%parent)) then
      call locate_photon_cell(r, z, c%parent, cout, found)
      return
    else
      cout => null()
      found = .false.
      return
    end if
  end if
end subroutine locate_photon_cell



subroutine calc_intersection_ray_cell(ray, c, length, r, z, eps, found, dirtype)
  type(type_ray), intent(in) :: ray
  type(type_cell), pointer, intent(in) :: c
  double precision, intent(out) :: length, r, z, eps
  logical, intent(out) :: found
  integer, intent(out) :: dirtype
  double precision rr, zz, A, B, C1, C2, D1, D2
  double precision, dimension(6) :: L
  double precision, parameter :: eps_ratio = 1D-10
  integer idx, i
  double precision, parameter :: FL = -1D0 ! False length
  !
  ! For intesection with top and bottom surfaces
  if (ray%vz .ne. 0D0) then
    L(1) = (c%ymax - ray%z) / ray%vz
    L(2) = (c%ymin - ray%z) / ray%vz
  else
    L(1) = FL
    L(2) = FL
  end if
  if (L(1) .ge. 0D0) then
    rr = (ray%x + L(1)*ray%vx)**2 + (ray%y + L(1)*ray%vy)**2
    if ((rr .lt. c%xmin*c%xmin) .or. (rr .gt. c%xmax*c%xmax)) then
      L(1) = FL
    end if
  end if
  if (L(2) .ge. 0D0) then
    rr = (ray%x + L(2)*ray%vx)**2 + (ray%y + L(2)*ray%vy)**2
    if ((rr .lt. c%xmin*c%xmin) .or. (rr .gt. c%xmax*c%xmax)) then
      L(2) = FL
    end if
  end if
  !
  ! For intersection with the inner and outer cylindrical boundaries
  A = ray%vx**2 + ray%vy**2
  B = 2D0 * (ray%x*ray%vx + ray%y*ray%vy)
  C1 = ray%x**2 + ray%y**2 - c%xmin**2
  C2 = ray%x**2 + ray%y**2 - c%xmax**2
  D1 = B*B - 4D0 * A * C1
  D2 = B*B - 4D0 * A * C2
  ! Inner cylinder
  if (D1 .gt. 0D0) then
    if (A .ne. 0D0) then
      L(3) = (-B + sqrt(D1)) / (2D0*A)
      zz = ray%z + ray%vz * L(3)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(3) = FL
      end if
      L(4) = (-B - sqrt(D1)) / (2D0*A)
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
    if (A .ne. 0D0) then
      L(5) = (-B + sqrt(D2)) / (2D0*A)
      zz = ray%z + ray%vz * L(5)
      if ((zz .lt. c%ymin) .or. (zz .gt. c%ymax)) then
        L(5) = FL
      end if
      L(6) = (-B - sqrt(D2)) / (2D0*A)
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
  rr = huge(0D0)
  idx = 0
  do i=1, 6
    if ((L(i) .ge. 0D0) .and. (L(i) .lt. rr)) then
      rr = L(i)
      idx = i
    end if
  end do
  if (idx .eq. 0) then
    found = .false.
    dirtype = -1
  else
    found = .true.
    length = L(idx)
    eps = eps_ratio * min(c%xmax-c%xmin, c%ymax-c%ymin)
    L(idx) = L(idx) + eps
    r = sqrt((ray%x + ray%vx * L(idx))**2 + (ray%y + ray%vy * L(idx))**2)
    z = ray%z + ray%vz * L(idx)
    if (.not. is_inside_cell(r, z, c)) then
      ! Exit the cell c instead of entering c
      ! 1: Exit through the top
      ! 2: Exit through the bottom
      ! 3: Exit through the inner edge
      ! 4: Exit through the inner edge also
      ! 5: Exit through the outer edge
      ! 6: Exit through the outer edge also
      dirtype = idx
    else
      dirtype = -1
    end if
  end if
end subroutine calc_intersection_ray_cell



recursive subroutine find_cell_by_pos(r, z, c, cout)
  ! Locate (r,z) in c; return the leaf cout that belongs to c.
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  type(type_cell), pointer, intent(out) :: cout
  integer i
  if (is_inside_cell(r, z, c)) then
    if (c%nChildren .eq. 0) then
      cout => c
      return
    else
      do i=1, c%nChildren
        if (is_inside_cell(r, z, c%children(i)%p)) then
          call find_cell_by_pos(r, z, c%children(i)%p, cout)
          exit
        end if
      end do
    end if
  else
    cout => null()
  end if
end subroutine find_cell_by_pos



function is_inside_cell(r, z, c)
  logical is_inside_cell
  double precision, intent(in) :: r, z
  type(type_cell), pointer, intent(in) :: c
  is_inside_cell = &
    (c%xmin .le. r) .and. (c%xmax .ge. r) .and. &
    (c%ymin .le. z) .and. (c%ymax .ge. z)
end function is_inside_cell



function get_idx_for_kappa(lam, dust)
  integer get_idx_for_kappa
  double precision, intent(in) :: lam
  type(type_optical_property), intent(in) :: dust
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  if ((lam .lt. dust%lam(1)) .or. (lam .gt. dust%lam(dust%n))) then
    get_idx_for_kappa = 0
    return
  else
    imin = 1
    imax = dust%n
    do i=1, dust%n
      if (imin .ge. imax-ITH) then
        do j=max(2,imin), imax
          if ((dust%lam(j-1) .le. lam) .and. (dust%lam(j) .gt. lam)) then
            get_idx_for_kappa = j-1
            return
          end if
        end do
        exit
      else
        imid = (imin + imax) / 2
        if (dust%lam(imid) .le. lam) then
          imin = imid
        else
          imax = imid
        end if
      end if
    end do
  end if
end function get_idx_for_kappa



subroutine get_next_lam(lamthis, idx, star, eph)
  ! star%lam(idx) <= lamthis < star%lam(idx+1)
  integer, intent(inout) :: idx
  double precision, intent(inout) :: lamthis
  double precision, intent(in) :: eph
  type(type_stellar_spectrum), intent(in) :: star
  integer i
  double precision val, tmp, v
  val = eph
  do i=idx, star%n-1
    v = (star%vals(i) + star%vals(i+1)) * 0.5D0
    tmp = v * (star%lam(i+1) - lamthis)
    if (tmp .ge. val) then
      lamthis = val/v + lamthis
      if (lamthis .ge. star%lam(i+1)) then
        idx = i + 1
      end if
      return
    else
      val = val - tmp
      lamthis = star%lam(i+1)
      idx = i + 1
    end if
  end do
end subroutine get_next_lam



function get_stellar_luminosity(star, lam1, lam2) result(lumi)
  double precision lumi
  type(type_stellar_spectrum), intent(in) :: star
  double precision, intent(in), optional :: lam1, lam2
  integer i
  lumi = 0D0
  do i=1, star%n-1
    if (present(lam1) .and. present(lam2)) then
      if ((star%lam(i)   .lt. lam1) .or. (star%lam(i)   .gt. lam2) .or. &
          (star%lam(i+1) .lt. lam1) .or. (star%lam(i+1) .gt. lam2)) then
        cycle
      end if
    end if
    lumi = lumi + &
      (star%vals(i+1) + star%vals(i)) * 0.5D0 * &
      (star%lam(i+1) - star%lam(i))
  end do
end function get_stellar_luminosity



subroutine load_stellar_spectrum(fname, star)
  character(len=*), intent(in) :: fname
  type(type_stellar_spectrum), intent(inout) :: star
  integer i, fU
  integer nrows, ios
  character(len=128) str
  !
  nrows = GetFileLen_comment_blank(fname, '!')
  star%n = nrows
  allocate(star%lam(star%n), star%vals(star%n))
  !
  call openFileSequentialRead(fU, fname, 99, 1)
  i = 0
  do
    call read_a_nonempty_row(fU, str, '(A128)', ios)
    if (ios .eq. 0) then
      i = i + 1
    else
      exit
    end if
    if (str(1:1) .eq. '!') then
      cycle
    end if
    read(str, '(F24.8, X, F24.8)') star%lam(i), star%vals(i)
  end do
  close(fU)
end subroutine load_stellar_spectrum




subroutine make_stellar_spectrum(lam0, lam1, nlam, star)
  ! lam must be in angstrom.
  double precision, intent(in) :: lam0, lam1
  integer, intent(in) :: nlam
  type(type_stellar_spectrum), intent(inout) :: star
  integer i
  double precision dlam, coeff
  !
  coeff = 4D0*phy_Pi**2 * (star%radius*phy_Rsun_CGS)**2
  star%n = nlam
  allocate(star%lam(star%n), star%vals(star%n))
  dlam = (lam1-lam0) / dble(nlam-1)
  do i=1, star%n
    if (i .eq. 1) then
      star%lam(i) = lam0
    else if (i .eq. star%n) then
      star%lam(i) = lam1
    else
      star%lam(i) = lam0 + dlam * dble(i-1)
    end if
    star%vals(i) = planck_B_lambda(star%T, star%lam(i)*phy_Angstrom2cm) * coeff * phy_Angstrom2cm
  end do
  !
end subroutine make_stellar_spectrum




function get_surf_max_angle()
  double precision get_surf_max_angle
  integer i, i0
  double precision r, z, w
  get_surf_max_angle = 0D0
  do i=1, surf_cells%nlen
    i0 = surf_cells%idx(i)
    r = cell_leaves%list(i0)%p%par%rmin
    z = cell_leaves%list(i0)%p%par%zmax
    w = z / sqrt(r*r + z*z)
    if (w .gt. get_surf_max_angle) then
      get_surf_max_angle = w
    end if
  end do
end function get_surf_max_angle




subroutine load_H2O_ab_crosssection(fname, wa)
  ! For the data format from Ted Bergin.
  character(len=*), intent(in) :: fname
  type(type_optical_property), intent(out) :: wa
  integer i, fU, ios
  double precision l1, l2, s
  !
  wa%n = GetFileLen(fname) - 2
  allocate(wa%lam(wa%n), wa%ab(wa%n), wa%sc(wa%n), wa%g(wa%n))
  wa%sc = 0D0
  wa%g = 0D0
  !
  call openFileSequentialRead(fU, fname, 99, 1)
  read(fU, *)
  read(fU, *)
  i = 0
  do
    read(fU, '(F6.1, 2X, F6.1, 2X, F6.1)', iostat=ios) l1, l2, s
    if (ios .ne. 0) then
      exit
    end if
    i = i + 1
    wa%lam(i) = (l1 + l2) * 0.5D0
    wa%ab(i) = s * 1D-18 ! Convert into cm2
    !write(*,*) i, wa%lam(i), wa%ab(i)
  end do
  close(fU)
end subroutine load_H2O_ab_crosssection



subroutine load_dust_data(fname, dust)
  character(len=*), intent(in) :: fname
  type(type_optical_property), intent(out) :: dust
  integer i, fU
  integer iformat, nrows, ios
  character(len=128) str
  call openFileSequentialRead(fU, fname, 299, 1)
  call read_a_nonempty_row(fU, str, '(A128)', ios)
  read(str, '(I16)') iformat
  call read_a_nonempty_row(fU, str, '(A128)', ios)
  read(str, '(I16)') nrows
  dust%n = nrows
  allocate(dust%lam(nrows), dust%ab(nrows), dust%sc(nrows), dust%g(nrows))
  dust%lam = 0D0
  dust%ab = 0D0
  dust%sc = 0D0
  dust%g = 0D0
  i = 0
  do
    call read_a_nonempty_row(fU, str, '(A128)', ios)
    if (ios .eq. 0) then
      i = i + 1
    else
      exit
    end if
    if ((str(1:1) .eq. '!') .or. (str(1:1) .eq. '#')) then
      cycle
    end if
    select case(iformat)
      case(1)
        read(str, '(2(F18.8, X))') dust%lam(i), dust%ab(i)
      case(2)
        read(str, '(3(F18.8, X))') dust%lam(i), dust%ab(i), dust%sc(i)
      case(3)
        read(str, '(4(F18.8, X))') dust%lam(i), dust%ab(i), dust%sc(i), dust%g(i)
    end select
  end do
  close(fU)
  dust%lam = dust%lam / phy_Angstrom2micron
end subroutine load_dust_data



subroutine make_H_Lya(T, hi)
  ! Zheng 2002
  ! A factor of c (speed of light) seems to be missing in that paper.
  type(type_optical_property), intent(out) :: hi
  double precision, intent(in) :: T
  double precision, parameter :: l0 = 1215.668D0
  double precision, parameter :: nu0 = 2.4660718D15
  double precision, parameter :: dnul = 9.938D7
  double precision, parameter :: f12 = 0.4162D0
  double precision dlam0, dlam, ratio
  double precision nu, x, dnu_th, a
  integer i, n2
  !
  hi%n = 51 ! Must be odd
  n2 = 25 ! Must be exactly (hi%n-1)/2
  !
  allocate(hi%lam(hi%n), hi%ab(hi%n), hi%sc(hi%n), hi%g(hi%n))
  hi%ab = 0D0
  hi%g = 0D0
  !
  hi%lam(hi%n/2+1) = l0
  dlam0 = (lam_range_LyA(2) - l0) * 1D-4
  ratio = get_ratio_of_interval_log(l0, lam_range_LyA(2), dlam0, n2)
  dlam = dlam0
  hi%lam(n2+1) = l0
  do i=n2+2, hi%n
    hi%lam(i) = hi%lam(i-1) + dlam
    hi%lam(hi%n - i + 1) = hi%lam(hi%n - i + 2) - dlam
    dlam = dlam * ratio
  end do
  !
  dnu_th = nu0 * sqrt(8D0*phy_kBoltzmann_SI*T/phy_Pi/phy_mProton_SI) / phy_SpeedOfLight_SI
  a = dnul / (2D0 * dnu_th)
  !
  do i=1, hi%n
    nu = phy_SpeedOfLight_SI / (hi%lam(i) * 1D-10)
    x = abs((nu - nu0)) / dnu_th
    hi%sc(i) = f12 * sqrt(phy_Pi) * phy_electronClassicalRadius_CGS * phy_SpeedOfLight_CGS &
      / dnu_th * Voigt(a, x)
    !write(*,*) i, hi%lam(i), x, a, hi%sc(i)
  end do
end subroutine make_H_Lya



function get_Tdust_Stefan_Boltzmann(flux)
  double precision get_Tdust_Stefan_Boltzmann
  double precision, intent(in) :: flux
  get_Tdust_Stefan_Boltzmann = sqrt(sqrt(flux / phy_StefanBoltzmann_CGS))
end function get_Tdust_Stefan_Boltzmann



function get_Tdust_from_LUT(val, lut, idx)
  double precision get_Tdust_from_LUT
  double precision, intent(in) :: val
  type(type_LUT_Tdust), intent(in) :: lut
  integer, intent(out) :: idx
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  !
  if (val .le. lut%vals(0)) then
    idx = 0
    get_Tdust_from_LUT = 0D0
    return
  else if (val .ge. lut%vals(lut%n)) then
    ! Using very crude extrapolation
    idx = lut%n
    get_Tdust_from_LUT = lut%Tds(lut%n) * val / lut%vals(lut%n)
    return
  else if (isnan(val)) then
    write(*,'(A//)') 'val is NaN!'
    return
  else
    imin = 1
    imax = lut%n
    do i=1, lut%n
      if (imin .ge. imax-ITH) then
        do j=imin, imax-1
          if ((lut%vals(j) .le. val) .and. (lut%vals(j+1) .ge. val)) then
            idx = j
            exit
          end if
        end do
        get_Tdust_from_LUT = lut%Tds(idx) + &
          (val - lut%vals(idx)) * &
          (lut%Tds(idx+1) - lut%Tds(idx)) / (lut%vals(idx+1) - lut%vals(idx))
        return
      else
        imid = (imin + imax) / 2
        if (lut%vals(imid) .le. val) then
          imin = imid
        else
          imax = imid
        end if
      end if
    end do
  end if
end function get_Tdust_from_LUT



subroutine make_LUT_Tdust(dust, lut)
  type(type_optical_property), intent(in) :: dust
  type(type_LUT_Tdust), intent(out) :: lut
  integer i, j
  double precision r, Tmin, Tmax, dT0, dT, tmp
  lut%n = 1028
  lut%m = dust%n - 1
  allocate(lut%Tds(0:lut%n), lut%vals(0:lut%n), lut%table(0:lut%m, 0:lut%n))
  lut%Tds(0) = 0D0
  lut%vals(0) = 0D0
  lut%table(0, :) = 0D0
  lut%table(:, 0) = 0D0
  Tmin = 1D0
  Tmax = 2D3
  dT0 = 1D-1
  dT = dT0
  r = get_ratio_of_interval_log(Tmin, Tmax, dT0, lut%n)
  do i=1, lut%n
    if (i .eq. 1) then
      lut%Tds(i) = Tmin
    else
      lut%Tds(i) = lut%Tds(i-1) + dT
    end if
    dT = dT * r
    lut%vals(i) = 0D0
    do j=1, lut%m
      tmp = &
        (dust%lam(j+1) - dust%lam(j)) * phy_Angstrom2cm * & ! dLambda
        (dust%ab(j+1) + dust%ab(j)) * 0.5D0 * & ! kappa
        planck_B_lambda(lut%Tds(i), & ! B_Lambda
                        (dust%lam(j+1)+dust%lam(j))*0.5D0*phy_Angstrom2cm)
      lut%table(j, i) = lut%table(j-1, i) + tmp
    end do
    lut%vals(i) = lut%table(lut%m, i)
  end do
end subroutine make_LUT_Tdust



subroutine convert_spec_to_distri(spec, distri, xtype)
  type(type_spectrum_generic), intent(in) :: spec
  type(type_distribution_table), intent(out) :: distri
  integer, intent(in) :: xtype
  integer i
  !
  distri%n = spec%n
  if (allocated(distri%pvals)) then
    deallocate(distri%pvals)
  end if
  allocate(distri%pvals(0:distri%n))
  if (xtype .eq. 1) then ! Using wavelength
    do i=1, distri%n
      distri%pvals(i) = spec%vals(i) * (spec%intervals(1, i) + spec%intervals(2, i))
    end do
  else if (xtype .eq. 2) then ! Using frequency
    do i=1, distri%n
      distri%pvals(i) = spec%vals(i) / (spec%intervals(1, i) + spec%intervals(2, i))
    end do
  else ! No unit conversion
    do i=1, distri%n
      distri%pvals(i) = spec%vals(i)
    end do
  end if
  distri%pvals(0) = 0D0
  do i=2, distri%n
    distri%pvals(i) = distri%pvals(i-1) + distri%pvals(i)
  end do
  distri%pvals = distri%pvals / distri%pvals(distri%n)
end subroutine convert_spec_to_distri




function get_a_sample(distri)
  integer get_a_sample
  type(type_distribution_table), intent(in) :: distri
  double precision r
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  ! The random number generator must be initialized somewhere else.
  call random_number(r)
  ! Binary search
  imin = 1
  imax = distri%n
  do i=1, distri%n
    if (imin .ge. imax-ITH) then
      do j=imin, imax
        if ((distri%pvals(j-1) .le. r) .and. (distri%pvals(j) .gt. r)) then
          get_a_sample = j
          return
        end if
      end do
      exit
    else
      imid = (imin + imax) / 2
      if (distri%pvals(imid) .le. r) then
        imin = imid
      else
        imax = imid
      end if
    end if
  end do
end function get_a_sample



function tau2frac(tau)
  double precision tau2frac
  double precision, intent(in) :: tau
  if (tau .le. 1D-2) then
    tau2frac = tau
  else
    tau2frac = 1D0 - exp(-tau)
  end if
end function tau2frac



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



function planck_B_nu(T, nu)
  double precision planck_B_nu
  double precision, intent(in) :: T, nu
  double precision tmp
  double precision, parameter :: TH = 1D-8
  tmp = (phy_hPlanck_CGS * nu) / (phy_kBoltzmann_CGS * T)
  if (tmp .gt. TH) then
    tmp = exp(tmp) - 1D0
  end if
  planck_B_nu = &
    2D0*phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / tmp
end function planck_B_nu




function planck_B_lambda(T, lambda_CGS)
  double precision planck_B_lambda
  double precision, intent(in) :: T, lambda_CGS
  double precision tmp
  double precision, parameter :: TH = 1D-8
  tmp = (phy_hPlanck_CGS * phy_SpeedOfLight_CGS) / (lambda_CGS * phy_kBoltzmann_CGS * T)
  if (tmp .gt. TH) then
    tmp = exp(tmp) - 1D0
  end if
  planck_B_lambda = &
    2D0*phy_hPlanck_CGS * phy_SpeedOfLight_CGS**2 / lambda_CGS**5 / tmp
end function planck_B_lambda



subroutine get_reemit_dir_uniform(ray)
  type(type_ray), intent(inout) :: ray
  double precision s, t
  double precision, dimension(2) :: x
  do
    call random_number(x)
    x(1) = x(1)*2D0 - 1D0
    x(2) = x(2)*2D0 - 1D0
    s = x(1) * x(1) + x(2) * x(2)
    if (s .le. 1D0) then
      t = sqrt(1D0 - s)
      ray%vx = 2D0 * x(1) * t
      ray%vy = 2D0 * x(2) * t
      ray%vz = 1D0 - 2D0 * s
      exit
    end if
  end do
end subroutine get_reemit_dir_uniform



subroutine get_emit_dir_uniform(ray, minw, maxw)
  type(type_ray), intent(inout) :: ray
  double precision, intent(in) :: minw, maxw
  double precision s, t
  double precision, dimension(2) :: x
  do
    call random_number(x)
    x(1) = x(1)*2D0 - 1D0
    x(2) = x(2)*2D0 - 1D0
    s = x(1) * x(1) + x(2) * x(2)
    if (s .le. 1D0) then
      ray%vz = 1D0 - 2D0 * s
      if ((ray%vz .le. maxw) .and. &
          (ray%vz .ge. minw)) then
        t = sqrt(1D0 - s)
        ray%vx = 2D0 * x(1) * t
        ray%vy = 2D0 * x(2) * t
        exit
      end if
    end if
  end do
end subroutine get_emit_dir_uniform



subroutine get_reemit_dir_HenyeyGreenstein(ray, g)
  use phy_const
  implicit none
  type(type_ray), intent(inout) :: ray
  double precision, intent(in) :: g
  double precision t, phi, costheta, sintheta
  double precision, dimension(2) :: p
  type(type_direction_cartesian) dir0, dir_rel
  double precision, parameter :: g_SMALL = 1D-2
  !
  dir0%u = ray%vx
  dir0%v = ray%vy
  dir0%w = ray%vz
  call random_number(p)
  if (abs(g) .gt. g_SMALL) then
    t = (1D0 - g*g) / (1D0 + g * (2D0 * p(1) - 1D0))
    costheta = 0.5D0 / g * &
               (1D0 + g * g - t * t)
  else
    costheta = p(1)*2D0 - 1D0
  end if
  sintheta = sqrt(1D0 - costheta * costheta)
  phi = (2D0 * phy_Pi) * p(2)
  dir_rel%u = sintheta * cos(phi)
  dir_rel%v = sintheta * sin(phi)
  dir_rel%w = costheta
  dir0 = get_anisotropic_dir(dir0, dir_rel)
  ray%vx = dir0%u
  ray%vy = dir0%v
  ray%vz = dir0%w
  ! Not sure about this.
  !if (ray%z .lt. 0D0) then
  !  ray%z  = -ray%z
  !  ray%vz = -ray%vz
  !end if
end subroutine get_reemit_dir_HenyeyGreenstein



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


  type(type_direction_cartesian) function get_anisotropic_dir(dir0, dir_rel)
    type(type_direction_cartesian) dir0, dir_rel
    type(type_sphere_coor_quat) d
    d = cartesian_to_spherical_quat(dir0)
    get_anisotropic_dir = rot_around_Y(dir_rel, d%costheta, d%sintheta)
    get_anisotropic_dir = rot_around_Z(get_anisotropic_dir, &
      d%cosphi, d%sinphi)
  end function get_anisotropic_dir


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




subroutine update_cell_and_photon(ph, c, cstart, itype, reemit, escaped)
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  type(type_cell), intent(in), pointer :: cstart
  integer, intent(in) :: itype
  logical, intent(out) :: reemit, escaped
  double precision Tdust_old
  double precision tau_TB, tau_IO, area_eff
  integer idx
  select case(itype)
    case(1) ! Absorption by dust
      ! Update dust temperature
      ! Re-emit a new photon
      c%optical%ab_count = c%optical%ab_count + 1
      c%optical%kph = c%optical%en_gain / ph%en
      !
      c%optical%en_gain = c%optical%en_gain + ph%en
      !if (c%optical%en_gain .ge. 1.1D0 * c%optical%en_prev) then
        c%optical%en_prev = c%optical%en_gain
        Tdust_old = c%par%Tdust1
        c%par%Tdust1 = get_Tdust_from_LUT(c%optical%en_gain/(4D0*phy_Pi*c%par%mdust_cell), lut_0, idx)
        if (ph%iKap .gt. 0) then
          tau_TB = c%par%dz * phy_AU2cm * c%optical%summed_ab(ph%iKap)
          tau_IO = c%par%dr * phy_AU2cm * c%optical%summed_ab(ph%iKap)
          if ((tau_TB .gt. 2D0) .and. (tau_IO .gt. 2D0)) then
            c%par%Tdust1 = get_Tdust_Stefan_Boltzmann(c%optical%en_gain/c%par%surf_area)
            !area_eff = 0D0
            !if (tau_TB .gt. 1D-3) then
            !  area_eff = area_eff + (c%par%area_T + c%par%area_B) * (1D0 - exp(-tau_TB))
            !else
            !  area_eff = area_eff + (c%par%area_T + c%par%area_B) * tau_TB
            !end if
            !if (tau_IO .gt. 1D-3) then
            !  area_eff = area_eff + (c%par%area_I + c%par%area_O) * (1D0 - exp(-tau_IO))
            !else
            !  area_eff = area_eff + (c%par%area_I + c%par%area_O) * tau_IO
            !end if
            !c%par%Tdust1 = get_Tdust_Stefan_Boltzmann(c%optical%en_gain / area_eff)
          end if
        end if
      !end if
      !if (max(c%xmax-c%xmin, c%ymax-c%ymin) * c%par%n_gas * phy_AU2cm .ge. 1D23) then
      !  c%par%Tdust1 = get_Tdust_Stefan_Boltzmann(c%optical%en_gain/c%par%surf_area)
      !end if
      !
      ph%lam = get_reemit_lam(Tdust_old, c%par%Tdust1, c%optical%kph, lut_0, dust_0, idx)
      ph%iKap = get_idx_for_kappa(ph%lam*phy_Angstrom2micron, dust_0)
      call get_reemit_dir_uniform(ph%ray)
      call send_photon_outof_cell(ph, c, cstart, escaped)
      reemit = .true.
    case(2) ! Scattering by dust
      ! Re-emit a new photon
      reemit = .true.
      call get_reemit_dir_HenyeyGreenstein(ph%ray, dust_0%g(ph%iKap))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !call send_photon_outof_cell(ph, c, cstart, escaped)
    !
    case(3) ! Absorption by H atom
      ! Do nothing
      reemit = .false.
    case(4) ! Scattering by H atom
      ! Re-emit a new photon
      call get_reemit_dir_uniform(ph%ray)
      reemit = .true.
    case(5) ! Absorption by H2O
      ! 
      reemit = .false.
    case(6) ! Scattering by H2O
      ! Do nothing
      call get_reemit_dir_uniform(ph%ray)
      reemit = .true.
  end select
end subroutine update_cell_and_photon



subroutine send_photon_outof_cell(ph, c, cstart, escaped)
  ! Before calling this subroutine ph should be inside c.
  ! This subroutine moves ph out of c according to the direction of ph, and
  ! repoint c to the new cell that ph resides in.
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer :: cnext
  double precision length, r, z, eps
  logical found, escaped
  integer dirtype
  !
  double precision tau_abso, frac_abso
  !
  escaped = .false.
  call calc_intersection_ray_cell(ph%ray, c, length, r, z, eps, found, dirtype)
  !
    if (c%using) then
      tau_abso = c%optical%summed_ab(ph%iKap) * length * phy_AU2cm
      if (tau_abso .le. 1D-4) then
        frac_abso = tau_abso
      else
        frac_abso = 1D0 - exp(-tau_abso)
      end if
      c%optical%en_gain_abso = c%optical%en_gain_abso + frac_abso * ph%en
    end if
  !
  if (found) then
    ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
    ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
    ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    ! Not sure about this.
    if (ph%ray%z .lt. 0D0) then
      ph%ray%z  = -ph%ray%z
      ph%ray%vz = -ph%ray%vz
    end if
    !
    call locate_photon_cell(r, z, c, cnext, found)
    if (.not. found) then! Not entering a neighboring cell
      call enter_the_domain(ph, cstart, cnext, found)
      if (.not. found) then ! Escape
        escaped = .true.
        return
      end if
    end if
    c => cnext
  else
    write(*,'(A/)') 'This should not happen!  In send_photon_outof_cell.'
  end if
end subroutine send_photon_outof_cell



subroutine copy_resample(n1, x1, y1, n2, x2, y2)
  ! copy the data in x1 and y1 into x2 and y2 with interpolation and resampling
  ! x1 and x2 must be sorted in ascending order
  integer, intent(in) :: n1, n2
  double precision, dimension(n1), intent(in)    :: x1, y1
  double precision, dimension(n2), intent(inout) :: x2, y2
  integer i, ibg, j
  ibg = 1
  do i=1, n2
    if (x2(i) .lt. x1(1)) then
      cycle
    end if
    if (x2(i) .gt. x1(n1)) then
      return
    end if
    do j=ibg, n1-1
      if ((x1(j) .le. x2(i)) .and. (x1(j+1) .ge. x2(i))) then
        y2(i) = y1(j) + (x2(i) - x1(j)) * (y1(j+1) - y1(j)) / (x1(j+1) - x1(j))
        ibg = j
        exit
      end if
    end do
  end do
end subroutine copy_resample



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


end module montecarlo
