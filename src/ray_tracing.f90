module ray_tracing

use data_struct
use grid
use chemistry
use statistic_equilibrium
use montecarlo
use lamda
use hitran
use cdms

implicit none

type :: type_cube
  integer nx, ny, nf
  double precision xmin, xmax, dx, ymin, ymax, dy
  double precision f0, fmin, fmax, df
  double precision view_theta
  integer itr
  type(type_rad_transition), pointer :: rapar => null()
  double precision, dimension(:,:,:), allocatable :: val
end type type_cube


type(type_mole_exc_conf) :: raytracing_conf

type(type_molecule_exc), private :: mole_exc

character(len=256) dir_name_log

namelist /raytracing_configure/ &
  raytracing_conf


contains


subroutine make_cubes_continuum
  integer i, j, k
  integer nx, ny, nf, nth
  !
  double precision :: xmin, xmax, ymin, ymax
  double precision tmp_coeff, lam_max, lam_min
  character(len=128) im_dir, fname
  type(type_cube) :: cube
  double precision, dimension(:), allocatable :: vec_flux
  !
  nf  = raytracing_conf%nlam
  nth = raytracing_conf%nth
  nx  = raytracing_conf%nx
  ny  = raytracing_conf%ny
  !
  if (abs(raytracing_conf%maxx) .le. 1D-8) then
    xmax = max(root%xmax, root%ymax)
  else
    xmax = raytracing_conf%maxx
  end if
  if (abs(raytracing_conf%maxy) .le. 1D-8) then
    ymax = max(root%xmax, root%ymax)
  else
    ymax = raytracing_conf%maxy
  end if
  xmin = -xmax
  ymin = -ymax
  !
  im_dir = trim(combine_dir_filename(raytracing_conf%dir_save_image, 'images/'))
  if (.not. dir_exist(im_dir)) then
    call my_mkdir(im_dir)
  end if
  !
  cube%nx = nx
  cube%ny = ny
  cube%nf = nf
  cube%xmin = xmin
  cube%xmax = xmax
  cube%ymin = ymin
  cube%ymax = ymax
  cube%dx = (xmax - xmin) / dble(nx-1)
  cube%dy = (ymax - ymin) / dble(ny-1)
  !
  allocate(cube%val(nx, ny, nf), vec_flux(nf))
  !
  tmp_coeff = cube%dx*cube%dy * (phy_AU2cm / (raytracing_conf%dist*phy_pc2cm))**2 / &
              phy_jansky2CGS
  !
  do i=1, raytracing_conf%nlam_window
    do j=1, nth ! Viewing angles
      cube%view_theta = raytracing_conf%view_thetas(j)
      !
      lam_min = raytracing_conf%lam_mins(i)
      lam_max = raytracing_conf%lam_maxs(i)
      !
      write(*, '(2(A, I6, A, I6, /), A, 2F12.2/, A, F7.2/)') &
        'lam', i, ' of ', raytracing_conf%nlam_window, &
        'angl', j, ' of ', nth, &
        'lam range = ', lam_min, lam_max, &
        'theta = ', cube%view_theta
      !
      cube%fmin = phy_SpeedOfLight_SI/(lam_max*1D-6)
      cube%fmax = phy_SpeedOfLight_SI/(lam_min*1D-6)
      cube%df = (cube%fmax - cube%fmin) / dble(nf-1)
      !
      call make_a_cube(cube, is_line=.false.)
      !
      do k=1, nf ! Frequency channels
        vec_flux(k) = sum(cube%val(:,:,k)) * tmp_coeff
      end do
      !
      write(fname, '(A, I0.5, "_", 2(F09.2, "_"), I0.5, "_", F09.2, ".fits")') &
        'cont_', i, lam_min, lam_max, j, cube%view_theta
      call dropout_char(fname, ' ')
      !
      if (raytracing_conf%save_spectrum_only) then
        call save_cube_to_fits_spec_only(trim(combine_dir_filename(im_dir, fname)), &
          cube, vec_flux, is_line=.false.)
      else
        call save_cube_to_fits(trim(combine_dir_filename(im_dir, fname)), &
          cube, vec_flux, is_line=.false.)
      end if
      !
    end do
  end do
end subroutine make_cubes_continuum



subroutine make_cubes_line
  integer itr
  integer i, j, k
  integer nx, ny, nf, nth
  !
  double precision :: xmin, xmax, ymin, ymax
  double precision VeloHalfWidth_this, freq_w
  double precision tmp_coeff
  character(len=128) im_dir, fname
  type(type_cube) :: cube
  double precision, dimension(:,:), allocatable :: &
    arr_tau, Ncol_up, Ncol_low
  double precision, dimension(:), allocatable :: vec_flux
  double precision p_flux
  !
  nf  = raytracing_conf%nf
  nth = raytracing_conf%nth
  nx  = raytracing_conf%nx
  ny  = raytracing_conf%ny
  !
  if (abs(raytracing_conf%maxx) .le. 1D-8) then
    xmax = max(root%xmax, root%ymax)
  else
    xmax = raytracing_conf%maxx
  end if
  if (abs(raytracing_conf%maxy) .le. 1D-8) then
    ymax = max(root%xmax, root%ymax)
  else
    ymax = raytracing_conf%maxy
  end if
  xmin = -xmax
  ymin = -ymax
  !
  im_dir = trim(combine_dir_filename(raytracing_conf%dir_save_image, 'images/'))
  if (.not. dir_exist(im_dir)) then
    call my_mkdir(im_dir)
  end if
  !
  cube%nx = nx
  cube%ny = ny
  cube%nf = nf
  cube%xmin = xmin
  cube%xmax = xmax
  cube%ymin = ymin
  cube%ymax = ymax
  cube%dx = (xmax - xmin) / dble(nx-1)
  cube%dy = (ymax - ymin) / dble(ny-1)
  !
  allocate(cube%val(nx, ny, nf), &
           arr_tau(nx, ny), &
           Ncol_up(nx, ny), &
           Ncol_low(nx, ny), &
           vec_flux(nf))
  !
  tmp_coeff = cube%dx*cube%dy * (phy_AU2cm / (raytracing_conf%dist*phy_pc2cm))**2 / &
              phy_jansky2CGS
  !
  do i=1, mole_exc%ntran_keep ! Transitions
    itr = mole_exc%itr_keep(i)
    cube%itr = itr
    cube%f0 = mole_exc%p%rad_data%list(itr)%freq
    cube%rapar => mole_exc%p%rad_data%list(itr)
    !
    do j=1, nth ! Viewing angles
      cube%view_theta = raytracing_conf%view_thetas(j)
      !
      write(*, '(2(A, I6, A, I6, /), A, ES16.8/, A, F7.2)') &
        'Tran', i, ' of ', mole_exc%ntran_keep, &
        'angl', j, ' of ', nth, &
        'freq (Hz) = ', cube%f0, &
        'theta (deg) = ', cube%view_theta
      !
      ! Kepler broadening + thermal/turbulent broadening
      VeloHalfWidth_this = raytracing_conf%VeloKepler * &
        sin(cube%view_theta * (phy_Pi / 180D0)) + &
        raytracing_conf%VeloTurb
      write(*,'(A, F9.2/)') 'HWHM = (km/s)', VeloHalfWidth_this/1D3
      freq_w = cube%f0 * VeloHalfWidth_this / phy_SpeedOfLight_SI
      cube%fmin = cube%f0 - freq_w
      cube%fmax = cube%f0 + freq_w
      cube%df = freq_w * 2D0 / dble(nf-1)
      !
      call make_a_cube(cube, arr_tau, Ncol_up, Ncol_low, is_line=.true.)
      !
      do k=1, nf ! Frequency channels
        vec_flux(k) = sum(cube%val(:,:,k)) * tmp_coeff
      end do
      !
      p_flux = maxval(vec_flux) - 0.5D0*(vec_flux(1)+vec_flux(nf))
      if (p_flux .ge. raytracing_conf%min_flux) then
        ! Only save strong lines, otherwise take too much storaage.
        write(fname, '(A, 3(I0.5,"_"), ES14.5, "_", F09.2, ".fits")') &
          'line_', i, itr, j, cube%f0, cube%view_theta
        call dropout_char(fname, ' ')
        !
        if (raytracing_conf%save_spectrum_only) then
          call save_cube_to_fits_spec_only(trim(combine_dir_filename(im_dir, fname)), &
            cube, vec_flux, arr_tau, Ncol_up, Ncol_low, is_line=.true.)
        else
          call save_cube_to_fits(trim(combine_dir_filename(im_dir, fname)), &
            cube, vec_flux, arr_tau, Ncol_up, Ncol_low, is_line=.true.)
        end if
      else
        write(*, '(A)') 'Line too weak!  Will not save it.'
        write(*, '(A, ES16.6/)') 'Peak flux (baseline removed):', p_flux
      end if
      !
    end do
  end do
end subroutine make_cubes_line



subroutine make_a_cube(cb, arr_tau, Ncol_up, Ncol_low, is_line)
  type(type_cube), intent(inout) :: cb
  double precision, dimension(:,:), intent(out), optional :: arr_tau, Ncol_up, Ncol_low
  logical, intent(in), optional :: is_line
  integer i, j, k
  double precision x, y, z
  double precision costheta, sintheta
  double precision tau
  double precision, dimension(:), allocatable :: Inu_0
  double precision Nup, Nlow
  type(type_photon_ray_multi) ph
  logical is_l
  !
  if (present(is_line)) then
    is_l = is_line
  else
    is_l = .true.
  end if
  !
  z = -max(root%xmax, root%ymax, &
           abs(cb%xmax), abs(cb%xmin), &
           abs(cb%ymax), abs(cb%ymin)) * 5D0
  !
  costheta = cos(cb%view_theta * (phy_Pi / 180D0))
  sintheta = sin(cb%view_theta * (phy_Pi / 180D0))
  !
  ph%ray%vx = 0D0
  ph%ray%vy = -sintheta
  ph%ray%vz =  costheta
  !
  ph%nf = cb%nf
  !
  if (is_l) then
    ph%iTran = cb%itr
  end if
  !
  allocate(ph%f(cb%nf), ph%lam(cb%nf), ph%iKap(cb%nf), ph%Inu(cb%nf), Inu_0(cb%nf))
  !
  do i=1, cb%nf
    ph%f(i) = cb%fmin + cb%df * dble(i-1)
    ph%lam(i) = phy_SpeedOfLight_CGS / phy_Angstrom2cm / ph%f(i) ! in Angstrom
    ph%iKap(i) = get_idx_for_kappa(ph%lam(i), dust_0)
    Inu_0(i) = planck_B_nu(phy_CMB_T, ph%f(i))
  end do
  !
  y = cb%ymin
  do j=1, cb%ny
    !
    write(*, '(A, 10X, A, I6, " of ", I6)') &
      CHAR(27)//'[A', 'Row ', j, cb%ny
    !
    x = cb%xmin
    do i=1, cb%nx
      !
      ph%ray%x =  x
      ph%ray%y =  y * costheta - z * sintheta
      ph%ray%z =  y * sintheta + z * costheta
      !
      ph%Inu = Inu_0
      !
      call integerate_a_ray(ph, tau, Nup, Nlow, is_l)
      !
      cb%val(i, j, :)   = ph%Inu
      !
      if (is_l .and. present(Ncol_up)) then
        Ncol_up(i, j)  = Nup
        Ncol_low(i, j) = Nlow
        arr_tau(i, j)  = tau
      end if
      !
      x = x + cb%dx
    end do
    !
    y = y + cb%dy
    !
  end do
end subroutine make_a_cube



subroutine integerate_a_ray(ph, tau, Nup, Nlow, is_line)
  ! ph must be guaranteed to be inside c.
  ! An intersection between ph and c must exist, unless there is a
  ! numerical error.
  type(type_photon_ray_multi), intent(inout) :: ph
  double precision, intent(out) :: tau
  double precision, intent(out) :: Nup, Nlow
  logical, intent(in), optional :: is_line
  type(type_cell), pointer :: c
  type(type_cell), pointer :: cnext
  logical found
  double precision length, r, z, eps
  integer dirtype
  integer itr, ilow, iup, iL, iU, iKp
  double precision ylow, yup, Aul, Bul, Blu
  double precision f0, width_nu, line_alpha, line_J
  double precision tau_this, t1
  double precision, dimension(:), allocatable :: tau_s
  logical is_l
  !
  double precision cont_alpha, cont_J
  !
  integer i, ii
  !
  tau = 0D0
  Nup = 0D0
  Nlow = 0D0
  !
  call enter_the_domain_mirror(ph%ray, root, c, found)
  if (.not. found) then
    return
  end if
  !
  if (present(is_line)) then
    is_l = is_line
  else
    is_l = .true.
  end if
  !
  allocate(tau_s(ph%nf))
  tau_s = 0D0
  !
  if (is_l) then
    itr  = ph%iTran
    ilow = mole_exc%p%rad_data%list(itr)%ilow
    iup  = mole_exc%p%rad_data%list(itr)%iup
    iL   = mole_exc%ilv_reverse(ilow)
    iU   = mole_exc%ilv_reverse(iup)
    f0   = mole_exc%p%rad_data%list(itr)%freq ! Rest frequency
    Aul  = mole_exc%p%rad_data%list(itr)%Aul
    Bul  = mole_exc%p%rad_data%list(itr)%Bul
    Blu  = mole_exc%p%rad_data%list(itr)%Blu
  end if
  !
  do ii=1, root%nOffspring*2
    ! Get the intersection between the photon ray and the boundary of the cell
    ! that this photon resides in
    call calc_intersection_ray_cell_mirror(ph%ray, c, &
      length, r, z, eps, found, dirtype)
    if (.not. found) then
      write(*,'(A, I6, 10ES10.2/)') 'In integerate_a_ray, ph not cross c: ', &
        ii, &
        ph%ray%x, ph%ray%y, ph%ray%z, &
        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
        c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    !
    if (c%using) then
      !
      if (is_l) then
        call set_using_mole_params(mole_exc%p, c)
        !
        ylow = c%focc%vals(iL)
        yup  = c%focc%vals(iU)
        if (raytracing_conf%adjust_yup_ylow_nonLTE) then
          call do_adjust_yup_ylow_nonLTE( &
            yup, ylow, &
            mole_exc%p%level_list(iup)%weight, &
            mole_exc%p%level_list(ilow)%weight, &
            f0, c%par%Tgas, c%par%n_gas, raytracing_conf%n_critical_CGS, &
            c%cont_lut%J(ph%iKap(1)))
        end if
        !if (mole_exc%conf%useLTE) then
        !  if (abs((yup/ylow) / &
        !      (Blu/Bul * exp(-(mole_exc%p%rad_data%list(itr)%Eup- &
        !          mole_exc%p%rad_data%list(itr)%Elow)/c%par%Tgas))-1D0) &
        !      .ge. 1D-2) then
        !    write(*, '(A, 11ES16.6E3/)') &
        !      'Wrong: ', ylow, yup, Blu, Bul, &
        !      mole_exc%p%rad_data%list(itr)%Elow, &
        !      mole_exc%p%rad_data%list(itr)%Eup, &
        !      c%par%Tgas, &
        !      (yup/ylow) / &
        !      (Blu/Bul * exp(-(mole_exc%p%rad_data%list(itr)%Eup- &
        !          mole_exc%p%rad_data%list(itr)%Elow)/c%par%Tgas)), &
        !      (yup/ylow) / &
        !      (Blu/Bul * exp(-(mole_exc%p%level_list(iup)%energy- &
        !          mole_exc%p%level_list(ilow)%energy)/c%par%Tgas)), &
        !      Blu/Bul, &
        !      mole_exc%p%level_list(iup)%weight/mole_exc%p%level_list(ilow)%weight
        !  end if
        !end if
        !
        Nup  = Nup  + (mole_exc%p%density_mol * length * phy_AU2cm) * yup
        Nlow = Nlow + (mole_exc%p%density_mol * length * phy_AU2cm) * ylow
      end if
      !
      do i=1, ph%nf
        !
        ! Get the continuum kappa and j.
        iKp = ph%iKap(i)
        if ((iKp .gt. 0) .and. c%using) then
          if (iKp .lt. dust_0%n) then
            ! Simple linear interpolation.
            t1 = (ph%lam(i) - dust_0%lam(iKp)) / &
                 (dust_0%lam(iKp+1) - dust_0%lam(iKp))
            cont_alpha = &
              (c%optical%ext_tot(iKp+1) - c%optical%ext_tot(iKp)) * t1 + &
              c%optical%ext_tot(iKp)
            cont_J = cont_alpha * ( &
              (c%cont_lut%J(iKp+1) - c%cont_lut%J(iKp)) * t1 + &
              c%cont_lut%J(iKp))
          else
            cont_alpha = c%optical%ext_tot(iKp)
            cont_J = c%cont_lut%J(iKp) * cont_alpha
          end if
        else
          cont_alpha = 0D0
          cont_J = 0D0
        end if
        !
        if (.not. is_l) then
          call integrate_one_step(ph%Inu(i), tau_this, &
            cont_J, cont_alpha, length*phy_AU2cm)
        else
          ! Get the line kappa and j.
          ! Local line broadening
          width_nu = f0 * mole_exc%p%dv / phy_SpeedOfLight_CGS
          !
          ! Rybicki & Lightman, p31
          t1 = phy_hPlanck_CGS * f0 / (4D0*phy_Pi) * &
               mole_exc%p%density_mol / (phy_sqrt2Pi * width_nu)
          line_alpha = t1 * (ylow * Blu - yup  * Bul)
          line_J     = t1 * yup * Aul
          !if (line_alpha .lt. 0D0) then
          !  write(*, '(A, 10ES20.8E3/)') &
          !    'Maser: ', ylow, yup, Blu, Bul, &
          !    mole_exc%p%rad_data%list(itr)%Elow, &
          !    mole_exc%p%rad_data%list(itr)%Eup, &
          !    c%par%Tgas, &
          !    (yup/ylow) / &
          !    (Blu/Bul * exp(-(mole_exc%p%rad_data%list(itr)%Eup- &
          !        mole_exc%p%rad_data%list(itr)%Elow)/c%par%Tgas)), &
          !    line_alpha, line_J
          !end if
          !
          call integrate_line_within_one_cell(ph, i, length, f0, width_nu, &
            line_alpha, line_J, cont_alpha, cont_J, tau_this)
        end if
        !
        tau_s(i) = tau_s(i) + tau_this
      end do
    end if
    !
    ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
    ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
    ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    !
    call locate_photon_cell_mirror(r, z, c, cnext, found)
    if (.not. found) then! Not entering a neighboring cell
      ! May be entering a non-neighboring cell?
      call enter_the_domain_mirror(ph%ray, root, cnext, found)
      if (.not. found) then ! Escape
        !
        if (raytracing_conf%subtract_cont_tau) then
          tau = maxval(tau_s) - 0.5D0*(tau_s(1) + tau_s(ph%nf))
        else
          tau = maxval(tau_s)
        end if
        !
        return
      end if
    end if
    !
    c => cnext
    !
  end do
  !
  write(*, '(A)') 'Error in integerate_a_ray:'
  write(*, '(A)') 'Should never reach here!'
  write(*,'(I6, 10ES10.2/)') &
    i, &
    ph%ray%x, ph%ray%y, ph%ray%z, &
    ph%ray%vx, ph%ray%vy, ph%ray%vz, &
    c%xmin, c%xmax, c%ymin, c%ymax
  write(*,'(4ES10.2, I4, L4/)') sqrt(r), z, length, eps, dirtype, found
  stop
end subroutine integerate_a_ray




pure subroutine integrate_line_within_one_cell(ph, ifr, length, f0, width_nu, &
        line_alpha, line_J, cont_alpha, cont_J, tau)
  type(type_photon_ray_multi), intent(inout) :: ph
  integer, intent(in) :: ifr
  double precision, intent(in) :: length, f0, width_nu, line_alpha, line_J, cont_alpha, cont_J
  double precision, intent(out) :: tau
  double precision nu, nu1, nu2, nu3, numin, numax, dnu, dl, dl_CGS, ltmp, x
  double precision jnu, knu, dtau
  type(type_ray) ray
  integer i, ndiv
  !
  nu1 = get_doppler_nu(a_star%mass, ph%f(ifr), ph%ray)
  !
  ray%vx = ph%ray%vx
  ray%vy = ph%ray%vy
  ray%vz = ph%ray%vz
  !
  ray%x = ph%ray%x + ph%ray%vx * length
  ray%y = ph%ray%y + ph%ray%vy * length
  ray%z = ph%ray%z + ph%ray%vz * length
  nu2 = get_doppler_nu(a_star%mass, ph%f(ifr), ray)
  !
  ray%x = ph%ray%x + ph%ray%vx * length * 0.5D0
  ray%y = ph%ray%y + ph%ray%vy * length * 0.5D0
  ray%z = ph%ray%z + ph%ray%vz * length * 0.5D0
  nu3 = get_doppler_nu(a_star%mass, ph%f(ifr), ray)
  !
  numin = min(nu1,nu2,nu3)
  numax = max(nu1,nu2,nu3)
  !
  x = min(abs(nu1-f0), abs(nu2-f0), abs(nu3-f0)) / width_nu
  if (x .gt. 10D0) then
    call integrate_one_step(ph%Inu(ifr), tau, cont_J, cont_alpha, length*phy_AU2cm)
  else
    ndiv = 2 + int(10D0 * (numax-numin) / width_nu)
    dl = length / dble(ndiv)
    dl_CGS = dl * phy_AU2cm
    ltmp = 0.5D0 * dl
    tau = 0D0
    do i=1, ndiv
      ray%x = ph%ray%x + ph%ray%vx * ltmp
      ray%y = ph%ray%y + ph%ray%vy * ltmp
      ray%z = ph%ray%z + ph%ray%vz * ltmp
      !
      nu = get_doppler_nu(a_star%mass, ph%f(ifr), ray)
      !
      call integrate_one_step_line(ph%Inu(ifr), dtau, (nu - f0) / width_nu, &
                 dl_CGS, line_J, line_alpha, cont_J, cont_alpha)
      tau = tau + dtau
      ltmp = ltmp + dl
    end do
  end if
end subroutine integrate_line_within_one_cell



pure subroutine integrate_one_step_line(I0, tau, x, dl, line_J, line_alpha, cont_J, cont_alpha)
  ! The medium within one step will be assumed to be uniform in every aspect.
  double precision, intent(inout) :: I0, tau
  double precision, intent(in) :: x, dl, line_J, line_alpha, cont_J, cont_alpha
  double precision t1, jnu, knu
  if ((x .gt. 10D0) .or. (x .lt. -10D0)) then
    call integrate_one_step(I0, tau, cont_J, cont_alpha, dl)
  else
    t1 = exp(-x*x*0.5D0)
    jnu = t1 * line_J + cont_J
    knu = t1 * line_alpha + cont_alpha
    call integrate_one_step(I0, tau, jnu, knu, dl)
  end if
end subroutine integrate_one_step_line



pure subroutine integrate_one_step(Inu, tau, jnu, knu, len_CGS)
  double precision, intent(inout) :: Inu, tau
  double precision, intent(in) :: jnu, knu, len_CGS
  double precision t1
  !
  tau = knu * len_CGS
  if (tau .ge. 1D-4) then
    if (tau .ge. 50D0) then
      Inu = jnu/knu
    else
      t1 = exp(-tau)
      Inu = Inu * t1 + jnu/knu * (1D0 - t1)
    end if
  else if (tau .lt. 0D0) then
    t1 = exp(-tau)
    Inu = Inu * t1 + jnu/knu * (1D0 - t1)
  else
    Inu = Inu * (1D0 - tau) + jnu * len_CGS
  end if
end subroutine integrate_one_step




subroutine save_cube_to_fits(filename, cube, vec_flux, arr_tau, Ncol_up, Ncol_low, is_line)
  use my_timer
  character(len=*), intent(in) :: filename
  type(type_cube), intent(in) :: cube
  double precision, dimension(:), intent(in) :: vec_flux
  double precision, dimension(:,:), intent(in), optional :: arr_tau, Ncol_up, Ncol_low
  logical, intent(in), optional :: is_line
  type(type_fits_par) :: fp
  logical is_l
  !
  double precision dv, vmax
  type(date_time) a_date_time
  !
  if (present(is_line)) then
    is_l = is_line
  else
    is_l = .true.
  end if
  !
  fp%filename = trim(filename)
  !
  fp%stat = 0
  fp%blocksize = 1
  fp%pcount = 0
  fp%gcount = 1
  fp%group=1
  fp%fpixel=1
  fp%decimals = 16
  fp%author = 'Fujun Du (fdu@umich.edu)'
  fp%user = ''
  fp%simple=.true.
  fp%extend=.true.
  fp%bitpix=-64 ! double
  !
  call ftgiou(fp%fU, fp%stat)
  call ftinit(fp%fU, fp%filename, fp%blocksize, fp%stat)
  !
  fp%naxis=3
  fp%naxes(1) = cube%nx
  fp%naxes(2) = cube%ny
  fp%naxes(3) = cube%nf
  fp%nelements = cube%nx * cube%ny * cube%nf
  !
  call ftphpr(fp%fU, fp%simple, fp%bitpix, fp%naxis, fp%naxes, &
              fp%pcount, fp%gcount, fp%extend, fp%stat)
  !
  ! The cube.
  call ftpprd(fp%fU, fp%group, fp%fpixel, fp%nelements, cube%val, fp%stat)
  !
  call ftpkyd(fp%fU, 'CDELT1', cube%dx,  fp%decimals, 'dx (AU)', fp%stat)
  call ftpkyd(fp%fU, 'CDELT2', cube%dy,  fp%decimals, 'dy (AU)', fp%stat)
  call ftpkyd(fp%fU, 'CDELT3', cube%df,  fp%decimals, 'df (Hz)', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX3', 1.0D0,    fp%decimals, 'k0', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL1', cube%xmin,  fp%decimals, 'xmin', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL2', cube%ymin,  fp%decimals, 'ymin', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL3', cube%fmin,  fp%decimals, 'fmin', fp%stat)
  call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
  call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
  call ftpkys(fp%fU, 'CTYPE3', 'F', 'Hz', fp%stat)
  !
  call ftpkyd(fp%fU, 'Dist',  raytracing_conf%dist, fp%decimals, 'pc', fp%stat)
  call ftpkyd(fp%fU, 'Theta', cube%view_theta, fp%decimals, 'deg', fp%stat)
  call ftpkyd(fp%fU, 'MaxFlux', maxval(vec_flux),  fp%decimals, 'jy', fp%stat)
  if (is_l .and. present(arr_tau)) then
    call ftpkyd(fp%fU, 'MaxTau',  maxval(arr_tau),   fp%decimals, '', fp%stat)
  end if
  if (is_l) then
    call ftpkys(fp%fU, 'ExtName', 'LineCube', '', fp%stat)
    !call ftpkyj(fp%fU, 'Itr',   cube%itr,        'trans num', fp%stat)
    call ftpkyd(fp%fU, 'F0',    cube%f0,     fp%decimals, 'Hz', fp%stat)
    call ftpkyd(fp%fU, 'lam0',  cube%rapar%lambda, fp%decimals, 'Angstrom', fp%stat)
    call ftpkyd(fp%fU, 'Eup',   cube%rapar%Eup,  fp%decimals, 'K', fp%stat)
    call ftpkyd(fp%fU, 'Elow',  cube%rapar%Elow,  fp%decimals, 'K', fp%stat)
    !call ftpkyj(fp%fU, 'iup',   cube%rapar%iup,  '', fp%stat)
    !call ftpkyj(fp%fU, 'ilow',  cube%rapar%ilow, '', fp%stat)
    call ftpkyd(fp%fU, 'Aul',   cube%rapar%Aul,  fp%decimals, 's-1', fp%stat)
    call ftpkyd(fp%fU, 'Bul',   cube%rapar%Bul,  fp%decimals, '', fp%stat)
    call ftpkyd(fp%fU, 'Blu',   cube%rapar%Blu,  fp%decimals, '', fp%stat)
    call ftpkys(fp%fU, 'Qnum',  trim(cube%rapar%qnum), '', fp%stat)
  else
    call ftpkys(fp%fU, 'ExtName', 'ContCube', '', fp%stat)
  end if
  !
  !call ftpkys(fp%fU, 'Author', fp%author, '', fp%stat)
  call ftpkys(fp%fU, 'User',   fp%user,   '', fp%stat)
  call ftpkys(fp%fU, 'SavedAt', trim(a_date_time%date_time_str()), '', fp%stat)
  !
  ! First extension: tau map
  if (is_l .and. present(arr_tau)) then
    call ftcrhd(fp%fU, fp%stat)
    fp%naxis = 2
    fp%naxes(1) = cube%nx
    fp%naxes(2) = cube%ny
    call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
    !
    call ftpprd(fp%fU, fp%group, fp%fpixel, cube%nx*cube%ny, arr_tau, fp%stat)
    !
    call ftpkyd(fp%fU, 'CDELT1', cube%dx,  fp%decimals, 'dx', fp%stat)
    call ftpkyd(fp%fU, 'CDELT2', cube%dy,  fp%decimals, 'dy', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL1', cube%xmin,  fp%decimals, 'xmin', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL2', cube%ymin,  fp%decimals, 'ymin', fp%stat)
    call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
    call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
    call ftpkys(fp%fU, 'ExtName', 'TauMap', 'peak values', fp%stat)
  end if
  !
  ! Second extension: integrated map
  call ftcrhd(fp%fU, fp%stat)
  fp%naxis = 2
  fp%naxes(1) = cube%nx
  fp%naxes(2) = cube%ny
  call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
  !
  call ftpprd(fp%fU, fp%group, fp%fpixel, cube%nx*cube%ny, sum(cube%val, 3) * cube%df, fp%stat)
  !
  call ftpkyd(fp%fU, 'CDELT1', cube%dx,  fp%decimals, 'dx', fp%stat)
  call ftpkyd(fp%fU, 'CDELT2', cube%dy,  fp%decimals, 'dy', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL1', cube%xmin,  fp%decimals, 'xmin', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL2', cube%ymin,  fp%decimals, 'ymin', fp%stat)
  call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
  call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
  call ftpkys(fp%fU, 'ExtName', 'IntMap', 'Int(I, nu)', fp%stat)
  !
  ! Third extension: upper column density
  if (is_l .and. present(Ncol_up)) then
    call ftcrhd(fp%fU, fp%stat)
    fp%naxis = 2
    fp%naxes(1) = cube%nx
    fp%naxes(2) = cube%ny
    call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
    !
    call ftpprd(fp%fU, fp%group, fp%fpixel, cube%nx*cube%ny, Ncol_up, fp%stat)
    !
    call ftpkyd(fp%fU, 'CDELT1', cube%dx,  fp%decimals, 'dx', fp%stat)
    call ftpkyd(fp%fU, 'CDELT2', cube%dy,  fp%decimals, 'dy', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL1', cube%xmin,  fp%decimals, 'xmin', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL2', cube%ymin,  fp%decimals, 'ymin', fp%stat)
    call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
    call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
    call ftpkys(fp%fU, 'ExtName', 'ColumnDensityUp', 'cm-2', fp%stat)
  end if
  !
  ! Fourth extension: lower column density
  if (is_l .and. present(Ncol_low)) then
    call ftcrhd(fp%fU, fp%stat)
    fp%naxis = 2
    fp%naxes(1) = cube%nx
    fp%naxes(2) = cube%ny
    call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
    !
    call ftpprd(fp%fU, fp%group, fp%fpixel, cube%nx*cube%ny, Ncol_low, fp%stat)
    !
    call ftpkyd(fp%fU, 'CDELT1', cube%dx,  fp%decimals, 'dx', fp%stat)
    call ftpkyd(fp%fU, 'CDELT2', cube%dy,  fp%decimals, 'dy', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
    call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL1', cube%xmin,  fp%decimals, 'xmin', fp%stat)
    call ftpkyd(fp%fU, 'CRVAL2', cube%ymin,  fp%decimals, 'ymin', fp%stat)
    call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
    call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
    call ftpkys(fp%fU, 'ExtName', 'ColumnDensityLow', 'cm-2', fp%stat)
  end if
  !
  ! Fifth extension: spectrum integrated over the whole region
  call ftcrhd(fp%fU, fp%stat)
  fp%naxis = 2
  fp%naxes(1) = cube%nf
  fp%naxes(2) = 1
  call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
  !
  call ftpprd(fp%fU, fp%group, fp%fpixel, cube%nf, vec_flux, fp%stat)
  !
  call ftpkyd(fp%fU, 'CDELT1', cube%df,  fp%decimals, 'df', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,  fp%decimals, 'k0', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL1', cube%fmin, fp%decimals, 'fmin', fp%stat)
  call ftpkys(fp%fU, 'CTYPE1', 'F', 'Hz', fp%stat)
  if (is_l) then
    call ftpkyd(fp%fU, 'F0', cube%f0, fp%decimals, 'Hz', fp%stat)
  end if
  !
  call ftpkyd(fp%fU, 'CDELT2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkys(fp%fU, 'CTYPE2', 'null', 'To make ds9 work.', fp%stat)
  call ftpkys(fp%fU, 'ExtName', 'FluxSpec', 'jy', fp%stat)
  !
  call ftclos(fp%fU, fp%stat)
  call ftfiou(fp%fU, fp%stat)
end subroutine save_cube_to_fits


subroutine save_cube_to_fits_spec_only(filename, cube, vec_flux, arr_tau, Ncol_up, Ncol_low, is_line)
  use my_timer
  character(len=*), intent(in) :: filename
  type(type_cube), intent(in) :: cube
  double precision, dimension(:), intent(in) :: vec_flux
  double precision, dimension(:,:), intent(in), optional :: arr_tau, Ncol_up, Ncol_low
  logical, intent(in), optional :: is_line
  type(type_fits_par) :: fp
  logical is_l
  !
  double precision dv, vmax
  type(date_time) a_date_time
  !
  if (present(is_line)) then
    is_l = is_line
  else
    is_l = .true.
  end if
  !
  fp%filename = trim(filename)
  !
  fp%stat = 0
  fp%blocksize = 1
  fp%pcount = 0
  fp%gcount = 1
  fp%group=1
  fp%fpixel=1
  fp%decimals = 16
  fp%author = 'Fujun Du (fdu@umich.edu)'
  fp%user = ''
  fp%simple=.true.
  fp%extend=.true.
  fp%bitpix=-64 ! double
  !
  call ftgiou(fp%fU, fp%stat)
  call ftinit(fp%fU, fp%filename, fp%blocksize, fp%stat)
  !
  fp%naxis=2
  fp%naxes(1) = cube%nf
  fp%naxes(2) = 1
  fp%nelements = cube%nf
  !
  call ftphpr(fp%fU, fp%simple, fp%bitpix, fp%naxis, fp%naxes, &
              fp%pcount, fp%gcount, fp%extend, fp%stat)
  !
  ! The cube.
  call ftpprd(fp%fU, fp%group, fp%fpixel, fp%nelements, vec_flux, fp%stat)
  !
  call ftpkyd(fp%fU, 'CDELT1', cube%df,  fp%decimals, 'df', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,  fp%decimals, 'k0', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL1', cube%fmin, fp%decimals, 'fmin', fp%stat)
  call ftpkys(fp%fU, 'CTYPE1', 'F', 'Hz', fp%stat)
  !
  call ftpkyd(fp%fU, 'CDELT2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkyd(fp%fU, 'CRPIX2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkyd(fp%fU, 'CRVAL2', 0D0, fp%decimals, 'null', fp%stat)
  call ftpkys(fp%fU, 'CTYPE2', 'null', 'To make ds9 work.', fp%stat)
  call ftpkys(fp%fU, 'ExtName', 'FluxSpec', 'jy', fp%stat)
  !
  call ftpkyd(fp%fU, 'Dist',  raytracing_conf%dist, fp%decimals, 'pc', fp%stat)
  call ftpkyd(fp%fU, 'Theta', cube%view_theta, fp%decimals, 'deg', fp%stat)
  call ftpkyd(fp%fU, 'MaxFlux', maxval(vec_flux),  fp%decimals, 'jy', fp%stat)
  if (is_l .and. present(arr_tau)) then
    call ftpkyd(fp%fU, 'MaxTau',  maxval(arr_tau),   fp%decimals, '', fp%stat)
  end if
  if (is_l) then
    !call ftpkyj(fp%fU, 'Itr',   cube%itr,    'trans num', fp%stat)
    call ftpkyd(fp%fU, 'F0',    cube%f0,     fp%decimals, 'Hz', fp%stat)
    call ftpkyd(fp%fU, 'lam0',  cube%rapar%lambda, fp%decimals, 'Angstrom', fp%stat)
    call ftpkyd(fp%fU, 'Eup',   cube%rapar%Eup,  fp%decimals, 'K', fp%stat)
    call ftpkyd(fp%fU, 'Elow',  cube%rapar%Elow, fp%decimals, 'K', fp%stat)
    !call ftpkyj(fp%fU, 'iup',   cube%rapar%iup,  '', fp%stat)
    !call ftpkyj(fp%fU, 'ilow',  cube%rapar%ilow, '', fp%stat)
    call ftpkyd(fp%fU, 'Aul',   cube%rapar%Aul,  fp%decimals, 's-1', fp%stat)
    call ftpkyd(fp%fU, 'Bul',   cube%rapar%Bul,  fp%decimals, '', fp%stat)
    call ftpkyd(fp%fU, 'Blu',   cube%rapar%Blu,  fp%decimals, '', fp%stat)
    call ftpkys(fp%fU, 'Qnum',  trim(cube%rapar%qnum), '', fp%stat)
  end if
  !
  !call ftpkys(fp%fU, 'Author', fp%author, '', fp%stat)
  call ftpkys(fp%fU, 'User',   fp%user,   '', fp%stat)
  call ftpkys(fp%fU, 'SavedAt', trim(a_date_time%date_time_str()), '', fp%stat)
  !
  call ftclos(fp%fU, fp%stat)
  call ftfiou(fp%fU, fp%stat)
end subroutine save_cube_to_fits_spec_only




subroutine line_excitation_do
  type(type_cell), pointer :: c
  integer i
  write(*, '(A/)') 'Doing energy level excitation calculation.'
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    call allocate_local_cont_lut(c)
    call make_local_cont_lut(c)
    call do_exc_calc(c)
  end do
end subroutine line_excitation_do


subroutine continuum_tran_prep
  integer i
  write(*, '(/A)') 'Preparing for the continuum radiative transfer.'
  do i=1, leaves%nlen
    call allocate_local_cont_lut(leaves%list(i)%p)
    call make_local_cont_lut(leaves%list(i)%p)
  end do
end subroutine continuum_tran_prep



subroutine line_tran_prep
  write(*, '(/A)') 'Preparing for the line radiative transfer.'
  call load_exc_molecule
  !
  if (.not. mole_exc%conf%useLTE) then
    call init_statistic_sol(mole_exc%p%n_level, raytracing_conf%solve_method)
  end if
end subroutine line_tran_prep



subroutine load_exc_molecule
  integer i, i0, i1, j
  character(len=const_len_species_name) str, str1
  integer, dimension(:), allocatable :: itmp, itmp1
  double precision freq, en, Qpart, Ttmp
  integer iup, ilow, fU
  logical in_freq_window
  !
  mole_exc%conf = raytracing_conf
  allocate(mole_exc%p)
  !
  mole_exc%p%abundance_factor = mole_exc%conf%abundance_factor
  !
  select case(mole_exc%conf%line_database)
    case ('lamda')
      call load_moldata_LAMDA(&
        combine_dir_filename(mole_exc%conf%dirname_mol_data, &
        mole_exc%conf%fname_mol_data), mole_exc%p)
    case ('hitran')
      call load_hitran_mol(&
        mole_exc%conf%dirname_mol_data, &
        mole_exc%conf%mole_name, &
        mole_exc%p) !, &
        !Elow_range=(/0D0, mole_exc%conf%E_max*5D0/)
        !lam_range=(/minval(mole_exc%conf%lam_mins(1:mole_exc%conf%nlam_window)), &
        !            maxval(mole_exc%conf%lam_maxs(1:mole_exc%conf%nlam_window))/), &
        !Elow_range=(/mole_exc%conf%E_min, mole_exc%conf%E_max/)
    case ('cdms')
      call load_cdms_mol(mole_exc%conf%dirname_mol_data, &
                         mole_exc%conf%fname_mol_data, &
                         mole_exc%conf%fname_parti_data, &
                         mole_exc%p)
      mole_exc%p%name_molecule = mole_exc%conf%mole_name
    case default
      write(*, '(A)') 'Unknown line excitation data format!'
      write(*, '(A)') 'Currently only support:'
      write(*, '(A)') '"lambda", "hitran", and "cdms" (case sensitive).'
      stop
  end select
  !
  mole_exc%p%name_surrogate = mole_exc%conf%mole_name_surrogate
  !
  mole_exc%p%iType = -1
  !
  i = index(mole_exc%p%name_molecule, '(')
  if (i .eq. 0) then
    str = mole_exc%p%name_molecule
    mole_exc%p%iType = 0
    str1 = ''
  else
    str = mole_exc%p%name_molecule(1:(i-1))
    i = index(mole_exc%p%name_molecule, 'ortho')
    if (i .ne. 0) then
      mole_exc%p%iType = 1
      str1 = 'ortho'
    else
      i = index(mole_exc%p%name_molecule, 'para')
      if (i .ne. 0) then
        mole_exc%p%iType = 2
        str1 = 'para'
      end if
    end if
  end if
  !
  if (len_trim(mole_exc%p%name_surrogate) .ge. 1) then
    str = adjustl(mole_exc%p%name_surrogate)
  end if
  !
  mole_exc%p%iSpe = -1
  !
  do i=1, chem_species%nSpecies
    if (str .eq. chem_species%names(i)) then
      mole_exc%p%iSpe = i
      exit
    end if
  end do
  if ((mole_exc%p%iSpe .eq. -1) .or. (mole_exc%p%iType .eq. -1)) then
    write(*, '(A)') 'In load_exc_molecule:'
    write(*, '(A)') 'Unidentified molecule name and/or type:'
    write(*, '(3A)') mole_exc%p%name_molecule, ' ', mole_exc%p%name_surrogate
    write(*, '(A)') 'In file:'
    write(*, '(A)') combine_dir_filename( &
      mole_exc%conf%dirname_mol_data, &
      mole_exc%conf%fname_mol_data)
    stop
  end if
  write(*, '(A, 2A16)') 'Molecule: ', trim(str), str1
  write(*, '(A, I6)') 'Total number of levels: ', mole_exc%p%n_level
  write(*, '(A, I6)') 'Total number of radiative transitions: ', mole_exc%p%rad_data%n_transition
  write(*, '(A, I6)') 'Total number of collisional partners: ', mole_exc%p%colli_data%n_partner
  do i=1, mole_exc%p%colli_data%n_partner
    write(*, '(I2, 2X, 2A)') i, 'Partner name: ', mole_exc%p%colli_data%list(i)%name_partner
    write(*, '(I2, 2X, A, I6)') i, 'Total number of collisional transitions: ', &
      mole_exc%p%colli_data%list(i)%n_transition
    write(*, '(I2, 2X, A, I6)') i, 'Total number of collisional temperatures: ', &
      mole_exc%p%colli_data%list(i)%n_T
  end do
  !write(*,*)
  !write(*, '(A, 2ES12.4)') 'Frequency range to consider: ', &
  !     mole_exc%conf%freq_min, mole_exc%conf%freq_max
  !
  call openFileSequentialWrite(fU, &
       combine_dir_filename(dir_name_log, 'energy_levels_all.dat'), &
       999, getu=1)
  do i=2,11
    Ttmp = 1D1**(dble(i)*0.3D0)
    mole_exc%p%f_occupation = mole_exc%p%level_list%weight * exp(-mole_exc%p%level_list%energy / Ttmp)
    Qpart = sum(mole_exc%p%f_occupation)
    write(fU, '(A, ES12.2, A, ES16.6)') 'Partition function for T = ', Ttmp, ' K: ', Qpart
  end do
  write(fU, '(A10, 3A19)') 'Num', 'E(K)', 'g', 'f(T=300K)'
  mole_exc%p%f_occupation = mole_exc%p%level_list%weight * exp(-mole_exc%p%level_list%energy / 3D2)
  do i=1, mole_exc%p%n_level
    write(fU, '(I10, 3ES19.10)') i, mole_exc%p%level_list(i)%energy, &
        mole_exc%p%level_list(i)%weight, mole_exc%p%f_occupation(i)/Qpart
  end do
  close(fU)
  !
  allocate(itmp(mole_exc%p%n_level), &
           itmp1(mole_exc%p%rad_data%n_transition), &
           mole_exc%ilv_reverse(mole_exc%p%n_level))
  mole_exc%ilv_reverse = 0
  i0 = 0
  i1 = 0
  do i=1, mole_exc%p%rad_data%n_transition
    freq = mole_exc%p%rad_data%list(i)%freq
    en   = mole_exc%p%rad_data%list(i)%Eup
    !
    in_freq_window = .false.
    do j=1, mole_exc%conf%nfreq_window
      if ((mole_exc%conf%freq_mins(j) .le. freq) .and. &
          (freq .le. mole_exc%conf%freq_maxs(j))) then
        in_freq_window = .true.
        exit
      end if
    end do
    !
    if (in_freq_window .and. &
        (en .ge. mole_exc%conf%E_min) .and. &
        (en .le. mole_exc%conf%E_max)) then
      i1 = i1 + 1
      itmp1(i1) = i
      iup = mole_exc%p%rad_data%list(i)%iup
      ilow = mole_exc%p%rad_data%list(i)%ilow
      if (.not. is_in_list_int(ilow, i0, itmp(1:i0))) then
        i0 = i0 + 1
        itmp(i0) = ilow
        mole_exc%ilv_reverse(ilow) = i0
      end if
      if (.not. is_in_list_int(iup, i0, itmp(1:i0))) then
        i0 = i0 + 1
        itmp(i0) = iup
        mole_exc%ilv_reverse(iup) = i0
      end if
    end if
  end do
  !
  mole_exc%nlevel_keep = i0
  mole_exc%ntran_keep  = i1
  allocate(mole_exc%ilv_keep(i0), &
           mole_exc%itr_keep(i1))
  mole_exc%ilv_keep = itmp(1:i0)
  mole_exc%itr_keep = itmp1(1:i1)
  !
  ! 2014-07-03 Thu 16:02:26
  ! Sort the transitions in increasing frequency
  do i=1, mole_exc%ntran_keep
    do j=1, i-1
      i0 = mole_exc%itr_keep(i)
      i1 = mole_exc%itr_keep(j)
      if (mole_exc%p%rad_data%list(i0)%freq .lt. &
          mole_exc%p%rad_data%list(i1)%freq) then
        mole_exc%itr_keep(i) = i1
        mole_exc%itr_keep(j) = i0
      end if
    end do
  end do
  !
  deallocate(itmp, itmp1)
  write(*, '(A, I6)') 'Number of levels to keep:', mole_exc%nlevel_keep
  write(*, '(A, I6)') 'Number of transitions to keep:', mole_exc%ntran_keep
  !do i=1, mole_exc%nlevel_keep
  !  i0 = mole_exc%ilv_keep(i)
  !  write(*, '(I4, ES12.4)') i, mole_exc%p%level_list(i0)%energy
  !end do
  write(*,*)
end subroutine load_exc_molecule


subroutine set_using_mole_params(mole, c)
  type(type_cell), intent(in), pointer :: c
  type(type_molecule_energy_set), intent(inout), pointer :: mole
  select case (mole%iType)
  case (0)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe)
  case (1)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe) * 0.75D0
  case (2)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe) * 0.25D0
  case default
    write(*, '(A)') 'In set_using_mole_params:'
    write(*, '(A, I4)') 'Unknown molecule type: ', mole%iType
    write(*, '(A)') 'Will use the full abundance.'
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe)
  end select
  !
  mole%density_mol = mole%density_mol * mole%abundance_factor
  ! Ad hoc
  !!mole%density_mol = mole%density_mol * &
  !!  calibrate_HD_abundance(c%abundances(mole%iSpe), &
  !!                         c%par%f_selfshielding_toStar_H2, &
  !!                         c%par%Av_toStar, &
  !!                         mole%abundance_factor)
  !write(*, '(A, ES12.4)') 'Applying abundance modification factor:', mole%abundance_factor
  !
  mole%Tkin = c%par%Tgas
  !mole%dv = c%par%velo_width_turb
  ! 2014-06-10 Tue 02:36:34
  mole%dv = sqrt(phy_kBoltzmann_CGS*c%par%Tgas/(chem_species%mass_num(mole%iSpe) * phy_mProton_CGS))
  mole%length_scale = c%par%coherent_length
end subroutine set_using_mole_params


pure function calibrate_HD_abundance(X_H2, f_H2, Av, D2H) result(ratio)
  ! Return HD/H2
  double precision ratio
  double precision, intent(in) :: X_H2, f_H2, Av, D2H
  double precision r
  if (Av .ge. 2D0) then
    ratio = 2D0 * D2H
  else
    r = max(0D0, (0.5D0 / X_H2 - 1D0) * sqrt(2D0) / f_H2)
    ratio = 2D0 * D2H * (0.5D0 / X_H2) / (1D0 + r)
  end if
end function calibrate_HD_abundance



subroutine do_exc_calc(c)
  type(type_cell), intent(inout), pointer :: c
  integer i, iLow, iUp, ic
  double precision Qpart, tmp
  !
  mol_sta_sol => mole_exc%p
  !
  call set_using_mole_params(mol_sta_sol, c)
  !
  mol_sta_sol%f_occupation = mol_sta_sol%level_list%weight * &
      exp(-mol_sta_sol%level_list%energy / mol_sta_sol%Tkin)
  Qpart = sum(mol_sta_sol%f_occupation)
  do i=1, mol_sta_sol%n_level
    mol_sta_sol%f_occupation(i) = &
      mol_sta_sol%f_occupation(i) / Qpart
  end do
  !
  if (.not. mole_exc%conf%useLTE) then
    do i=1, mol_sta_sol%colli_data%n_partner
      select case (mol_sta_sol%colli_data%list(i)%name_partner)
      case ('H2')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_H2
      case ('o-H2')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          0.75D0 * c%par%n_gas * c%par%X_H2
      case ('p-H2')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          0.25D0 * c%par%n_gas * c%par%X_H2
      case ('H')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_HI
      case ('H+')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_Hplus
      case ('e')
        mol_sta_sol%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_E
      case default
        write(*, '(A)') 'In do_exc_calc:'
        write(*, '(A)') 'Unknown collision partner:'
        write(*, '(A)') mol_sta_sol%colli_data%list(i)%name_partner
        write(*, '(A)') 'Will use zero abundance for this partner.'
        mol_sta_sol%colli_data%list(i)%dens_partner = 0D0
      end select
    end do
    !
    current_cell_ptr => c
    !
    call statistic_equil_solve
  end if
  !
  if (.not. allocated(c%focc)) then
    allocate(c%focc)
    c%focc%nlevels = mole_exc%nlevel_keep
    allocate(c%focc%vals(c%focc%nlevels))
  end if
  !
  c%focc%vals = mol_sta_sol%f_occupation(mole_exc%ilv_keep)
  !
  !if (mole_exc%conf%useLTE) then
  !  ic = 0
  !  do i=1, mole_exc%nlevel_keep
  !    if (i .gt. 100) then
  !      iLow = mole_exc%ilv_keep(i-100)
  !      iUp  = mole_exc%ilv_keep(i)
  !      tmp = mol_sta_sol%level_list(iUp)%weight / &
  !            mol_sta_sol%level_list(iLow)%weight * &
  !            exp(-(mol_sta_sol%level_list(iUp)%energy - &
  !            mol_sta_sol%level_list(iLow)%energy) / mol_sta_sol%Tkin)
  !      if (abs(c%focc%vals(i)/c%focc%vals(i-100)/tmp-1D0) .ge. 1D-6) then
  !        write(*, '(A, 3I6, 8ES16.6)') 'Population error: ', &
  !            i, iLow, iUp, &
  !            c%focc%vals(i), c%focc%vals(i-100), &
  !            c%focc%vals(i)/c%focc%vals(i-100)/tmp, &
  !            mol_sta_sol%level_list(iUp)%energy, &
  !            mol_sta_sol%level_list(iLow)%energy, &
  !            mol_sta_sol%level_list(iUp)%weight, &
  !            mol_sta_sol%level_list(iLow)%weight, &
  !            c%par%Tgas
  !        ic = ic + 1
  !      end if
  !    end if
  !  end do
  !  if (ic .gt. 0) then
  !    stop
  !  end if
  !end if
  !
  !nullify(mol_sta_sol)
end subroutine do_exc_calc



subroutine make_local_cont_lut(c)
  type(type_cell), intent(inout), pointer :: c
  integer i
  double precision dlam, lam
  !
  do i=1, c%cont_lut%n
    if (i .lt. c%cont_lut%n) then
      dlam = dust_0%lam(i+1) - dust_0%lam(i)
      lam = (dust_0%lam(i+1) + dust_0%lam(i)) * 0.5D0
      ! Energy per unit area per unit frequency per second per sqradian
      c%cont_lut%J(i) = c%optical%flux(i) &
        / dlam * lam * lam * phy_Angstrom2cm / phy_SpeedOfLight_CGS &
        / (4D0 * phy_Pi)
    else
      c%cont_lut%J(i) = c%cont_lut%J(i-1)
    end if
  end do
end subroutine make_local_cont_lut



subroutine allocate_local_cont_lut(c)
  type(type_cell), intent(inout), pointer :: c
  !
  if (.not. allocated(c%cont_lut)) then
    allocate(c%cont_lut)
  end if
  !
  if (.not. allocated(c%cont_lut%J)) then
    c%cont_lut%n = dust_0%n
    allocate(c%cont_lut%J(dust_0%n))
  end if
end subroutine allocate_local_cont_lut


function get_ave_temperature() result(Tave)
  double precision Tave, m
  integer i
  Tave = 0D0
  m = 0D0
  do i=1, leaves%nlen
    m    = m    + leaves%list(i)%p%par%mgas_cell
    Tave = Tave + leaves%list(i)%p%par%mgas_cell * leaves%list(i)%p%par%Tgas
  end do
  Tave = Tave / m
end function get_ave_temperature


subroutine do_adjust_yup_ylow_nonLTE(yup, ylow, gu, gl, nu, Tgas, n_H, n_crit, Jnu)
  double precision, intent(inout) :: yup, ylow
  double precision, intent(in) :: gu, gl, nu, Tgas, n_H, n_crit, Jnu
  double precision r1, r2, r3, r4, r, t
  r1 = gu / gl
  r2 = n_H / n_crit
  r3 = Jnu / ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * nu**3)
  r4 = exp(-phy_hPlanck_CGS*nu/(phy_kBoltzmann_CGS*Tgas))
  r = r1 * (r3 + r2 * r4) / (1D0 + r3 + r2)
  t = yup + ylow
  yup  = t * r / (1D0 + r)
  ylow = t     / (1D0 + r)
end subroutine do_adjust_yup_ylow_nonLTE


end module ray_tracing
