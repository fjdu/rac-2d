module montecarlo

use phy_const
use trivials
use grid
use data_struct
use ray_propagating

implicit none

type :: type_dust_optical_collection
  integer n
  type(type_optical_property), dimension(:), allocatable :: list
end type type_dust_optical_collection


type :: type_dust_lut_collection
  integer n
  type(type_LUT_Tdust), dimension(:), allocatable :: list
end type type_dust_lut_collection


type(type_optical_property) dust_0, HI_0, water_0
type(type_stellar_params) a_star
type(type_montecarlo_config) mc_conf

type(type_distribution_table) p4lam
type(type_global_material_collection) opmaterials

type(type_dust_optical_collection) dusts
type(type_dust_lut_collection) luts

type(type_photon_collector) collector

double precision, dimension(2), parameter :: lam_range_Xray = (/0.1D0, 1D2/)
double precision, dimension(2), parameter :: lam_range_UV   = (/9D2, 3D3/)
double precision, dimension(2), parameter :: lam_range_LyA  = (/1210D0, 1220D0/)
double precision, dimension(2), parameter :: lam_range_Vis  = (/3D3, 8D3/)
double precision, dimension(2), parameter :: lam_range_NIR  = (/8D3, 5D4/)
double precision, dimension(2), parameter :: lam_range_MIR  = (/5D4, 3D5/)
double precision, dimension(2), parameter :: lam_range_FIR  = (/3D5, 2D6/)
double precision, dimension(2), parameter :: lam_range_LyA_ext  = (/1100D0, 1300D0/)

integer, parameter :: icl_HI   = 1, &
                      icl_H2O  = 2, &
                      icl_dust = 3, &
                      ncl_nondust = 2

namelist /montecarlo_configure/ mc_conf


contains


!subroutine save_photon(ph, mc)
!  type(type_photon_packet), intent(in) :: ph
!  type(type_montecarlo_config), intent(in) :: mc
!  if (mc%savephoton) then
!    write(mc%fU, '(8ES19.10, I10)') & 
!        ph%ray%x, ph%ray%y, ph%ray%z, &
!        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
!        ph%lam, ph%en, ph%e_count
!  end if
!end subroutine save_photon


subroutine make_luts
  integer i
  luts%n = dusts%n
  allocate(luts%list(luts%n))
  do i=1, luts%n
    call make_LUT_Tdust(dusts%list(i), luts%list(i), mc_conf%nlen_lut, &
        mc_conf%TdustMin, mc_conf%TdustMax)
  end do
end subroutine make_luts


subroutine get_mc_stellar_par(star, mc)
  type(type_stellar_params), intent(inout) :: star
  type(type_montecarlo_config), intent(inout) :: mc
  !
  if (mc_conf%ph_init_symmetric) then
    ! - Since each phantom cell with z<0 is identified with a corresponding cell
    !   with z>0, in this case we need to reduce the received radiation energy of
    !    each cell by a factor of 2.
    ! - A really robust implementation should have a complete disk model with
    !   both the z>=0 and z<0 parts, but that requires too much memory.
    ! - Or we could use a complete disk model only for the radiative transfer,
    !   while for chemistry and for storing the optical data we only use one half.
    star%vals = star%vals0 * (mc%maxw - mc%minw) / 4D0
    star%lumi = star%lumi0 * (mc%maxw - mc%minw) / 4D0
  else
    star%vals = star%vals0 * (mc%maxw - mc%minw) / 2D0
    star%lumi = star%lumi0 * (mc%maxw - mc%minw) / 2D0
  end if
  !
  mc%eph = star%lumi / dble(mc%nph)
  !
  write(*,'(A, ES16.6, A)') 'Stellar luminosity within (minw,maxw): ', &
    star%lumi, ' erg s-1.'
  write(*,'(A, ES16.6, A)') 'Lumi per photon: ', mc%eph, ' erg s-1.'
end subroutine get_mc_stellar_par



subroutine align_optical_data
  ! The imported optical data for dust, HI, and water are not aligned.
  ! So here I will align them to the same lambda vector.
  double precision, dimension(:), allocatable :: v, v1
  integer i, j, n, n1, n_using
  double precision lam, en, mu_median
  !
  n1 = HI_0%n + water_0%n
  n = dust_0%n + HI_0%n + water_0%n
  allocate(v1(n1), v(n))
  call merge_vec(water_0%n, water_0%lam, HI_0%n, HI_0%lam, n1, v1, n_using)
  call merge_vec(dust_0%n, dust_0%lam, n1, v1, n, v, n_using)
  !
  call reassign_optical(dust_0, n, v)
  do i=1, dusts%n
    call reassign_optical(dusts%list(i), n, v)
    !
    do j=1, dusts%list(i)%n
      lam = dusts%list(i)%lam(j)
      if ((lam .lt. lam_range_Xray(1)) .or. (lam .gt. lam_range_Xray(2))) then
        cycle
      end if
      en = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
           / (lam * 1D-8) / phy_eV2erg / 1D3 ! in keV
      ! Draine 2003, equation 9
      ! 0.1 deg = 360 arcsec
      mu_median = cos(min(1D0, 0.1D0/180D0 / en) * phy_Pi)
      dusts%list(i)%g(j) = 1D0 - mu_median / sqrt(2D0)
    end do
    !
  end do
  call reassign_optical(HI_0, n, v)
  call reassign_optical(water_0, n, v)
  !
  deallocate(v1, v)
end subroutine align_optical_data



pure subroutine reassign_optical(op, n, v)
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
end subroutine reassign_optical



subroutine make_global_coll
  !
  opmaterials%ntype = ncl_nondust + dusts%n
  !
  allocate(opmaterials%list(opmaterials%ntype))
  !
  opmaterials%list(icl_HI) = HI_0
  opmaterials%list(icl_H2O) = water_0
  opmaterials%list(icl_dust:) = dusts%list
  !
end subroutine make_global_coll



subroutine make_Xray_abs_sca(c)
  ! Make the local X-ray absorption and scattering parameters.
  ! The cross sections are on the per H atom basis.
  use load_bethell_xray_cross
  type(type_cell), intent(in), pointer :: c
  double precision lam, en
  integer i, n
  !
  n = HI_0%n
  !
  ! opmaterials acts as a global storage here.
  !
  if (.not. allocated(opmaterials%Xray_gas_abs)) then
    allocate(opmaterials%Xray_gas_abs(n), &
             opmaterials%Xray_gas_sca(n), &
             opmaterials%Xray_dus_abs(n), &
             opmaterials%Xray_dus_sca(n))
  end if
  do i=1, n
    lam = HI_0%lam(i)
    if ((lam .lt. lam_range_Xray(1)) .or. &
        (lam .gt. lam_range_Xray(2))) then
      opmaterials%Xray_gas_abs(i) = 0D0
      opmaterials%Xray_gas_sca(i) = 0D0
      opmaterials%Xray_dus_abs(i) = 0D0
      opmaterials%Xray_dus_sca(i) = 0D0
      cycle
    end if
    !
    en = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
         / (lam * 1D-8) / phy_eV2erg / 1D3 ! in keV
    !
    opmaterials%Xray_gas_abs(i) = sigma_Xray_Bethell_gas(en)
    ! H + He
    opmaterials%Xray_gas_sca(i) = phy_ThomsonScatterCross_CGS * (1D0 + 1D0/6D0)
    ! Assume the Thomson scattering is isotropic, though it is not.
    !
    opmaterials%Xray_dus_abs(i) = &
      sigma_Xray_Bethell_dust(en, c%par%dust_depletion, &
        c%par%ratioDust2HnucNum, c%par%GrainRadius_CGS)
    ! The scattering cross section is an analytical fitting based on the table
    ! of Draine 2003.
    opmaterials%Xray_dus_sca(i) = c%par%dust_depletion * &
                                  1.3D-22 / (en**1.8D0 + 0.4D0)
  end do
  !
end subroutine make_Xray_abs_sca




subroutine update_gl_optical_OTF(T)
  ! The HI scattering dependends on temperature, so it needs to be
  ! recalculated for each cell.
  use voigt
  double precision, intent(in) :: T
  integer i
  double precision dnu_th, a, nu, x, tmp
  !
  dnu_th = phy_LyAlpha_nu0 * &
    sqrt(8D0*phy_kBoltzmann_SI*T/phy_Pi/phy_mProton_SI) &
      / phy_SpeedOfLight_SI
  a = phy_LyAlpha_dnul / (2D0 * dnu_th)
  tmp = phy_LyAlpha_f12 * sqrt(phy_Pi) * &
        phy_electronClassicalRadius_CGS * &
        phy_SpeedOfLight_CGS / dnu_th
  !
  do i=1, opmaterials%list(icl_HI)%n
    if ((HI_0%lam(i) .ge. lam_range_LyA_ext(1)) .and. &
        (HI_0%lam(i) .le. lam_range_LyA_ext(2))) then
      nu = phy_SpeedOfLight_SI / (HI_0%lam(i) * 1D-10)
      x = abs((nu - phy_LyAlpha_nu0)) / dnu_th
      HI_0%sc(i) = tmp * max(0D0, voigt_scalar(x, a))
      !write(*, '(F10.5, I5, 2ES15.6)') T, i, HI_0%lam(i), HI_0%sc(i)
    end if
  end do
  !
  ! Apply the changes.
  opmaterials%list(icl_HI) = HI_0
  !
end subroutine update_gl_optical_OTF



pure subroutine allocate_local_optics(c, ntype, nlam)
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in) :: ntype, nlam
  if (.not. c%using) then
    return
  end if
  if (allocated(c%optical%X)) then
    return
  end if
  c%optical%ntype = ntype*2
  c%optical%nlam = nlam
  ! ! <timestamp>2013-11-07 Thu 15:49:33</timestamp>
  ! An error in previous version?
  ! X should have the dimension equal to the number of optical materials, not
  ! twice this number.
  allocate(c%optical%X(ntype), &
           c%optical%acc(c%optical%nlam, c%optical%ntype), &
           c%optical%summed(c%optical%nlam), &
           c%optical%summed_ab(c%optical%nlam), &
           c%optical%summed_sc(c%optical%nlam), &
           c%optical%flux(c%optical%nlam), &
           c%optical%phc(c%optical%nlam), &
           c%optical%dir_wei(c%optical%nlam))
end subroutine allocate_local_optics



subroutine reset_local_optics(c)
  type(type_cell), intent(inout), pointer :: c
  integer i
  if (.not. c%using) then
    return
  end if
  !
  c%par%en_gains = 0D0
  c%par%en_prevs = 0D0
  c%par%en_gains_abso = 0D0
  c%par%kphs = 0D0
  !
  c%optical%cr_count = 0
  c%par%sc_count_HI = 0
  c%par%ab_count_water = 0
  c%par%ab_en_water = 0D0
  c%par%ab_count_dust = 0
  !
  do i=1, c%optical%nlam
    c%optical%dir_wei(i)%u = 0D0
    c%optical%dir_wei(i)%v = 0D0
    c%optical%dir_wei(i)%w = 0D0
    c%optical%flux(i) = 0D0
    c%optical%phc(i) = 0
  end do
  !
  call update_local_opticalX(c)
  call update_gl_optical_OTF(c%par%Tgas)
  call make_Xray_abs_sca(c)
  call make_local_optics(c, opmaterials)
end subroutine reset_local_optics



pure subroutine update_local_opticalX(c)
  type(type_cell), intent(inout), pointer :: c
  if (c%using) then
    c%optical%X(icl_HI)    = c%par%n_gas * c%par%X_HI
    c%optical%X(icl_H2O)   = c%par%n_gas * c%par%X_H2O
    c%optical%X(icl_dust:) = c%par%rho_dusts(1:c%par%ndustcompo) ! Mass density
  else
    c%optical%X = 0D0
  end if
end subroutine update_local_opticalX



subroutine montecarlo_do(mc, cstart)
  type(type_montecarlo_config), intent(inout) :: mc
  type(type_photon_packet) ph, ph0
  type(type_cell), intent(in), pointer :: cstart
  type(type_cell), pointer :: cthis
  integer(kind=LongInt) i, cPrema
  logical found, escaped, destructed
  double precision eph_acc
  integer fU
  !
  call init_random_seed
  !
  eph_acc = 0D0
  cPrema = 0
  !
  !if (mc%savephoton) then
  !  call openFileSequentialWrite(fU, &
  !    combine_dir_filename(mc%mc_dir_out, mc%fname_photons), 512)
  !end if
  !
  ph0%lam = a_star%lam(1)
  ph0%iSpec = 1
  !
  write(*,*)
  !
  i = 0
  do
    i = i + 1
    ! The exact number of photons will not be exactly equal to mc%nph, because
    ! each photon may have a different energy due to refinement at certain
    ! wavelength (Lya or any spectral feature you are interested in).
    !
    ! Emit a photon based on the stellar spectrum and the current photon
    ! wavelength.
    ! The wavelength of the photon to be emitted depends on the wavelength of
    ! the previous photon.
    call emit_a_photon(mc, ph0) ! ph%lam is in angstrom
    ph = ph0
    !
    eph_acc = eph_acc + ph%en
    if (eph_acc .gt. a_star%lumi) then
      mc_conf%icount = i
      exit
    end if
    !
    if (mod(i, 1000) .eq. 0) then
      write (*, &
        '(A, 4X, A, I12, 2X, "(", F0.4, "%)", &
          & 4X, A, ES14.6, I6, " of ", I6)') &
        CHAR(27)//'[A', "Monte Carlo...  Photon ", i, &
        eph_acc*1D2/a_star%lumi, "lam = ", ph%lam, ph%iSpec, a_star%n
    end if
    !
    ! Get the index for accessing the optical data
    ! Here we don't do the doppler shift, since for a single step that is not
    ! important, and furthermore, the photon velocity is perpendicular to the
    ! Kepplerian velocity.
    ph%iKap = get_idx_for_kappa(ph%lam, dust_0)
    if (ph%iKap .eq. 0) then ! No optical data for this lambda.
      !call save_photon(ph, mc)
      cycle
    end if
    !
    ! Enter the disk domain
    if (mc_conf%ph_init_symmetric) then
      call enter_the_domain_mirror(ph%ray, cstart, cthis, found)
    else
      call enter_the_domain(ph%ray, cstart, cthis, found)
    end if
    if (.not. found) then
      write(*, '(A, 6ES16.6, /, 4X, 4ES16.6)') 'Photon did not enter the domain:', &
                      ph%ray%x, ph%ray%y, ph%ray%z, &
                      ph%ray%vx, ph%ray%vy, ph%ray%vz, &
                      cstart%xmin, cstart%xmax, cstart%ymin, cstart%ymax
      cycle
    end if
    !
    ! Increase the crossing count
    cthis%optical%cr_count = cthis%optical%cr_count + 1
    !
    ! Track this photon until it is destroyed (and not reemitted) or has
    ! escaped the domain.
    call walk_scatter_absorb_reemit(ph, cthis, cstart, mc%nmax_cross, &
        escaped, destructed)
    !
    if (escaped) then ! Photon escaped from the disk domain
      !
      !call save_photon(ph, mc)
      !
      if (mc_conf%collect_photon) then
        call collect_photon_do(collector, ph)
      end if
      !
    else if (destructed) then
      ! Do nothing
    else
      cPrema = cPrema + 1
      write(*,'(I8, A/)') cPrema, 'Premature end of photon transport!'
    end if
  end do
  !
  !if (mc%savephoton) then
  !  close(fU)
  !end if
  !
end subroutine montecarlo_do



subroutine emit_a_photon(mc, ph)
  type(type_montecarlo_config), intent(inout) :: mc
  type(type_photon_packet), intent(inout) :: ph
  !
  if ((ph%lam .lt. lam_range_UV(1)) .or. &
      (ph%lam .gt. lam_range_UV(2)) .or. &
      mc%use_blackbody_star) then
    ! Not UV
    if ((ph%lam .lt. lam_range_Xray(1)) .or. &
        (ph%lam .gt. lam_range_Xray(2))) then
      ! Not X-ray
      ph%en = mc%eph
    else
      ! X-ray
      ph%en = mc%eph * mc%refine_Xray
    end if
  else if ((ph%lam .lt. lam_range_LyA(1)) .or. &
           (ph%lam .gt. lam_range_LyA(2))) then
    ! UV but not LyA
    ph%en = mc%eph * mc%refine_UV
  else
    ! LyA
    ph%en = mc%eph * mc%refine_LyA
  end if
  call get_next_lam(ph%lam, ph%iSpec, a_star, ph%en)
  ph%ray%x = mc%starpos_r
  ph%ray%y = 0D0
  ph%ray%z = mc%starpos_z
  ph%e_count = 0
  call get_emit_dir_uniform(ph%ray, mc%minw, mc%maxw)
end subroutine emit_a_photon



pure subroutine get_next_lam(lamthis, idx, star, eph)
  ! star%lam(idx) <= lamthis < star%lam(idx+1)
  integer, intent(inout) :: idx
  double precision, intent(inout) :: lamthis
  double precision, intent(in) :: eph
  type(type_stellar_params), intent(in) :: star
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



subroutine walk_scatter_absorb_reemit(ph, c, cstart, imax, &
  escaped, destructed)
  ! ph must be guaranteed to be inside c.
  ! An intersection between ph and c must exist, unless there is a
  ! numerical error.
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  type(type_cell), intent(in), pointer :: cstart
  integer(kind=LongInt), intent(in) :: imax
  logical, intent(out) :: escaped, destructed
  logical found, encountered
  type(type_cell), pointer :: cnext
  double precision tau_this, frac_abso
  integer(kind=LongInt) i
  double precision length, r, z, eps
  double precision rnd, tau, albedo, t
  integer dirtype
  integer itype, idust
  !
  escaped = .false.
  destructed = .false.
  !
  call random_number(rnd)
  tau = -log(rnd)
  !
  ! imax is usually set to a large number
  do i=1, imax
    ! Get the intersection between the photon ray and the boundary of the cell
    ! that this photon resides in
    if (mc_conf%ph_init_symmetric) then
      call calc_intersection_ray_cell_mirror(ph%ray, c, &
        length, r, z, eps, found, dirtype)
    else
      call calc_intersection_ray_cell(ph%ray, c, &
        length, r, z, eps, found, dirtype)
    end if
    if (.not. found) then
      write(*,'(A, I6, 9ES15.6/)') 'ph does not cross c: ', &
        i, sqrt(ph%ray%x**2+ph%ray%y**2), ph%ray%z, &
        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
        c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    !
    ! Get the index of the lambda in the array of optical data.  Even if lambda
    ! is not changed, the index may change from cell to cell due to change in
    ! the cell velocity.
    ! Note that ph%lam is always the photon wavelength as seen in the global
    ! rest frame, namely the irrotational frame centered on the star, hence we
    ! need to first doppler-shift the lambda to the local rest frame of the
    ! cell, then find the index.
    ! 
    ph%iKap = get_idx_for_kappa( &
        get_doppler_lam(a_star%mass, ph%lam, ph%ray), &
        dust_0)
    !
    if ((ph%iKap .gt. 0) .and. c%using) then
      tau_this = c%optical%summed(ph%iKap) * length * phy_AU2cm
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
    !
    if (c%using) then
      ! Todo
      if (c%optical%summed(ph%iKap) .gt. 0D0) then
        albedo = c%optical%summed_sc(ph%iKap) / &
                 c%optical%summed(ph%iKap)
        frac_abso = tau2frac(tau_this) * (1D0 - albedo)
      else
        frac_abso = 0D0
      end if
      !
      c%par%en_prevs = c%par%en_gains
      ! Distribute the energy into different dust species.
      ! Maybe better to use absorption cross section as the weight?
      c%par%en_gains = c%par%en_gains + frac_abso * ph%en * c%par%abso_wei
      !
      c%optical%flux(ph%iKap) = c%optical%flux(ph%iKap) + length * ph%en
      c%optical%phc(ph%iKap) = c%optical%phc(ph%iKap) + 1
      !
      ! Project the photon direction to local radial frame, and save this
      ! information for describing the radiation field.
      t = sqrt(ph%ray%x**2 + ph%ray%y**2) + 1D-100
      c%optical%dir_wei(ph%iKap)%u = c%optical%dir_wei(ph%iKap)%u + &
        length * ph%en * (ph%ray%x * ph%ray%vx + ph%ray%y * ph%ray%vy) / t
      c%optical%dir_wei(ph%iKap)%v = c%optical%dir_wei(ph%iKap)%v + &
        length * ph%en * (ph%ray%x * ph%ray%vy - ph%ray%y * ph%ray%vx) / t
      c%optical%dir_wei(ph%iKap)%w = c%optical%dir_wei(ph%iKap)%w + &
        length * ph%en * ph%ray%vz
    end if
    !
    if (encountered) then
      ! Something must be happening within this cell.
      call find_encounter_type(ph%iKap, c%optical, itype)
      !
      select case (itype)
        case (1)
          ! Gas absorption of X-ray
          destructed = .true.
          return
          !write(*,'(A)') 'In walk_scatter_absorb_reemit:'
          !write(*, '(A/)') 'H absorption: not possible!'
          !stop
        case (2) ! H scattering of Lya or X-ray
          ! Project the lambda into the local rest frame
          ph%lam = get_doppler_lam(a_star%mass, ph%lam, ph%ray)
          !
          c%par%sc_count_HI = c%par%sc_count_HI + 1
          call get_reemit_dir_uniform(ph%ray)
          ! write(*,'(A/)') 'Scattered by H atom!'
        case (3) ! Water absorption; no reemission
          destructed = .true.
          c%par%ab_count_water = c%par%ab_count_water + 1
          c%par%ab_en_water = c%par%ab_en_water + ph%en
          ! write(*,'(A/)') 'Absorbed by water!'
          return
        case (4)
          write(*,'(A)') 'In walk_scatter_absorb_reemit:'
          write(*, '(A/)') 'Should not happen: water scattering!'
          stop
        case default
          idust = (itype+1)/2 - ncl_nondust
          if (mod(itype, 2) .eq. 1) then ! Dust absorption
            if ((ph%lam .ge. lam_range_Xray(1)) .and. &
                (ph%lam .le. lam_range_Xray(2))) then
              ! X-ray absorption
              destructed = .true.
              return
            end if
            c%par%ab_count_dust = c%par%ab_count_dust + 1
            c%par%en_gains_abso(idust) = c%par%en_gains_abso(idust) + ph%en
            call dust_reemit(ph, c, idust)
            if (ph%lam .le. 0D0) then
              ! Absorbed by a dust grain that still owes energy to the gas.
              destructed = .true.
              return
            end if
          else ! Dust scattering, including X-ray
            ! Project the lambda into the local rest frame
            ph%lam = get_doppler_lam(a_star%mass, ph%lam, ph%ray)
            call get_reemit_dir_HenyeyGreenstein(ph%ray, &
                 dusts%list(idust)%g(ph%iKap))
          end if
          !if (itype .gt. c%optical%ntype) then
          !  write(*,'(A)') 'In walk_scatter_absorb_reemit:'
          !  write(*,'(A/)') 'Should not have this case.'
          !  stop
          !end if
      end select
      !
      ! Project the photon lambda back into the global rest frame.
      ph%lam = project_doppler_lam(a_star%mass, ph%lam, ph%ray)
      !
      ! Encounter has occurred, so we need to generate a new tau
      call random_number(rnd)
      tau = -log(rnd)
      !
    else
      if (mc_conf%ph_init_symmetric) then
        call locate_photon_cell_alt_mirror(r, z, c, dirtype, cnext, found)
      else
        call locate_photon_cell_alt(r, z, c, dirtype, cnext, found)
      end if
      if (.not. found) then! Not entering a neighboring cell
        ! May be entering a non-neighboring cell?
        if (mc_conf%ph_init_symmetric) then
          call enter_the_domain_mirror(ph%ray, cstart, cnext, found)
        else
          call enter_the_domain(ph%ray, cstart, cnext, found)
        end if
        if (.not. found) then ! Escape
          escaped = .true.
          return
        end if
        !!!! Todo
        !if (ph%ray%z .lt. 0D0) then
        !  ph%ray%z  = -ph%ray%z
        !  ph%ray%vz = -ph%ray%vz
        !end if
      end if
      c => cnext
      ! Each time entering a cell, the cell will gain a crossing count.
      c%optical%cr_count = c%optical%cr_count + 1
    end if
  end do
end subroutine walk_scatter_absorb_reemit



subroutine dust_reemit(ph, c, idust)
  type(type_photon_packet), intent(inout) :: ph
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in) :: idust
  double precision Tdust_old
  integer idx
  !
  if (c%par%n_dusts(idust) .le. 1D-100) then
    ph%lam = -1D0
    return
  end if
  Tdust_old = c%par%Tdusts(idust)
  !
  c%par%Tdusts(idust) = get_Tdust_from_LUT( &
    (c%par%en_gains(idust) + c%par%en_exchange(idust)) &
    / (4D0*phy_Pi * c%par%mdusts_cell(idust)), &
    luts%list(idust), idx)
  !if (c%par%Tdusts(idust) .le. 0D0) then
  !  write(*, '(/A)') 'In dust_reemit:'
  !  write(*, '(A, I4)') 'Tdust(i)<=0: ', idust
  !  write(*, *) c%par%en_gains
  !  write(*, *) c%par%en_exchange
  !  write(*, *)
  !end if
  if ((Tdust_old .ge. c%par%Tdusts(idust)) .or. &
      (Tdust_old .le. 0D0) .or. &
      (c%par%Tdusts(idust) .le. 0D0) .or. &
      isnan(c%par%Tdusts(idust))) then
    ph%lam = -1D0
    return
  end if
  !
  c%par%kphs(idust) = &
    (c%par%en_prevs(idust) + c%par%en_exchange(idust)) / ph%en
  !
  if (c%par%kphs(idust) .lt. 0D0) then
    ph%lam = -1D0
    return
  end if
  !
  ph%lam = get_reemit_lam(Tdust_old, &
                          c%par%Tdusts(idust), &
                          c%par%kphs(idust), &
                          luts%list(idust), &
                          dusts%list(idust), &
                          idx)
  !
  call get_reemit_dir_uniform(ph%ray)
end subroutine dust_reemit



function get_Tdust_from_LUT(val, lut, idx)
  double precision get_Tdust_from_LUT
  double precision, intent(in) :: val
  type(type_LUT_Tdust), intent(in) :: lut
  integer, intent(out) :: idx
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  logical found
  !
  found = .false.
  if (val .le. lut%vals(0)) then
    idx = 0
    get_Tdust_from_LUT = 0D0
    return
  else if (val .le. lut%vals(1)) then
    idx = 0
    get_Tdust_from_LUT = lut%Tds(idx) + &
      (val - lut%vals(idx)) * &
      (lut%Tds(idx+1) - lut%Tds(idx)) / &
      (lut%vals(idx+1) - lut%vals(idx))
    return
  else if (val .ge. lut%vals(lut%n)) then
    ! Using very crude extrapolation
    idx = lut%n
    get_Tdust_from_LUT = lut%Tds(lut%n) !* sqrt(sqrt(val / lut%vals(lut%n)))
    if (get_Tdust_from_LUT .ge. 1D4) then
      write(*, '(/A)') 'In get_Tdust_from_LUT:'
      write(*, '(A)') 'Abnormal Tdust!'
      write(*, '(ES16.6)') get_Tdust_from_LUT
      write(*, '(3ES16.6, I10/)') lut%Tds(lut%n), val, lut%vals(lut%n), idx
    end if
    return
  else if (isnan(val)) then
    write(*,'(/A)') 'In get_Tdust_from_LUT:'
    write(*,'(A/)') 'val is NaN!'
    !get_Tdust_from_LUT = phy_NaN
    stop
    !idx = -1
    !return
  else
    imin = 1
    imax = lut%n
    do i=1, lut%n
      if (imin .ge. imax-ITH) then
        do j=imin, imax-1
          if ((lut%vals(j) .le. val) .and. (lut%vals(j+1) .ge. val)) then
            idx = j
            found = .true.
            exit
          end if
        end do
        if (found) then
          get_Tdust_from_LUT = lut%Tds(idx) + &
            (val - lut%vals(idx)) * &
            (lut%Tds(idx+1) - lut%Tds(idx)) / &
            (lut%vals(idx+1) - lut%vals(idx))
          return
        else
          write(*,'(/A)') 'In get_Tdust_from_LUT:'
          write(*,'(A)') 'Cannot found idx:'
          write(*,'(2I5)') imin, imax
          write(*,'(3ES12.4/)') val, lut%vals(imin), lut%vals(imax)
          stop
        end if
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
  ilam = get_a_sample(p4lam, r)
  if (ilam .eq. 0) then
    write(*,'(A)') 'In get_reemit_lam:'
    write(*,'(A)') 'ilam = 0:'
    write(*,'(3ES16.6)') T0, T1, kph
    write(*,'(2I16)') idx0, idx1
    write(*,'(4ES12.4,/)') a, b, c0, c1
    write(*,'(10ES16.4)') p4lam%pvals(0:9)
    stop
  end if
  get_reemit_lam = dust%lam(ilam)
end function get_reemit_lam



pure function get_doppler_lam(M, lam0, ray) result(lam)
  double precision lam
  double precision, intent(in) :: lam0, M
  type(type_ray), intent(in) :: ray
  !
  double precision v, r, vd, tmp
  !
  tmp = ray%x*ray%x + ray%y*ray%y
  r = sqrt(tmp + ray%z*ray%z)
  v = sqrt((phy_GravitationConst_CGS * phy_Msun_CGS / phy_AU2cm) &
           * M / r)
  !
  r = v / sqrt(tmp)
  vd = (-ray%y * ray%vx + ray%x * ray%vy) * r
  !
  lam = lam0 * (1D0 + vd/phy_SpeedOfLight_CGS)
end function get_doppler_lam



pure function get_doppler_nu(M, nu0, ray) result(nu)
  ! The disk is rotating around the z axis in counterclockwise direction.
  ! The output nu is the frequency as seen by the fluid element.
  double precision nu
  double precision, intent(in) :: nu0
  type(type_ray), intent(in) :: ray
  double precision, intent(in) :: M
  !
  double precision v, r, vd
  !
  r = sqrt(ray%x*ray%x + ray%y*ray%y + ray%z*ray%z)
  v = sqrt(phy_GravitationConst_CGS * &
           M * phy_Msun_CGS / (r * phy_AU2cm))
  r = v / sqrt(ray%x*ray%x + ray%y*ray%y)
  !
  vd = (-ray%y * ray%vx + ray%x * ray%vy) * r
  !
  nu = nu0 * (1D0 - vd/phy_SpeedOfLight_CGS)
end function get_doppler_nu



pure function project_doppler_lam(M, lam0, ray) result(lam)
  double precision lam
  double precision, intent(in) :: lam0
  type(type_ray), intent(in) :: ray
  double precision, intent(in) :: M
  !
  double precision v, r, vd, tmp
  !
  tmp = ray%x*ray%x + ray%y*ray%y
  r = sqrt(tmp + ray%z*ray%z)
  v = sqrt((phy_GravitationConst_CGS * phy_Msun_CGS / phy_AU2cm) &
           * M / r)
  r = v / sqrt(tmp)
  !
  vd = (-ray%y * ray%vx + ray%x * ray%vy) * r
  !
  lam = lam0 * (1D0 - vd/phy_SpeedOfLight_CGS)
end function project_doppler_lam



subroutine find_encounter_type(iKap, coll, itype)
  integer, intent(in) :: iKap
  type(type_local_encounter_collection), intent(in) :: coll
  integer, intent(out) :: itype
  double precision r
  integer i
  !
  call random_number(r)
  !
  do i=1, coll%ntype
    if (r .lt. coll%acc(iKap, i)) then
      itype = i
      return
    end if
  end do
  !
  itype = coll%ntype
  !
end subroutine find_encounter_type



pure function get_idx_for_kappa(lam, dust)
  ! Actually only the dust%lam is needed.
  integer get_idx_for_kappa
  double precision, intent(in) :: lam
  type(type_optical_property), intent(in) :: dust
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  if ((lam .lt. dust%lam(1)) .or. (lam .ge. dust%lam(dust%n))) then
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
  ! Should never reach here.
  get_idx_for_kappa = 0
end function get_idx_for_kappa



subroutine make_local_optics(c, glo)
  ! This subroutine does not depend on the order of HI, water, and dust.
  ! acc: accumulative (along the axis of the index of optical materials)
  ! scattering + absorption opacity at each wavelength.
  ! summed_ab: total absorption opacity
  ! summed_sc: total scattering opacity
  ! summed: summed_ab + summed_sc
  type(type_cell), pointer, intent(inout) :: c
  type(type_global_material_collection), intent(in) :: glo
  type(type_local_encounter_collection), pointer :: loc
  integer i, i0
  !
  loc => c%optical
  !
  loc%acc = 0D0
  !
  ! X-ray cross sections for gas attributed to H.
  loc%acc(:, 1) = glo%list(1)%ab * loc%X(1) + glo%Xray_gas_abs * c%par%n_gas
  loc%acc(:, 2) = loc%acc(:, 1) + &
                  glo%list(1)%sc * loc%X(1) + glo%Xray_gas_sca * c%par%n_gas
  loc%summed_ab = glo%list(1)%ab * loc%X(1) + glo%Xray_gas_abs * c%par%n_gas
  loc%summed_sc = glo%list(1)%sc * loc%X(1) + glo%Xray_gas_sca * c%par%n_gas
  !
  do i=4, glo%ntype*2, 2
    i0 = i/2
    loc%acc(:, i-1) = loc%acc(:, i-2) + glo%list(i0)%ab * loc%X(i0)
    loc%acc(:, i)   = loc%acc(:, i-1) + glo%list(i0)%sc * loc%X(i0)
    loc%summed_ab = loc%summed_ab + glo%list(i0)%ab * loc%X(i0)
    loc%summed_sc = loc%summed_sc + glo%list(i0)%sc * loc%X(i0)
  end do
  !
  ! X-ray cross sections for dust attributed to the last type of material
  loc%acc(:, glo%ntype*2-1) = loc%acc(:, glo%ntype*2-1) + &
                  glo%Xray_dus_abs * c%par%n_gas * c%par%dust_depletion
  loc%acc(:, glo%ntype*2)   = loc%acc(:, glo%ntype*2)   + &
                  glo%Xray_dus_abs * c%par%n_gas * c%par%dust_depletion + &
                  glo%Xray_dus_sca * c%par%n_gas * c%par%dust_depletion
  loc%summed_ab = loc%summed_ab + &
                  glo%Xray_dus_abs * c%par%n_gas * c%par%dust_depletion
  loc%summed_sc = loc%summed_sc + &
                  glo%Xray_dus_sca * c%par%n_gas * c%par%dust_depletion
  !
  loc%summed = loc%acc(:, glo%ntype*2)
  !
  ! Normalize: this is because acc is used for finding out the encounter type
  do i=1, glo%ntype*2
    loc%acc(:, i) = loc%acc(:, i) / (loc%summed + 1D-100)
  end do
end subroutine make_local_optics



function get_stellar_luminosity(star, lam1, lam2) result(lumi)
  double precision lumi
  type(type_stellar_params), intent(in) :: star
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
  type(type_stellar_params), intent(inout) :: star
  integer i, fU
  integer nrows, ios
  character(len=128) str
  !
  nrows = GetFileLen_comment_blank(fname, '!')
  star%n = nrows
  if (allocated(star%lam)) then
    deallocate(star%lam, star%vals)
  end if
  allocate(star%lam(star%n), star%vals(star%n))
  !
  call openFileSequentialRead(fU, fname, 99, getu=1)
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
  type(type_stellar_params), intent(inout) :: star
  integer i
  double precision coeff
  !
  coeff = 4D0*phy_Pi**2 * (star%radius*phy_Rsun_CGS)**2
  star%n = nlam
  if (allocated(star%lam)) then
    deallocate(star%lam, star%vals)
  end if
  allocate(star%lam(star%n), star%vals(star%n))
  call logspace(star%lam, log10(lam0), log10(lam1), star%n)
  do i=1, star%n
    star%vals(i) = planck_B_lambda(star%T, &
      star%lam(i)*phy_Angstrom2cm) * coeff * phy_Angstrom2cm
  end do
  !
end subroutine make_stellar_spectrum




subroutine make_stellar_spectrum_Xray(nlam, star)
  ! Energy in keV; E0 < E1
  integer, intent(in) :: nlam
  type(type_stellar_params), intent(inout) :: star
  integer i
  double precision E, E0, E1, lam_min, lam_max
  !
  star%n = nlam
  if (allocated(star%lam)) then
    deallocate(star%lam, star%vals)
  end if
  allocate(star%lam(star%n), star%vals(star%n))
  !
  E0 = star%E0_Xray
  E1 = star%E1_Xray
  lam_min = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
            / (E1*1D3*phy_eV2erg) * 1D8 ! angstrom
  lam_max = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
            / (E0*1D3*phy_eV2erg) * 1D8 ! angstrom
  !
  call logspace(star%lam, log10(lam_min), log10(lam_max), star%n)
  !
  do i=1, star%n
    E = phy_hPlanck_CGS * phy_SpeedOfLight_CGS &
        / (star%lam(i) / 1D8) / (1D3*phy_eV2erg)
    star%vals(i) = exp(-(E*1D3*phy_eV2erg) / &
                        (phy_kBoltzmann_CGS*star%T_Xray))&
                   / (star%lam(i)**2)
  end do
  !
  star%vals = star%vals * &
              (star%lumi_Xray / &
               get_stellar_luminosity(star, lam_min, lam_max))
  !
end subroutine make_stellar_spectrum_Xray




function get_surf_max_angle(r0, z0)
  double precision get_surf_max_angle
  double precision, intent(in), optional :: r0, z0
  integer i, i0
  double precision r, z, w
  double precision r1, z1
  if (present(r0) .and. present(z0)) then
    r1 = r0
    z1 = z0
  else
    r1 = 0D0
    z1 = 0D0
  end if
  get_surf_max_angle = 0D0
  do i=1, surf_cells%nlen
    i0 = surf_cells%idx(i)
    r = leaves%list(i0)%p%xmin - r1
    z = leaves%list(i0)%p%ymax - z1
    w = z / sqrt(r*r + z*z)
    if (w .gt. get_surf_max_angle) then
      get_surf_max_angle = w
    end if
  end do
end function get_surf_max_angle




function get_bott_min_angle(r0, z0)
  double precision get_bott_min_angle
  double precision, intent(in), optional :: r0, z0
  integer i, i0
  double precision r, z, w
  double precision r1, z1
  if (present(r0) .and. present(z0)) then
    r1 = r0
    z1 = z0
  else
    r1 = 0D0
    z1 = 0D0
  end if
  get_bott_min_angle = 1D99
  do i=1, bott_cells%nlen
    i0 = bott_cells%idx(i)
    r = leaves%list(i0)%p%par%rmin - r1
    z = leaves%list(i0)%p%par%zmin - z1
    w = z / sqrt(r*r + z*z)
    if (w .lt. get_bott_min_angle) then
      get_bott_min_angle = w
    end if
  end do
end function get_bott_min_angle




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
  call openFileSequentialRead(fU, fname, 99, getu=1)
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



!subroutine load_dust_data(fname, dust)
!  character(len=*), intent(in) :: fname
!  type(type_optical_property), intent(out) :: dust
!  integer i, fU
!  integer iformat, nrows, ios
!  character(len=128) str
!  call openFileSequentialRead(fU, fname, 299, getu=1)
!  !
!  i = 0
!  do
!    call read_a_nonempty_row(fU, str, '(A128)', ios)
!    if (ios .ne. 0) then
!      exit
!    end if
!    ios = index(str, 'iformat')
!    if (ios .gt. 0) then
!      read(str((ios + len('iformat') + 3):), '(I16)') iformat
!      i = i + 1
!    end if
!    ios = index(str, 'nrows')
!    if (ios .gt. 0) then
!      read(str((ios + len('nrows') + 3):), '(I16)') nrows
!      i = i + 1
!    end if
!    if (i .eq. 2) then
!      exit
!    end if
!  end do
!  if (i .ne. 2) then
!    write(*,'(A)')  'In load_dust_data:'
!    write(*,'(A)') 'Cannot get format information!'
!    write(*,'(A/)') 'Will assume old format.'
!    !
!    rewind(fU)
!    call read_a_nonempty_row(fU, str, '(A128)', ios)
!    read(str, '(I16)') iformat
!    call read_a_nonempty_row(fU, str, '(A128)', ios)
!    read(str, '(I16)') nrows
!  end if
!  !
!  dust%n = nrows
!  write(*, '(A, I6)') "Number of lam: ", dust%n
!  !
!  allocate(dust%lam(nrows), dust%ab(nrows), dust%sc(nrows), dust%g(nrows))
!  dust%lam = 0D0
!  dust%ab = 0D0
!  dust%sc = 0D0
!  dust%g = 0D0
!  !
!  i = 0
!  do
!    call read_a_nonempty_row(fU, str, '(A128)', ios)
!    if (ios .ne. 0) then
!      exit
!    end if
!    if ((str(1:1) .eq. '!') .or. (str(1:1) .eq. '#')) then
!      cycle
!    end if
!    i = i + 1
!    select case(iformat)
!      case(1)
!        read(str, '(2(F18.8, X))') dust%lam(i), dust%ab(i)
!      case(2)
!        read(str, '(3(F18.8, X))') dust%lam(i), dust%ab(i), dust%sc(i)
!      case(3)
!        read(str, '(4(F18.8, X))') dust%lam(i), dust%ab(i), dust%sc(i), &
!            dust%g(i)
!    end select
!    ! write(*,*) i, dust%lam(i), dust%ab(i), dust%sc(i), dust%g(i)
!  end do
!  close(fU)
!  !
!  dust%lam = dust%lam / phy_Angstrom2micron
!end subroutine load_dust_data



subroutine make_H_Lya(T, hi)
  ! Zheng 2002
  ! A factor of c (speed of light) seems to be missing in that paper.
  use voigt
  type(type_optical_property), intent(out) :: hi
  double precision, intent(in) :: T
  double precision, parameter :: f12 = 0.4162D0
  double precision ratio, coeff
  double precision nu, nu_, x, x_, dx, dnu_th, a
  double precision, parameter :: xmax=5D3
  integer i, n2, idx, idx_
  !
  dnu_th = phy_LyAlpha_nu0 * &
    sqrt(8D0*phy_kBoltzmann_SI*T/phy_Pi/phy_mProton_SI) &
    / phy_SpeedOfLight_SI
  a = phy_LyAlpha_dnul / (2D0 * dnu_th)
  coeff = phy_LyAlpha_f12 * sqrt(phy_Pi) * &
          phy_electronClassicalRadius_CGS * &
          phy_SpeedOfLight_CGS / dnu_th
  !
  hi%n = 200
  n2 = (hi%n)/2
  !
  if (allocated(hi%lam)) then
    deallocate(hi%lam, hi%ab, hi%sc, hi%g)
  end if
  allocate(hi%lam(hi%n), hi%ab(hi%n), hi%sc(hi%n), hi%g(hi%n))
  hi%ab = 0D0
  hi%g = 0D0
  !
  x  =  1D-3
  x_ = -1D-3
  dx =  1D-2
  ratio = get_ratio_of_interval_log(x, xmax, dx, n2)
  do i=1, n2
    nu   = x  * dnu_th + phy_LyAlpha_nu0
    nu_  = x_ * dnu_th + phy_LyAlpha_nu0
    idx  = n2 - i + 1
    idx_ = n2 + i
    hi%lam(idx)  = phy_SpeedOfLight_SI / nu  * 1D10
    hi%lam(idx_) = phy_SpeedOfLight_SI / nu_ * 1D10
    !
    hi%sc(idx)  = coeff * max(0D0, voigt_scalar(abs(x),  a))
    hi%sc(idx_) = coeff * max(0D0, voigt_scalar(abs(x_), a))
    !
    !write(*, '(I4, 7ES14.4)') &
    !    i, a, x, hi%lam(idx), hi%sc(idx), x_, hi%lam(idx_), hi%sc(idx_)
    dx = dx * ratio
    x  = x  + dx
    x_ = x_ - dx
  end do
  !
end subroutine make_H_Lya



function get_Tdust_Stefan_Boltzmann(flux)
  double precision get_Tdust_Stefan_Boltzmann
  double precision, intent(in) :: flux
  get_Tdust_Stefan_Boltzmann = sqrt(sqrt(flux / phy_StefanBoltzmann_CGS))
end function get_Tdust_Stefan_Boltzmann



subroutine make_LUT_Tdust(dust, lut, nlen_lut, Tmin, Tmax)
  type(type_optical_property), intent(in) :: dust
  type(type_LUT_Tdust), intent(out) :: lut
  integer, intent(in) :: nlen_lut
  double precision, intent(in) :: Tmin, Tmax
  integer i, j
  double precision r, dT0, dT, tmp
  lut%n = nlen_lut ! Number of temperatures
  lut%m = dust%n - 1 ! Number of lambdas
  allocate(lut%Tds(0:lut%n), &
           lut%vals(0:lut%n), &
           lut%table(0:lut%m, 0:lut%n))
  lut%Tds(0) = 0D0
  lut%vals(0) = 0D0
  lut%table(0, :) = 0D0
  lut%table(:, 0) = 0D0
  dT0 = 1D-1 * Tmin
  dT = dT0
  r = get_ratio_of_interval_log(Tmin, Tmax, dT0, lut%n)
  do i=1, lut%n
    if (i .eq. 1) then
      lut%Tds(i) = Tmin
    else if (i .eq. lut%n) then
      lut%Tds(i) = max(lut%Tds(i-1)+dT*0.1D0, Tmax)
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




function get_a_sample(distri, r)
  integer get_a_sample
  type(type_distribution_table), intent(in) :: distri
  double precision, intent(in) :: r
  integer i, j, imin, imax, imid
  integer, parameter :: ITH = 5
  !
  get_a_sample = 0
  !
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



pure function tau2frac(tau)
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
        t = 2D0 * sqrt(1D0 - s)
        ray%vx = t * x(1)
        ray%vy = t * x(2)
        return
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
  phi = phy_2Pi * p(2)
  dir_rel%u = sintheta * cos(phi)
  dir_rel%v = sintheta * sin(phi)
  dir_rel%w = costheta
  dir0 = combine_dir(dir0, dir_rel)
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


type(type_direction_cartesian) function combine_dir(dir0, dir_rel)
  type(type_direction_cartesian) dir0, dir_rel
  type(type_sphere_coor_quat) d
  d = cartesian_to_spherical_quat(dir0)
  combine_dir = rot_around_Y(dir_rel, d%costheta, d%sintheta)
  combine_dir = rot_around_Z(combine_dir, &
    d%cosphi, d%sinphi)
end function combine_dir


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


subroutine set_up_collector(collector, minlam, maxlam, dmu, nmu, nr, nphi, ang_mins, ang_maxs)
  type(type_photon_collector), intent(inout) :: collector
  double precision, intent(in), optional :: maxlam, minlam, dmu
  integer, intent(in), optional :: nmu, nr, nphi
  double precision, dimension(:), intent(in), optional :: ang_mins, ang_maxs
  double precision dmu_size, delmu, delphi
  integer i
  !
  if (present(dmu)) then
    dmu_size = dmu
  else
    dmu_size = 0.1D0
  end if
  !
  if (present(nmu)) then
    collector%nmu = nmu
  else
    collector%nmu = 5
  end if
  !
  allocate(collector%mu_min(collector%nmu), collector%mu_max(collector%nmu))
  if (present(ang_mins) .and. present(ang_maxs)) then
    do i=1, collector%nmu
      collector%mu_min(i) = cos(ang_maxs(i) * phy_Pi / 180D0)
      collector%mu_max(i) = cos(ang_mins(i) * phy_Pi / 180D0)
    end do
  else
    delmu = (1D0-dmu_size) / dble(collector%nmu-1)
    do i=1, collector%nmu
      collector%mu_min(i) = dble(i-1) * delmu
      collector%mu_max(i) = collector%mu_min(i) + dmu_size
    end do
  end if
  !
  if (present(maxlam) .and. present(minlam)) then
    collector%iKap0 = max(1, get_idx_for_kappa(minlam, dust_0))
    if (maxlam .ge. dust_0%lam(dust_0%n-1)) then
      collector%iKap1 = dust_0%n - 1
    else
      collector%iKap1 = max(1, get_idx_for_kappa(maxlam, dust_0))
    end if
  else
    collector%iKap0 = 1
    collector%iKap1 = dust_0%n - 1
  end if
  !
  collector%nlam = collector%iKap1 - collector%iKap0 + 1
  !
  if (present(nr)) then
    collector%nr = nr
  else
    collector%nr = 50
  end if
  !
  if (present(nphi)) then
    collector%nphi = nphi
  else
    collector%nphi = 50
  end if
  !
  allocate(collector%r_min(collector%nr), collector%r_max(collector%nr), &
           collector%phi_min(collector%nphi), collector%phi_max(collector%nphi))
  collector%r_min(1) = 0D0
  collector%r_max(1) = root%xmin * 0.3D0
  !
  call logspace(collector%r_min(2:collector%nr), &
      log10(collector%r_max(1)), log10(root%xmax), collector%nr-1)
  !
  do i=2, collector%nr-1
    collector%r_max(i) = collector%r_min(i+1)
  end do
  collector%r_max(collector%nr) = &
    collector%r_min(collector%nr)*2D0 - collector%r_min(collector%nr-1)
  !
  delphi = phy_2Pi / dble(collector%nphi)
  collector%phi_min(1) = -phy_Pi
  collector%phi_max(1) = collector%phi_min(1) + delphi
  do i=2, collector%nphi
    collector%phi_min(i) = collector%phi_min(i-1) + delphi
    collector%phi_max(i) = collector%phi_min(i)   + delphi
  end do
  !
  allocate(collector%energy(collector%nmu, collector%nphi, collector%nr, collector%nlam), &
           collector%counts(collector%nmu, collector%nphi, collector%nr, collector%nlam))
  !
  collector%energy = 0D0
  collector%counts = 0
end subroutine set_up_collector



pure subroutine collect_photon_do(collector, ph)
  type(type_photon_collector), intent(inout) :: collector
  type(type_photon_packet), intent(in) :: ph
  type(type_position_cartesian) r_ortho
  type(type_direction_cartesian) uz0, ux, uy, uz
  double precision r_o_x, r_o_y, mu, phi, r
  integer i, imu, iphi, ir
  !
  if ((ph%iKap .lt. collector%iKap0) .or. &
      (ph%iKap .gt. collector%iKap1)) then
    return
  end if
  !
  ! First filter the elevation angle
  mu = abs(ph%ray%vz)
  imu = 0
  do i=1, collector%nmu
    if ((collector%mu_min(i) .le. mu) .and. &
        (mu .le. collector%mu_max(i))) then
      imu = i
      exit
    end if
  end do
  if (imu .eq. 0) then
    return
  end if
  !
  ! Get the component of the position vector orthogonal to the photon ray
  ! direction.  This is the displacement vector relative to the image center.
  call get_r_ortho(ph%ray, r_ortho)
  !
  ! Set up a local image coordinate frame.
  ! The photon ray direction is set to be the z axis.
  uz%u = ph%ray%vx
  uz%v = ph%ray%vy
  uz%w = ph%ray%vz
  if (abs(uz%w) .ge. 0.99D0) then ! Degenerate with uz0
    ux%u = 1D0; ux%v = 0D0; ux%w = 0D0
    uy%u = 0D0; uy%v = 1D0; uy%w = 0D0
  else
    ! Use the global z axis to set up the local axis.
    uz0%u = 0D0
    uz0%v = 0D0
    uz0%w = 1D0
    call cross_product(uz0, uz, ux, normalize=.true.)
    call cross_product(uz, ux, uy)
  end if
  !
  r_o_x = r_ortho%x * ux%u + r_ortho%y * ux%v + r_ortho%z * ux%w
  r_o_y = r_ortho%x * uy%u + r_ortho%y * uy%v + r_ortho%z * uy%w
  r = sqrt(r_o_x**2 + r_o_y**2)
  phi = atan2(r_o_y, r_o_x)
  !
  ! Find the radial and angular index
  ir = 0
  do i=1, collector%nr
    if ((collector%r_min(i) .le. r) .and. &
        (r .le. collector%r_max(i))) then
      ir = i
      exit
    end if
  end do
  if (ir .eq. 0) then
    return
  end if
  !
  iphi = 0
  do i=1, collector%nphi
    if ((collector%phi_min(i) .le. phi) .and. &
        (phi .le. collector%phi_max(i))) then
      iphi = i
      exit
    end if
  end do
  if (iphi .eq. 0) then
    return
  end if
  !
  collector%energy(imu, iphi, ir, ph%iKap-collector%iKap0+1) = &
    collector%energy(imu, iphi, ir, ph%iKap-collector%iKap0+1) + ph%en
  collector%counts(imu, iphi, ir, ph%iKap-collector%iKap0+1) = &
    collector%counts(imu, iphi, ir, ph%iKap-collector%iKap0+1) + 1
  !
end subroutine collect_photon_do


subroutine save_collected_photons(fname, collector, imu)
  character(len=*), intent(in) :: fname
  type(type_photon_collector), intent(inout) :: collector
  integer, intent(in) :: imu
  !double precision r, dr, dphi, dArea
  double precision dOmega, flux, f1, f2
  integer i, j, k, fU, ct
  !
  call openFileSequentialWrite(fU, fname, 99, getu=1)
  !
  dOmega = 2D0*phy_2Pi * (collector%mu_max(imu) - collector%mu_min(imu))
  !
  do i=1, collector%nlam
    flux = 0D0
    ct = 0
    do j=1, collector%nr
      !dr = collector%r_max(j) - collector%r_min(j)
      !r = 0.5D0 * (collector%r_max(j) + collector%r_min(j))
      do k=1, collector%nphi
        !dphi = collector%phi_max(k) - collector%phi_min(k)
        !dArea = dr * r * dphi * phy_AU2cm**2
        ! The factor 2 is to account for the fact that the mirror symmetry was
        ! used to count the photons.
        !flux = flux + collector%energy(imu, j, k, i) / (dArea * dOmega)
        flux = flux + collector%energy(imu, j, k, i)
        ct = ct + collector%counts(imu, j, k, i)
      end do
    end do
    f1 = phy_SpeedOfLight_SI / (dust_0%lam(i-1+collector%iKap0) * 1D-10)
    f2 = phy_SpeedOfLight_SI / (dust_0%lam(i+collector%iKap0) * 1D-10)
    flux = flux / dOmega / (mc_conf%dist * phy_pc2cm)**2 / phy_jansky2CGS / (f1-f2)
    write(fU, '(ES16.6, ES16.6, I10)') &
        dust_0%lam(i-1+collector%iKap0), flux, ct
  end do
  close(fU)
end subroutine save_collected_photons


subroutine save_collected_photons_iter(iiter)
  integer, intent(in) :: iiter
  character(len=128) :: fname
  integer i
  do i=1, collector%nmu
    write(fname, '(2(A, F0.2), A, I0.4, A)') 'spec_ang_', &
      acos(collector%mu_max(i))*180D0/phy_Pi, '_', &
      acos(collector%mu_min(i))*180D0/phy_Pi, &
      '_iter_', iiter, '.dat'
    call save_collected_photons( &
      combine_dir_filename(mc_conf%mc_dir_out, fname), collector, i)
  end do
end subroutine save_collected_photons_iter



pure subroutine cross_product(A, B, C, normalize)
  type(type_direction_cartesian), intent(in) :: A, B
  type(type_direction_cartesian), intent(out) :: C
  logical, intent(in), optional :: normalize
  double precision norm
  C%u = A%v * B%w - A%w * B%v
  C%v = A%w * B%u - A%u * B%w
  C%w = A%u * B%v - A%v * B%u
  if (present(normalize)) then
    if (normalize) then
      norm = sqrt(C%u*C%u + C%v*C%v + C%w*C%w)
      C%u = C%u / norm
      C%v = C%v / norm
      C%w = C%w / norm
    end if
  end if
end subroutine cross_product


pure subroutine get_r_ortho(ray, r_ortho)
  type(type_ray), intent(in) :: ray
  type(type_position_cartesian), intent(out) :: r_ortho
  double precision rdotv
  rdotv = ray%x * ray%vx + ray%y * ray%vy + ray%z * ray%vz
  r_ortho%x = ray%x - rdotv * ray%vx
  r_ortho%y = ray%y - rdotv * ray%vy
  r_ortho%z = ray%z - rdotv * ray%vz
end subroutine get_r_ortho


end module montecarlo
