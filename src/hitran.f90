! Ref: 2012-Rothman-The HITRAN2012 molecular spectroscopic database
! 2014-04-13 Sun 17:49:43
! Fujun Du

module hitran

use trivials
use data_struct
use phy_const
use quick_sort

implicit none

integer, parameter, private :: n_hitran_mol = 19
character(len=12), dimension(2, n_hitran_mol), private :: &
  hitran_mol_fnames = reshape( &
  (/'H2O         ',   '01_hit12.par', &
    'CO2         ',   '02_hit12.par', &
    'O3          ',   '03_hit12.par', &
    'CO          ',   '05_hit12.par', &
    'CH4         ',   '06_hit12.par', &
    'O2          ',   '07_hit12.par', &
    'NO          ',   '08_hit08.par', &
    'NH3         ',   '11_hit12.par', &
    'OH          ',   '13_hit12.par', &
    'N2          ',   '22_hit12.par', &
    'H2O2        ',   '25_hit08.par', &
    'C2H2        ',   '26_hit12.par', &
    'H2S         ',   '31_hit12.par', &
    'HCOOH       ',   '32_hit08.par', &
    'O           ',   '34_hit08.par', &
    'CH3CN       ',   '41_hit08.par', &
    'HC3N        ',   '44_hit12.par', &
    'H2          ',   '45_hit12.par', &
    'CS          ',   '46_hit12.par'  &
  /), (/2, n_hitran_mol/))


contains


subroutine load_hitran_mol(dir_name, mol_name, mol_data, &
  lam_range, Elow_range, tau_min, N_estimate, orthopara)
  ! Purpose: load the HITRAN molecular transition data
  ! The data structure will be stored in mol_data.
  ! lam_range should be in micron.
  ! Elow_range should be in K.
  ! tau_min and N_estimate are used to limit the line intensity.
  ! orthopara should be lower cases of 'ortho', 'para', or 'all'.
  !
  character(len=*), intent(in) :: dir_name, mol_name
  type(type_molecule_energy_set), pointer, intent(out) :: mol_data
  double precision, intent(in), dimension(2), optional :: lam_range, Elow_range
  double precision, intent(in), optional :: tau_min, N_estimate
  character(len=*), intent(in), optional :: orthopara
  !
  character(len=256) filename
  integer i, fU, idxmol, flen
  double precision, dimension(2) :: lam_r, Elow_r
  double precision :: tau_m, N_est, tau
  character(len=8) :: op
  !
  ! Following Table 1 of Rothman 2012
  integer, dimension(:), allocatable :: imol, iiso
  double precision, dimension(:), allocatable :: wavnum, inten, Acoeff, &
            halfWidthAir, halfWidthSelf, Elow, TcoeffAir, pShiftAir
  character(len=15), dimension(:), allocatable :: qUpperGl, qLowerGl, &
            qUpperLoc, qLowerLoc
  integer, dimension(:), allocatable :: iErr
  character(len=12), dimension(:), allocatable :: iRef
  character, dimension(:), allocatable :: flag
  double precision, dimension(:), allocatable :: gWeiUpper, gWeiLower
  !
  integer n_keep, n_unique, iOrthoPara, i0, i1, i2
  integer, dimension(:), allocatable :: idx_keep, idx_unique, idx_reverse
  double precision, dimension(:), allocatable :: Eup, Eall, gWeiAll
  double precision lam_micron, Elow_K
  !
  idxmol = 0
  do i=1, n_hitran_mol
    if (hitran_mol_fnames(1, i) .eq. mol_name) then
      idxmol = i
      exit
    end if
  end do
  if (idxmol .le. 0) then
    write(*, '(/A)') 'In load_hitran_mol:'
    write(*, '(2A)') mol_name, ' is not implemented in the database.'
    write(*, '(A/)') 'Contact the author for a solution.'
    return
  end if
  !
  filename = combine_dir_filename(dir_name, hitran_mol_fnames(2,idxmol))
  flen = GetFileLen(filename)
  !
  write(*, '(/2A)') 'Loading file: ', trim(filename)
  write(*, '(A, I10/)') 'Number of rows: ', flen
  !
  allocate(imol(flen), iiso(flen), wavnum(flen), inten(flen), &
    Acoeff(flen), halfWidthAir(flen), halfWidthSelf(flen), &
    Elow(flen), TcoeffAir(flen), pShiftAir(flen), &
    qUpperGl(flen), qLowerGl(flen), qUpperLoc(flen), qLowerLoc(flen), &
    iErr(flen), iRef(flen), flag(flen), gWeiUpper(flen), gWeiLower(flen))
  !
  allocate(idx_keep(flen), Eup(flen))
  !
  if (present(lam_range)) then
    lam_r = lam_range
  else
    lam_r = (/0D0, 1D99/)
  end if
  if (present(Elow_range)) then
    Elow_r = Elow_range
  else
    Elow_r = (/0D0, 1D99/)
  end if
  if (present(orthopara)) then
    op = orthopara
  else
    op = 'all'
  end if
  if (present(tau_min)) then
    tau_m = tau_min
  else
    tau_m = 0D0
  end if
  if (present(N_estimate)) then
    N_est = N_estimate
  else
    N_est = 1D25
  end if
  !
  call openFileSequentialRead(fU, filename, 256, getu=1)
  if (associated(mol_data)) then
    mol_data => null()
  end if
  allocate(mol_data)
  mol_data%name_molecule = mol_name
  !
  n_keep = 0
  do i=1, flen
    call read_a_line_hitran2012(fU, &
      imol(i), iiso(i), wavnum(i), inten(i), Acoeff(i), halfWidthAir(i), &
      halfWidthSelf(i), Elow(i), &
      TcoeffAir(i), pShiftAir(i), qUpperGl(i), qLowerGl(i), qUpperLoc(i), &
      qLowerLoc(i), &
      iErr(i), iRef(i), flag(i), gWeiUpper(i), gWeiLower(i))
    Eup(i) = Elow(i) + wavnum(i)
    !
    Elow_K = phy_cm_1_2K * Elow(i)
    !
    lam_micron = 1D4 / wavnum(i)
    !
    ! Assume dv = 1 km s-1, then c/dv = 3e5
    tau = inten(i) * N_est / wavnum(i) * 3D5
    !
    if ((lam_r(1) .le. lam_micron) .and. (lam_micron .le. lam_r(2)) .and. &
        (Elow_r(1) .le. Elow_K) .and. (Elow_K .le. Elow_r(2)) .and. &
        (tau .ge. tau_m)) then
      select case(op)
        case ('all')
          n_keep = n_keep + 1
          idx_keep(n_keep) = i
        case ('ortho', 'para')
          iOrthoPara = get_ortho_para(qUpperGl(i), qUpperLoc(i))
          if (((op .eq. 'ortho') .and. (iOrthoPara .eq. 1)) .or. &
              ((op .eq. 'para')  .and. (iOrthoPara .eq. 0))) then
            n_keep = n_keep + 1
            idx_keep(n_keep) = i
          end if
        case default
          write(*, '(2A)') 'Unknown ortho/para option: ', op
          write(*, '(A)') 'Must be lower case.'
      end select
    end if
    !
  end do
  !
  close(fU)
  !
  write(*, '(A)') 'Making the energy level structure.'
  !
  allocate(Eall(n_keep*2), gWeiAll(n_keep*2), idx_unique(n_keep*2), &
    idx_reverse(n_keep*2))
  !
  do i=1, n_keep
    i0 = idx_keep(i)
    Eall(2*i-1) = Elow(i0)
    Eall(2*i)   = Eup(i0)
    gWeiAll(2*i-1) = gWeiLower(i0)
    gWeiAll(2*i)   = gWeiUpper(i0)
  end do
  !
  call unique_vector_idx(Eall, 2*n_keep, idx_unique, n_unique, &
         1D-8, 1D-20, idx_reverse)
  mol_data%n_level = n_unique
  !
  allocate(mol_data%level_list(n_unique), &
           mol_data%f_occupation(n_unique))
  !
  mol_data%level_list%energy = Eall(idx_unique(1:n_unique)) * phy_cm_1_2K
  mol_data%level_list%weight = gWeiAll(idx_unique(1:n_unique))
  !
  allocate(mol_data%rad_data)
  !
  mol_data%rad_data%n_transition = n_keep
  !
  allocate(mol_data%rad_data%list(n_keep))
  !
  do i=1, n_keep
    i0 = idx_keep(i)
    mol_data%rad_data%list(i)%Aul  = Acoeff(i0)
    mol_data%rad_data%list(i)%freq = wavnum(i0) * phy_SpeedOfLight_CGS
    mol_data%rad_data%list(i)%lambda = & ! In angstrom
      phy_SpeedOfLight_SI / mol_data%rad_data%list(i)%freq * 1D10
    mol_data%rad_data%list(i)%Bul = mol_data%rad_data%list(i)%Aul / &
      ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * &
       (mol_data%rad_data%list(i)%freq)**3)
    !
    mol_data%rad_data%list(i)%Eup  = Eup(i0) * phy_cm_1_2K
    mol_data%rad_data%list(i)%Elow = Elow(i0) * phy_cm_1_2K
    !
    ! A bit tricky to get the index in the unique energy level list.
    mol_data%rad_data%list(i)%iup = idx_reverse(2*i-1)
      !binary_search(mol_data%level_list%energy, n_unique, &
      !  mol_data%rad_data%list(i)%Eup, 1)
    mol_data%rad_data%list(i)%ilow = idx_reverse(2*i)
      !binary_search(mol_data%level_list%energy, n_unique, &
      !  mol_data%rad_data%list(i)%Elow, 1)
    !
    i2 = mol_data%rad_data%list(i)%iup
    i1  = mol_data%rad_data%list(i)%ilow
    mol_data%rad_data%list(i)%Blu = mol_data%rad_data%list(i)%Bul * &
        mol_data%level_list(i2)%weight / mol_data%level_list(i1)%weight
  end do
  !
  write(*, '(/2A)') 'Molecule: ', mol_data%name_molecule
  write(*, '(A, I10)') 'Number of energy levels: ', mol_data%n_level
  write(*, '(A, 2F12.2)') 'Min,max (K): ', &
    mol_data%level_list(1)%energy, mol_data%level_list(mol_data%n_level)%energy
  write(*, '(A, I10)') 'Number of transitions: ', mol_data%rad_data%n_transition
  write(*, '(A, 2ES12.2)') 'Min,max (Hz): ', &
    mol_data%rad_data%list(1)%freq, &
    mol_data%rad_data%list(mol_data%rad_data%n_transition)%freq
  !
  allocate(mol_data%colli_data)
  ! No collisional data in the HITRAN database.
  mol_data%colli_data%n_partner = 0
  !
end subroutine load_hitran_mol


function get_ortho_para(qnum_gl, qnum_loc)
  ! Following radlite (hitran_extract.pro)
  ! 0: para
  ! 1: ortho
  integer get_ortho_para
  character(len=*), intent(in) :: qnum_gl, qnum_loc
  integer v3, ka, kc
  !
  ! Table 3 and 4 of the HITRAN04paper.
  read(qnum_gl, '(13X, I2)') v3
  read(qnum_loc, '(3X, I3, I3)') ka, kc
  !
  if (mod(ka+kc+v3, 2) .eq. 1) then
    ! Ortho
    get_ortho_para = 1
  else
    ! Para
    get_ortho_para = 0
  end if
end function get_ortho_para


subroutine read_a_line_hitran2012(fU, &
  imol, iiso, wavnum, inten, Acoeff, halfWidthAir, halfWidthSelf, Elow, &
  TcoeffAir, pShiftAir, qUpperGl, qLowerGl, qUpperLoc, qLowerLoc, &
  iErr, iRef, flag, gWeiUpper, gWeiLower)
  integer, intent(in) :: fU
  integer, intent(out) :: imol, iiso
  double precision, intent(out) :: wavnum, inten, Acoeff, &
            halfWidthAir, halfWidthSelf, Elow, TcoeffAir, pShiftAir
  character(len=15), intent(out) :: qUpperGl, qLowerGl, qUpperLoc, qLowerLoc
  integer, intent(out) :: iErr
  character(len=12), intent(out) :: iRef
  character, intent(out) :: flag
  double precision, intent(out) :: gWeiUpper, gWeiLower
  !
  read(fU, '(I2, I1, F12.0, F10.0, F10.0, F5.0, F5.0, F10.0, F4.0, F8.0, &
           & 4A15, I6, A12, A1, F7.0, F7.0)') &
    imol, iiso, wavnum, inten, Acoeff, halfWidthAir, halfWidthSelf, Elow, &
    TcoeffAir, pShiftAir, qUpperGl, qLowerGl, qUpperLoc, qLowerLoc, &
    iErr, iRef, flag, gWeiUpper, gWeiLower
  !
end subroutine read_a_line_hitran2012



subroutine test_load_ratran_adhoc
  type(type_molecule_energy_set), pointer :: molecule
  call load_hitran_mol(&
  '/Users/fdu/not_synced/work/from_others/hitran/HITRAN2012/HITRAN2012/By-Molecule/Uncompressed-files/', &
  'CO', molecule, &
  lam_range=(/1D0, 1D3/), &
  Elow_range=(/0.0D0, 5D3/), &
  tau_min=1D-6, N_estimate=1D22, orthopara='ortho')
  write(*,*) molecule%level_list(1:10)%energy
  write(*,*) molecule%level_list(1:10)%weight
  write(*,*) molecule%rad_data%list(1:10)%Aul
  write(*,*) molecule%rad_data%list(1:10)%Eup
  write(*,*) molecule%rad_data%list(1:10)%Elow
  stop
end subroutine test_load_ratran_adhoc

end module hitran