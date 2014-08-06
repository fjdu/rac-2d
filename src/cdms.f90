module cdms

use trivials
use data_struct
use phy_const
use quick_sort

implicit none

integer, private, parameter :: nT = 11
double precision, dimension(nT), private, parameter :: T_s = &
    (/1D3, 5D2, 3D2, 2.25D2, 1.5D2, 7.5D1, 3.75D1, 1.875D1, 9.375D0, 5D0, 2.725D0/)
double precision, dimension(nT), private :: lg10Q

contains

subroutine load_cdms_mol(dir_name, fname, fname_part, mol_data)
  character(len=*), intent(in) :: dir_name, fname, fname_part
  type(type_molecule_energy_set), pointer, intent(inout) :: mol_data
  !
  character(len=256) filename
  integer i, i1, i2, flen, fU
  double precision, dimension(:), allocatable :: freq, uncer, intens, Elow, Eup
  integer, dimension(:), allocatable :: dof, gup, glow, tag, cquan
  integer, dimension(:,:), allocatable :: quanup, quanlow
  double precision, dimension(:), allocatable :: Eall, gWeiAll
  integer, dimension(:), allocatable :: idx_unique, idx_reverse
  integer n_keep, n_unique
  double precision lowest_freq, atol
  !
  filename = combine_dir_filename(dir_name, fname)
  flen = GetFileLen(filename)
  !
  write(*, '(/2A)') 'Loading file: ', trim(filename)
  write(*, '(A, I10/)') 'Number of rows: ', flen
  !
  allocate(freq(flen), uncer(flen), intens(flen), dof(flen), &
           Elow(flen), Eup(flen), gup(flen), glow(flen), tag(flen), cquan(flen), &
           quanup(6, flen), quanlow(6, flen))
  !
  call openFileSequentialRead(fU, filename, 256, getu=1)
  !
  lowest_freq = 1D99
  do i=1, flen
    call read_a_line_cdms(fU, freq(i), uncer(i), intens(i), dof(i), &
        Elow(i), gup(i), tag(i), cquan(i), quanup(:,i), quanlow(:,i))
    !write(*,*) freq(i), uncer(i), intens(i), dof(i), &
    !    Elow(i), gup(i), tag(i), cquan(i), quanup(:,i), quanlow(:,i)
    freq(i) = freq(i) * 1D6
    Eup(i) = Elow(i) + freq(i) / phy_SpeedOfLight_CGS
    if ((2*quanup(1,i)+1) .ne. gup(i)) then
      write(*, '(A)') 'I do not know how to calcualte the statistical weight!'
      write(*, '(I6, 13I4)') i, gup(i), quanup(:,i), quanlow(:,i)
      stop
    end if
    glow(i) = 2*quanlow(1,i) + 1
    lowest_freq = min(lowest_freq, freq(i))
    !write(*,'(I4, 2ES20.12, 2I4)') i, Eup(i), Elow(i), gup(i), glow(i)
  end do
  !
  close(fU)
  !
  write(*, '(A)') 'Loading the cdms partition function.'
  call load_cdms_partition(dir_name, fname_part, abs(tag(1)))
  !
  write(*, '(A)') 'Making the energy level structure.'
  !
  n_keep = flen
  !
  allocate(Eall(n_keep*2), gWeiAll(n_keep*2), idx_unique(n_keep*2), &
    idx_reverse(n_keep*2))
  !
  do i=1, flen
    Eall(2*i-1) = Elow(i)
    Eall(2*i  ) = Eup(i)
    gWeiAll(2*i-1) = dble(glow(i))
    gWeiAll(2*i  ) = dble(gup(i))
  end do
  !
  atol = 1D-4 * lowest_freq / phy_SpeedOfLight_CGS
  !
  call unique_vector_idx(Eall, 2*n_keep, idx_unique, n_unique, &
         1D-10, atol, idx_reverse)
  !
  if (.not. associated(mol_data)) then
    allocate(mol_data)
  end if
  !
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
    mol_data%rad_data%list(i)%freq = freq(i)
    mol_data%rad_data%list(i)%lambda = & ! In angstrom
      phy_SpeedOfLight_SI / mol_data%rad_data%list(i)%freq * 1D10
    !
    mol_data%rad_data%list(i)%Eup  = Eup(i) * phy_cm_1_2K
    mol_data%rad_data%list(i)%Elow = Elow(i) * phy_cm_1_2K
    !
    write(mol_data%rad_data%list(i)%qnum, '(I4, ": ", 6I4, " -> ", 6I4)') &
        cquan(i), quanup(:, i), quanlow(:, i)
    !
    mol_data%rad_data%list(i)%ilow = idx_reverse(2*i-1)
    mol_data%rad_data%list(i)%iup = idx_reverse(2*i)
    !
    mol_data%rad_data%list(i)%Aul = &
      cdms_intensity2Aul(exp(log(10D0)*intens(i)), freq(i), &
      mol_data%rad_data%list(i)%Elow , mol_data%rad_data%list(i)%Eup, &
      dble(gup(i)), 3D2)
    !
    mol_data%rad_data%list(i)%Bul = mol_data%rad_data%list(i)%Aul / &
      ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * &
       (mol_data%rad_data%list(i)%freq)**3)
    !
    i2 = mol_data%rad_data%list(i)%iup
    i1  = mol_data%rad_data%list(i)%ilow
    mol_data%rad_data%list(i)%Blu = mol_data%rad_data%list(i)%Bul * &
      (mol_data%level_list(i2)%weight / mol_data%level_list(i1)%weight)
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
  ! No collisional data in the cdms database.
  mol_data%colli_data%n_partner = 0
  !
end subroutine load_cdms_mol


subroutine read_a_line_cdms(fU, freq, uncer, intens, dof, Elow, gup, tag, cquan, quanup, quanlow)
  integer, intent(in) :: fU
  double precision, intent(out) :: freq, uncer, intens, Elow
  integer, intent(out) :: dof, gup, tag, cquan
  integer, dimension(6), intent(out) :: quanup, quanlow
  !
  integer ios
  !
  read(fU, '(F13.0, F8.0, F8.0, I2, F10.0, I3, I7, I4, 6I2, 6I2)', IOSTAT=ios) &
    freq,  uncer, intens, dof, Elow, gup, tag, cquan, quanup, quanlow
end subroutine read_a_line_cdms


subroutine load_cdms_partition(dir_name, fname, moltag)
  character(len=*), intent(in) :: dir_name, fname
  integer, intent(in) :: moltag
  !
  character(len=256) filename
  character(len=256) srow
  integer i, i1, i2, flen, fU
  logical found
  !
  filename = combine_dir_filename(dir_name, fname)
  flen = GetFileLen(filename)
  !
  write(*, '(/2A)') 'Loading file: ', trim(filename)
  write(*, '(A, I10/)') 'Number of rows: ', flen
  !
  call openFileSequentialRead(fU, filename, 256, getu=1)
  found = .false.
  read(fU, '(A)') srow ! Skip the first row
  do i=2, flen
    read(fU, '(A)') srow
    read(srow, '(I7)') i1
    if (i1 .eq. moltag) then
      found = .true.
      write(*, '(A, I6)') 'Found in row ', i1
      exit
    end if
  end do
  close(fU)
  !
  if (.not. found) then
    write(*, '(A)') 'In loading cdms:'
    write(*, '(A)') 'Cannot find an entry in the partition function file!'
    write(*, '(A, I16)') 'moltag = ', moltag
    stop
  end if
  !
  read(srow, '(I7, 24X, I7, 11F13.0)') i1, i2, lg10Q
end subroutine load_cdms_partition



function cdms_calc_partition(T) result(Q)
  double precision Q
  double precision, intent(in) :: T
  integer i
  if (T .lt. T_s(1)) then
    Q = lg10Q(1)
  else if (T .gt. T_s(nT)) then
    Q = lg10Q(nT)
  else
    do i=1, nT-1
      if ((T_s(i) .le. T) .and. (T .le. T_s(i+1))) then
        Q = (lg10Q(i+1) - lg10Q(i)) / (T_s(i+1) - T_s(i)) &
            * (T - T_s(i)) + lg10Q(i)
        exit
      end if
    end do
  end if
  Q = exp(log(10D0) * Q)
end function cdms_calc_partition



function cdms_calc_partition_my(mol_data, T) result(Q)
  double precision Q
  type(type_molecule_energy_set), intent(in) :: mol_data
  double precision, intent(in) :: T
  integer i
  Q = 0D0
  do i=1, mol_data%n_level
    Q = Q + mol_data%level_list(i)%weight * exp(-mol_data%level_list(i)%energy/T)
  end do
end function cdms_calc_partition_my


function cdms_intensity2Aul(intens, freq_Hz, Elow, Eup, gup, T) result(Aul)
  double precision Aul
  double precision, intent(in) :: intens, freq_Hz, Elow, Eup, gup, T
  double precision T0
  ! eq 9 of 1998-Pickett-Submillimeter, millimeter and microwave spectral line catalog
  Aul = intens * (freq_Hz*1D-6)**2 * cdms_calc_partition(T) / gup &
        / (exp(-Elow/T) - exp(-Eup/T)) * 2.7964D-16
end function cdms_intensity2Aul



subroutine test_load_cdms_adhoc
  type(type_molecule_energy_set), pointer :: molecule
  integer i, n
  !
  call load_cdms_mol('/n/Users/fdu/now/', 'cdms_HD.cat', 'cdms_partition_functions.dat', molecule)
  !
  write(*, '(A)') 'Energy levels'
  n = min(50, molecule%n_level)
  do i=1, n
    write(*, '(I4, 2ES16.8)') &
        i, molecule%level_list(i)%energy, &
        molecule%level_list(i)%weight
  end do
  write(*, '(A)') 'Transitions'
  n = min(50, molecule%rad_data%n_transition)
  do i=1, n
    write(*, '(I4, 7ES16.8)') i, &
      molecule%rad_data%list(i)%freq, &
      molecule%rad_data%list(i)%lambda, &
      molecule%rad_data%list(i)%Eup, &
      molecule%rad_data%list(i)%Elow, &
      molecule%rad_data%list(i)%Aul, &
      molecule%rad_data%list(i)%Bul, &
      molecule%rad_data%list(i)%Blu
  end do
  write(*, '(A)') 'Partition function'
  do i=1, nT
    write(*, '(I4, 3ES16.8)') i, T_s(i), &
      lg10Q(i), log10(cdms_calc_partition_my(molecule, T_s(i)))
    ! The output of the above line demonstrate that the partition function Q of
    ! the CDMS database is simply calculated with
    ! $$ Q(T) = \Sigma g_i \exp(-E_i/T). $$
  end do
end subroutine test_load_cdms_adhoc

end module cdms
