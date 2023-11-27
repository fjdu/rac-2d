module lamda

use data_struct
use trivials
use phy_const

implicit none

contains

subroutine load_moldata_LAMDA(filename, mol_data)
  character(len=*), intent(in) :: filename
  type(type_molecule_energy_set), pointer, intent(inout) :: mol_data
  character(len=512) strtmp
  integer, parameter :: nstr_split = 64
  character(len=32), dimension(nstr_split) :: str_split
  character(len=8), parameter :: strfmt_row = '(A512)'
  character(len=8), parameter :: strfmt_float = '(F16.0)'
  character(len=8), parameter :: strfmt_int = '(I6)'
  integer i, j, k, fU, nout !, iup_, ilow_
  integer ios, iup, ilow
  integer n_T_, n_transition_
  character(len=const_len_qnum*2), dimension(:), allocatable :: qnum
  !
  write(*, '(2A)') 'Reading from file (LAMDA format): ', filename
  call openFileSequentialRead(fU, filename, 99999, getu=1)
  ! Get molecule name.
  read(fU,'(A1)') strtmp
  read(fU, strfmt_row) strtmp
  call split_str_by_space(strtmp, str_split, nstr_split, nout)
  mol_data%name_molecule = trim(str_split(1))
  ! Get energy level list
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU,'(I32)') mol_data%n_level
  read(fU,'(A1)') strtmp
  allocate(mol_data%level_list(mol_data%n_level), &
           mol_data%f_occupation(mol_data%n_level))
  allocate(qnum(mol_data%n_level))
  do i=1, mol_data%n_level
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_float, iostat=ios) mol_data%level_list(i)%energy
    if (ios .ne. 0) then
      write(*, '(A)') 'In load_moldata_LAMDA:'
      write(*, '(I6, 2X, A)') i, str_split(2)
      write(*, '(A, I6)') 'IOSTAT = ', ios
      call error_stop()
    end if
    read(str_split(3), strfmt_float, iostat=ios) mol_data%level_list(i)%weight
    if (ios .ne. 0) then
      write(*, '(A)') 'In load_moldata_LAMDA:'
      write(*, '(I6, 2X, A)') i, str_split(3)
      write(*, '(A, I6)') 'IOSTAT = ', ios
      call error_stop()
    end if
    !
    nout = len_trim(strtmp)
    qnum(i) = strtmp(max(1, nout-const_len_qnum*2+3):nout)
    !
  end do
  !
  ! Get radiative transitions
  allocate(mol_data%rad_data)
  read(fU,'(A1)') strtmp
  read(fU,'(I8)') mol_data%rad_data%n_transition
  read(fU,'(A1)') strtmp
  allocate(mol_data%rad_data%list(mol_data%rad_data%n_transition))
  do i=1, mol_data%rad_data%n_transition
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_int) mol_data%rad_data%list(i)%iup
    read(str_split(3), strfmt_int) mol_data%rad_data%list(i)%ilow
    read(str_split(4), strfmt_float) mol_data%rad_data%list(i)%Aul
    !read(str_split(5), strfmt_float) mol_data%rad_data%list(i)%freq
    !read(str_split(6), strfmt_float) mol_data%rad_data%list(i)%Eup
    !
    ! The frequency in the LAMDA database may be incorrect, so here I recompute from the
    ! energy difference.  The result is in Hz.
    iup  = mol_data%rad_data%list(i)%iup
    ilow = mol_data%rad_data%list(i)%ilow
    !
    mol_data%rad_data%list(i)%freq = phy_SpeedOfLight_CGS * &
      (mol_data%level_list(iup)%energy - &
       mol_data%level_list(ilow)%energy)
    !
    mol_data%rad_data%list(i)%Eup  = mol_data%level_list(iup)%energy * phy_cm_1_2K
    mol_data%rad_data%list(i)%Elow = mol_data%level_list(ilow)%energy * phy_cm_1_2K
    !
    mol_data%rad_data%list(i)%qnum = trim(adjustl(qnum(iup))) // &
        ' -> ' // trim(adjustl(qnum(ilow)))
  end do
  deallocate(qnum)
  !
  ! Convert the energy unit into Kelvin from cm-1
  mol_data%level_list%energy = mol_data%level_list%energy * phy_cm_1_2K
  !
  ! Lambda in Angstrom
  mol_data%rad_data%list%lambda = &
      phy_SpeedOfLight_SI/mol_data%rad_data%list%freq*1D10
  !
  mol_data%rad_data%list%Bul = mol_data%rad_data%list%Aul / &
    ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * &
     (mol_data%rad_data%list%freq)**3)
  do i=1, mol_data%rad_data%n_transition
    j = mol_data%rad_data%list(i)%iup
    k = mol_data%rad_data%list(i)%ilow
    mol_data%rad_data%list(i)%Blu = mol_data%rad_data%list(i)%Bul * &
        mol_data%level_list(j)%weight / mol_data%level_list(k)%weight
  end do
  !
  ! Collisional transitions
  allocate(mol_data%colli_data)
  ! Get the number of collisional partners
  read(fU,'(A1)') strtmp
  read(fU,'(I4)') mol_data%colli_data%n_partner
  allocate(mol_data%colli_data%list(mol_data%colli_data%n_partner))
  do i=1, mol_data%colli_data%n_partner
    ! Get the name of partner
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    mol_data%colli_data%list(i)%name_partner = trim(str_split(4))
    if (mol_data%colli_data%list(i)%name_partner .eq. 'electron') then
      mol_data%colli_data%list(i)%name_partner = 'e'
    end if
    ! Get the number of transitions and temperatures
    read(fU,'(A1)') strtmp
    read(fU,'(I8)') mol_data%colli_data%list(i)%n_transition
    read(fU,'(A1)') strtmp
    read(fU,'(I4)') mol_data%colli_data%list(i)%n_T
    !
    ! Name too long...
    n_transition_ = mol_data%colli_data%list(i)%n_transition
    n_T_ = mol_data%colli_data%list(i)%n_T
    if ((n_T_+3) .gt. nstr_split) then
      write(*,*) 'The number of different temperatures is too large!'
      write(*,*) 'nstr_split = ', nstr_split
      write(*,*) 'Change nstr_split of the source code to a higher value.'
      call error_stop()
    end if
    !
    allocate(mol_data%colli_data%list(i)%iup(n_transition_), &
             mol_data%colli_data%list(i)%ilow(n_transition_), &
             mol_data%colli_data%list(i)%T_coll(n_T_), &
             mol_data%colli_data%list(i)%Cul(n_T_, n_transition_))
    !
    ! Get the list of temperatures
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    do j=1, n_T_
      read(str_split(j), strfmt_float) mol_data%colli_data%list(i)%T_coll(j)
    end do
    ! Get the collision coefficients
    read(fU,'(A1)') strtmp
    do j=1, n_transition_
      read(fU, strfmt_row) strtmp
      call split_str_by_space(strtmp, str_split, nstr_split, nout)
      read(str_split(2), strfmt_int) mol_data%colli_data%list(i)%iup(j)
      read(str_split(3), strfmt_int) mol_data%colli_data%list(i)%ilow(j)
      do k=1, n_T_
        read(str_split(3+k), strfmt_float) mol_data%colli_data%list(i)%Cul(k, j)
        !iup_ = mol_data%colli_data%list(i)%iup(j)
        !ilow_ = mol_data%colli_data%list(i)%ilow(j)
      end do
    end do
  end do
  !
  close(fU)
  ! Test the results
  ! write(*,*) mol_data%name_molecule
  ! do i=1, mol_data%n_level
  !   write(*,*) i, mol_data%level_list(i)%energy, mol_data%level_list(i)%weight
  ! end do
  ! do i=1, mol_data%rad_data%n_transition
  !   write(*,*) i, mol_data%rad_data%list(i)%iup, &
  !     mol_data%rad_data%list(i)%ilow, &
  !     mol_data%rad_data%list(i)%Aul, &
  !     mol_data%rad_data%list(i)%Bul, &
  !     mol_data%rad_data%list(i)%Blu, &
  !     mol_data%rad_data%list(i)%freq, &
  !     mol_data%rad_data%list(i)%Eup
  ! end do
  ! do i=1, mol_data%colli_data%n_partner
  !   write(*,*) i, mol_data%colli_data%list(i)%name_partner
  !   write(*,*) mol_data%colli_data%list(i)%T_coll
  !   do j=1, mol_data%colli_data%list(i)%n_transition
  !     write(*,*) i, j, mol_data%colli_data%list(i)%iup(j), &
  !       mol_data%colli_data%list(i)%ilow(j), &
  !       mol_data%colli_data%list(i)%Cul(:, j)
  !   end do
  ! end do
end subroutine load_moldata_LAMDA


end module lamda
