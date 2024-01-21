subroutine read15(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF15 from existing ENDF-6 data library
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
  use A1_error_handling_mod
!
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   adoptfile    ! name of library for MF information (per MT)
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF12_15
!   E15          ! incident energy (in ENDF - 6 format)
!   E15ge        ! secondary energy
!   EPy          ! incident energy for probabilities (in ENDF - 6 format)
!   ge           ! gamma distribution
!   INTER15      ! interpolation scheme
!   INTER15g     ! interpolation scheme
!   INTER15ge    ! interpolation scheme
!   NBT15        ! separation value for interpolation scheme
!   NBT15g       ! separation value for interpolation scheme
!   NBT15ge      ! separation value for interpolation scheme
!   NE15g        ! number of incident energies for distribution
!   NP15         ! number of incident energies
!   NP15ge       ! number of secondary energy point
!   NR15         ! number of interpolation ranges
!   NR15g        ! number of interpolation ranges
!   NR15ge       ! number of interpolation ranges
!   Pg           ! probability (in ENDF - 6 format)
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: iE        ! energy counter
  integer           :: ii        ! counter
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: k         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: nlin      ! number of lines
!
! ********* Read data from particular MT-numbers from MF15 files *******
!
  MF = 15
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
!
! Read until MT is found
!
! Incident energies and probabilities
!
  MTstr = '     '
  write(MTstr(1:2), '(i2)') MF
  write(MTstr(3:5), '(i3)') MT
  do
    read(3, '(a80)', iostat = istat) string
    if (istat /= 0) then
      close (unit = 3)
      return
    endif
    if (string(71:75) == MTstr) exit
  enddo
  read(string(45:55), '(i11)', iostat = istat) NK(MF, MT)
  do k = 1, NK(MF, MT)
    read(3, '(44x, 2i11)', iostat = istat) NR15(k), NP15(k)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    nlin = 1 + (NR15(k) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBT15(k, j), INTER15(k, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    nlin = 1 + (NP15(k) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (Epy(k, j), Pg(k, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
!
! Secondary gamma distribution
!
    read(3, '(44x, 2i11)', iostat = istat) NR15g(k), NE15g(k)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    nlin = 1 + (NR15g(k) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBT15g(k, j), INTER15g(k, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    do iE = 1, NE15g(k)
      read(3, '(11x, e11.6, 22x, 2i11)', iostat = istat) E15(k, iE), NR15ge(k, iE), NP15ge(k, iE)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      nlin = 1 + (NR15ge(k, iE) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6i11)', iostat = istat) (NBT15ge(k, iE, j), INTER15ge(k, iE, j), j = ii + 1, ii + 3)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      enddo
      nlin = 1 + (NP15ge(k, iE) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6e11.6)', iostat = istat) (E15ge(k, iE, j), ge(k, iE, j), j = ii+1, ii+3)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      enddo
    enddo
  enddo
  close (unit = 3)
  return
end subroutine read15
! Copyright A.J. Koning 2021
