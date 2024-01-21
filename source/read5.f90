subroutine read5(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF5 from existing ENDF-6 data library
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
! Variables for initialization of ENDF format
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF5
!   E5           ! incident energy for MF5 (in ENDF - 6 format)
!   E5p          ! incident energy for which tabulated distribution is given
!   EFH          ! constant used in the energy - dependent fission neutron spectrum
!   EFL          ! constant used in the energy - dependent fission neutron spectrum
!   gE5          ! energy - spectrum values
!   INTER5       ! interpolation scheme
!   INTER5e      ! interpolation scheme
!   INTER5e2     ! interpolation scheme
!   LF           ! flag for energy distribution law
!   NBT5         ! separation value for interpolation scheme
!   NBT5e        ! separation value for interpolation scheme
!   NBT5e2       ! separation value for interpolation scheme
!   NE5e         ! number of incident energies for distribution
!   NF           ! number of secondary energy points
!   NP5          ! number of incident energies
!   NR5          ! number of interpolation ranges
!   NR5e         ! number of interpolation ranges
!   NR5e2        ! number of interpolation ranges
!   pE           ! fractional part of cross section
!   TM5          ! Effective (7) or maximum (12) temperature FNS parameter, T(E)
!   U            ! constant for upper energy limit
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
! ********* Read data from particular MT-numbers from MF5 files ********
!
  MF = 5
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
!
! Read until MT is found
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
    read(3, '(e11.6, 22x, 3i11)', iostat = istat) U(k), LF(k), NR5(k), NP5(k)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    nlin = 1 + (NR5(k) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBT5(k, j), INTER5(k, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    nlin = 1 + (NP5(k) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (E5p(k, j), pE(k, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
!
! LF=1: Arbitrary tabulated function
!
    if (LF(k) == 1) then
      read(3, '(44x, 2i11)', iostat = istat) NR5e(k), NE5e(k)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      nlin = 1 + (NR5e(k) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6i11)', iostat = istat) (NBT5e(k, j), INTER5e(k, j), j = ii+1, ii+3)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      enddo
      do iE = 1, NE5e(k)
        read(3, '(11x, e11.6, 22x, 2i11)', iostat = istat) E5(k, iE), NR5e2(k, iE), NF(k, iE)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        nlin = 1 + (NR5e2(k, iE) - 1) / 3
        do i = 1, nlin
          ii = 3 * (i - 1)
          read(3, '(6i11)', iostat = istat) (NBT5e2(k, iE, j), INTER5e2(k, iE, j), j = ii + 1, ii + 3)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        enddo
        nlin = 1 + (2 * NF(k, iE) - 1) / 6
        do i = 1, nlin
          ii = 6 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (gE5(k, iE, j), j = ii+1, ii+6)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        enddo
      enddo
    endif
    if (LF(k) == 5 .or. LF(k) == 7 .or. LF(k) == 9 .or. LF(k) == 12) then
      if (LF(k) == 5) read(3, '(//)')
      read(3, '(2e11.6, 22x, 2i11)', iostat = istat) EFL, EFH, NR5e(k), NE5e(k)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      nlin = 1 + (NR5e(k) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6i11)', iostat = istat) (NBT5e(k, j), INTER5e(k, j), j = ii+1, ii+3)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      enddo
      nlin = 1 + (2 * NE5e(k) - 1) / 6
      do i = 1, nlin
        ii = 6 * (i - 1)
        read(3, '(6e11.6)', iostat = istat) (TM5(k, j), j = ii+1, ii+6)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      enddo
    endif
  enddo
  mtexist(MF, MT) = .true.
  mfexist(MF) = .true.
  close (unit = 3)
  return
end subroutine read5
! Copyright A.J. Koning 2021
