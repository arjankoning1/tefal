subroutine read31
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF31 from existing ENDF-6 data library
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
! Definition of single and double precision variables
!   dbl           ! double precision kind
! All global variables
!   numencovtot   ! number of energies for covariances
! Variables for input of specific ENDF data
!   adopt         ! logical for existence of MF information (per MT)
!   adoptfile     ! name of library for MF information (per MT)
! Variables for initialization of ENDF format
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
! Variables for MF31_40
!   b31           ! covariance matrix element
!   LB31          ! flag for meaning of numbers
!   LS31          ! symmetry flag
!   MTL           ! lumped reaction identifier
!   NC31          ! number of NC - type sub - subsections
!   NE31          ! number of energies in energy array
!   NI31          ! number of NI - type sub - subsections
!   NL31          ! number of subsections
!   NT31          ! total number of entries
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: ii        ! counter
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: nlin      ! number of lines
  real(dbl)         :: bdb(6)    ! covariance matrix element
!
! ****************** Read covariance parameters from MF31 **************
!
  MF = 31
  MT = 452
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
!
! Read until MF31 is found
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
  read(string, '(33x, i11, 11x, i11)') MTL(MT), NL31
  read(3, '(44x, 2i11)', iostat = istat) NC31, NI31
  if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  read(3, '(22x, 4i11)', iostat = istat) LS31, LB31, NT31, NE31
  if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  if (NT31 > 10 * numencovtot) then
    write(*, '(" TEFAL-error: MF31 Increase numencovtot")')
    stop
  endif
  nlin = 1 + (NT31 - 1) / 6
  do i = 1, nlin
    ii = 6 * (i - 1)
    read(3, '(6e11.6)', iostat = istat) (bdb(j), j = 1, 6)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    do j = ii + 1, ii + 6
      if (abs(bdb(j - ii)) >= 1.e-38) then
        b31(j) = real(bdb(j - ii))
      else
        b31(j) = 0.
      endif
    enddo
  enddo
  close (unit = 3)
  mfexist(MF) = .true.
  mtexist(MF, MT) = .true.
  return
end subroutine read31
! Copyright A.J. Koning 2021
