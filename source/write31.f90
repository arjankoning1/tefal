subroutine write31
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF31
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! Variables for initialization of ENDF format
!   AWR        ! standard mass parameter
!   blank2     ! blank string
!   MAT        ! MAT number
!   mfexist    ! flag for existence of MF - number
!   mtexist    ! flag for existence of MT - number
!   FEND       ! ENDF - 6 format
!   SEND       ! ENDF - 6 format
!   ZA         ! standard charge parameter
! All global variables
!   b31           ! covariance matrix element
!   LB31          ! flag for meaning of numbers
!   LS31          ! symmetry flag
!   MTL           ! lumped reaction identifier
!   NC31          ! number of NC - type sub - subsections
!   NE31          ! number of energies in energy array
!   NI31          ! number of NI - type sub - subsections
!   NL31          ! number of subsections
!   NT31          ! total number of entries
!
! *** Declaration of local data
!
  implicit none
  integer :: MF              ! MF-number
  integer :: MT              ! MT-number
  integer :: NS              ! line number
!
! ***************************** Write MF31 *****************************
!
! hrwrite: subroutine to write header with real values
! xwrite : subroutine to write real value block
!
  MF = 31
  MT = 452
  if ( .not. mtexist(MF, MT)) return
  NS = 0
  open (unit = 2, file = 'MF31', status = 'replace')
  call hrwrite(ZA, AWR, 0, MTL(MT), 0, NL31, MAT, MF, MT, NS)
  call hrwrite(0., 0., 0, MT, NC31, NI31, MAT, MF, MT, NS)
  call hrwrite(0., 0., LS31, LB31, NT31, NE31, MAT, MF, MT, NS)
  call xwrite(NT31, b31, MAT, MF, MT, NS)
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write31
! Copyright A.J. Koning 2021
