subroutine write35
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF35
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
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numencovtot    ! number of energies for covariances
! Variables for initialization of ENDF format
!   AWR            ! standard mass parameter
!   blank2         ! blank string
!   FEND           ! ENDF - 6 format
!   MAT            ! MAT number
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
!   SEND           ! ENDF - 6 format
!   ZA             ! standard charge parameter
! Variables for ENDF format
!   NK             ! number of subsections
! Variables for MF31_40
!   b35            ! covariance matrix element
!   E35b           ! start energy of block
!   E35e           ! end energy of block
!   LB35           ! flag for meaning of numbers
!   LS35           ! symmetry flag
!   NE35           ! number of energies in energy array
!   NT35           ! total number of entries
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: k                         ! counter
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  real(sgl) :: x(numencovtot)            ! help variable
!
! ***************************** Write MF31 *****************************
!
! hrwrite: subroutine to write header with real values
! xwrite : subroutine to write real value block
!
  MF = 35
  MT = 18
  if ( .not. mtexist(MF, MT)) return
  NS = 0
  open (unit = 2, file = 'MF35', status = 'replace')
  call hrwrite(ZA, AWR, 0, 0, NK(MF, MT), 0, MAT, MF, MT, NS)
  do k = 1, NK(MF, MT)
    call hrwrite(E35b(k), E35e(k), LS35, LB35, NT35(k), NE35(k), MAT, MF, MT, NS)
    N = NT35(k)
    do i = 1, N
      x(i) = b35(k, i)
    enddo
    call xwrite(N, x, MAT, MF, MT, NS)
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write35
! Copyright A.J. Koning 2021
