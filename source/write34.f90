subroutine write34
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF34
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
!   AWR        ! standard mass parameter
!   MAT        ! MAT number
!   blank2     ! blank string
!   FEND       ! ENDF - 6 format
!   mtexist    ! flag for existence of MT - number
!   SEND       ! ENDF - 6 format
!   ZA         ! standard charge parameter
! All global variables
!   b34            ! covariance matrix element
!   LB34           ! flag for meaning of numbers
!   LS34           ! symmetry flag
!   LTT34          ! representation
!   LVT34          ! specification of transformation matrix
!   MAT34          ! material number for MF34
!   MT34           ! MT number for MF34
!   NE34           ! number of energies in energy array
!   NI34           ! number of NI - type sub - subsections
!   NL34           ! number of Legendre coefficients with covariances
!   NL341          ! number of Legendre coefficients with covariances
!   NMT34          ! total number of MT sections
!   NT34           ! total number of entries
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: L1                        ! integration limits
  integer   :: L2                        ! integration limits
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  real(sgl) :: x(numencovtot)            ! help variable
!
! ***************************** Write MF34 *****************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 34
  NS = 0
  open (unit = 2, file = 'MF34', status = 'replace')
  MT = 2
  if ( .not. mtexist(MF, MT)) return
  call hrwrite(ZA, AWR, LVT34, LTT34, 0, NMT34, MAT, MF, MT, NS)
  call hrwrite(0., 0., MAT34, MT34, NL34, NL341, MAT, MF, MT, NS)
  do L1 = 1, NL34
    do L2 = L1, NL341
      call hrwrite(0., 0., L1, L2, 0, NI34(L1, L2), MAT, MF, MT, NS)
      call hrwrite(0., 0., LS34, LB34, NT34(L1, L2), NE34(L1, L2), MAT, MF, MT, NS)
!
! 1. Covariance data per L1,L2 set
!
! xwrite: subroutine to write real value block
!
      N = NT34(L1, L2)
      do i = 1, N
        x(i) = b34(L1, L2, i)
      enddo
      call xwrite(N, x, MAT, MF, MT, NS)
    enddo
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write34
! Copyright A.J. Koning 2021
