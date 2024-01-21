subroutine make35
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF35
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
! All global variables
!   nummt       ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt       ! logical for existence of MF information (per MT)
! Variables for reaction initialization
!   enincmax    ! maximal incident energy
! Variables for initialization of ENDF format
!   mfexist     ! flag for existence of MF - number
!   mtexist     ! flag for existence of MT - number
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF1
!   EMAX        ! upper limit of energy range for evaluation
! Variables for MF31_40
!   b35         ! covariance matrix element
!   E35b        ! start energy of block
!   E35e        ! end energy of block
!   NE35        ! number of energies in energy array
!   NT35        ! total number of entries
!
! *** Declaration of local data
!
  implicit none
  integer :: j               ! counter
  integer :: k               ! counter
  integer :: MF              ! MF-number
  integer :: MT              ! MT-number
!
! ****************** Make covariance parameters for MF35 ***************
!
! read35  : subroutine to read MF35 from existing ENDF-6 data library
!
  MF = 35
  do MT = 1, nummt
    if ( .not. adopt(MF, MT)) cycle
    call read35(MT)
    if (enincmax > 20.) then
      NK(MF, MT) = NK(MF, MT) + 1
      k = NK(MF, MT)
      E35b(k) = E35e(k - 1)
      E35e(k) = EMAX
      NT35(k) = NT35(k - 1)
      NE35(k) = NE35(k - 1)
      do j = 1, NT35(k)
        b35(k, j) = b35(k - 1, j)
      enddo
    endif
    mtexist(MF, MT) = .true.
    mfexist(MF) = .true.
  enddo
  return
end subroutine make35
! Copyright A.J. Koning 2021
