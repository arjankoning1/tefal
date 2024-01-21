subroutine make4fission
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF4 for fission
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
!   mfexist    ! flag for existence of MF - number
!   mtexist    ! flag for existence of MT - number
! Variables for MF4
!   LCT        ! LAB / CM flag
!   LI4        ! isotropy flag
!   LTT        ! representation
!   LVT        ! specification of transformation matrix
!
! *** Declaration of local data
!
  implicit none
  integer :: MT              ! MT-number
!
! **************************** Make MF4 ********************************
!
  MT = 18
!
! ENDF-6 parameters
!
  LVT = 0
  LTT = 0
  LI4 = 1
  LCT = 1
  mtexist(4, MT) = .true.
  mfexist(4) = .true.
  return
end subroutine make4fission
! Copyright A.J. Koning 2021
