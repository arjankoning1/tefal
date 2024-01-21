subroutine write14(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF14
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
!   AWR       ! standard mass parameter
!   blank2    ! blank string
!   MAT       ! MAT number
!   SEND      ! ENDF - 6 format
!   ZA        ! standard charge parameter
! Variables for ENDF format
!   NK        ! number of subsections
! Variables for MF12_15
!   LI14      ! isotropy flag
!
! *** Declaration of local data
!
  implicit none
  integer :: MT              ! MT-number
  integer :: MF              ! MF-number
  integer :: NS              ! line number
!
! ***************************** Write MF14 *****************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 14
  NS = 0
  call hrwrite(ZA, AWR, LI14(MT), 0, NK(MF, MT), 0, MAT, MF, MT, NS)
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write14
! Copyright A.J. Koning 2021
