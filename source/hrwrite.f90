subroutine hrwrite(r1, r2, i1, i2, i3, i4, mAtNUM, MF, MT, NS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write header with real values
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
!   sgl     ! single precision kind
! Variables for initialization of ENDF format
!   HEAD    ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  character(len=11) :: c1        ! help variable
  character(len=11) :: c2        ! help variable
  character(len=11) :: endf      ! function for a number in ENDF-6 format
  integer           :: i1        ! value
  integer           :: i2        ! value
  integer           :: i3        ! value
  integer           :: i4        ! value
  integer           :: MATnum    ! material number
  integer           :: MF        ! MF-number
  integer           :: NS        ! line number
  integer           :: MT        ! MT-number
  real(sgl)         :: r1        ! value
  real(sgl)         :: r2        ! value
!
! ***************************** Write block ****************************
!
! Write values
!
  c1 = endf(r1)
  c2 = endf(r2)
  write(2, fmt = HEAD) c1, c2, i1, i2, i3, i4, MATnum, MF, MT, NS
  return
end subroutine hrwrite
! Copyright A.J. Koning 2021
