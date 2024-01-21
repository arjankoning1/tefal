subroutine rwrite(r1, r2, r3, r4, r5, r6, MATNUM, MF, MT, NS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write line with real values
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
!   sgl    ! single precision kind
! Variables for initialization of ENDF format
!   VALUE  ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  character(len=11) :: c1        ! help variable
  character(len=11) :: c2        ! help variable
  character(len=11) :: c3        ! character
  character(len=11) :: c4        ! character
  character(len=11) :: c5        ! character
  character(len=11) :: c6        ! character
  character(len=11) :: endf      ! function for a number in ENDF-6 format
  integer           :: MATnum    ! material number
  integer           :: MF        ! MF-number
  integer           :: NS        ! line number
  integer           :: MT        ! MT-number
  real(sgl)         :: r1        ! value
  real(sgl)         :: r2        ! value
  real(sgl)         :: r3        ! value
  real(sgl)         :: r4        ! value
  real(sgl)         :: r5        ! value
  real(sgl)         :: r6        ! value
!
! ***************************** Write block ****************************
!
! Write values
!
  c1 = endf(r1)
  c2 = endf(r2)
  c3 = endf(r3)
  c4 = endf(r4)
  c5 = endf(r5)
  c6 = endf(r6)
  write(2, fmt = VALUE) c1, c2, c3, c4, c5, c6, MATnum, MF, MT, NS
  return
end subroutine rwrite
! Copyright A.J. Koning 2021
