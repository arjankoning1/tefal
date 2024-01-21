subroutine tefalinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! arrayinitial: subroutine for initialization of arrays
! reacinitial : subroutine for initialization of nuclear reaction info
! endfinitial : subroutine for initialization of ENDF-6 formats
! mainout     : subroutine with main output
!
  call arrayinitial
  call reacinitial
  call endfinitial
  call mainout
  return
end subroutine tefalinitial
! Copyright A.J. Koning 2021
