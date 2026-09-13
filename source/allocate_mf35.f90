subroutine allocate_mf35
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF35
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: nblock35
!
! One additional block may be added by make35 for energies above 20 MeV.
!
  nblock35 = numencov35 + 1
!
! MF35 arrays
!
  allocate(NE35(nblock35))
  allocate(NT35(nblock35))

  allocate(E35b(nblock35))
  allocate(E35e(nblock35))

  allocate(b35(nblock35,numencov35))
!
! Initialization
!
  NE35 = 0
  NT35 = 0

  E35b = 0.
  E35e = 0.
  b35 = 0.

  LB35 = 0
  LS35 = 0

  return
end subroutine allocate_mf35
