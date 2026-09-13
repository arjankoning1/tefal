subroutine allocate_ddx
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for double-differential spectra
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  Nddx = 0
  if (.not. flagtabddx) return

  allocate(ddxemis(1:2,1:Nenspec,1:numddx,0:numen2))
  allocate(f0ddx(1:2,1:Nenspec,1:numddx,0:numen2))

  ddxemis = 0.
  f0ddx = 0.

  return
end subroutine allocate_ddx
