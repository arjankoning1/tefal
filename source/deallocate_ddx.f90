subroutine deallocate_ddx
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for double-differential spectra
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(ddxemis)) deallocate(ddxemis)
  if (allocated(f0ddx)) deallocate(f0ddx)

  return
end subroutine deallocate_ddx
