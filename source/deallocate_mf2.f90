subroutine deallocate_mf2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate large arrays for MF2
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(AJ)) deallocate(AJ)
  if (allocated(Er)) deallocate(Er)
  if (allocated(ER7)) deallocate(ER7)

  if (allocated(GF)) deallocate(GF)
  if (allocated(GFA)) deallocate(GFA)
  if (allocated(GFB)) deallocate(GFB)
  if (allocated(GG)) deallocate(GG)
  if (allocated(GN)) deallocate(GN)
  if (allocated(GT)) deallocate(GT)

  if (allocated(D)) deallocate(D)
  if (allocated(Es)) deallocate(Es)
  if (allocated(GFu)) deallocate(GFu)
  if (allocated(GGu)) deallocate(GGu)
  if (allocated(GN0)) deallocate(GN0)
  if (allocated(GX)) deallocate(GX)

  if (allocated(GAM7)) deallocate(GAM7)

  return
end subroutine deallocate_mf2
