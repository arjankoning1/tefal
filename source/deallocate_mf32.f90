subroutine deallocate_mf32
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF32
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(covdigit)) deallocate(covdigit)
  if (allocated(covix32)) deallocate(covix32)

  if (allocated(ISR)) deallocate(ISR)
  if (allocated(LCOMP)) deallocate(LCOMP)
  if (allocated(MLS)) deallocate(MLS)
  if (allocated(NJS32)) deallocate(NJS32)
  if (allocated(NLS32)) deallocate(NLS32)

  if (allocated(AJ32)) deallocate(AJ32)
  if (allocated(D32)) deallocate(D32)
  if (allocated(DAP)) deallocate(DAP)
  if (allocated(GF32)) deallocate(GF32)
  if (allocated(GG32)) deallocate(GG32)
  if (allocated(GNO32)) deallocate(GNO32)
  if (allocated(GX32)) deallocate(GX32)

  if (allocated(b32)) deallocate(b32)
  if (allocated(b32URR)) deallocate(b32URR)

  return
end subroutine deallocate_mf32
