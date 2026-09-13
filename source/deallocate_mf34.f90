subroutine deallocate_mf34
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF34
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(NE34)) deallocate(NE34)
  if (allocated(NI34)) deallocate(NI34)
  if (allocated(NT34)) deallocate(NT34)
  if (allocated(b34)) deallocate(b34)

  return
end subroutine deallocate_mf34
