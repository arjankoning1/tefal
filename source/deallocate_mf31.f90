subroutine deallocate_mf31
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF31
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(b31)) deallocate(b31)

  return
end subroutine deallocate_mf31
