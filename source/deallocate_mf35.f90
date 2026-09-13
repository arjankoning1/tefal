subroutine deallocate_mf35
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF35
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(NE35)) deallocate(NE35)
  if (allocated(NT35)) deallocate(NT35)

  if (allocated(E35b)) deallocate(E35b)
  if (allocated(E35e)) deallocate(E35e)
  if (allocated(b35)) deallocate(b35)

  return
end subroutine deallocate_mf35
