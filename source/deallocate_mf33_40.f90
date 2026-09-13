subroutine deallocate_mf33_40
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF33 and MF40
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(b33read)) deallocate(b33read)
  if (allocated(b33MTread)) deallocate(b33MTread)
  if (allocated(b8read)) deallocate(b8read)

  if (allocated(b33)) deallocate(b33)
  if (allocated(b33MT)) deallocate(b33MT)
  if (allocated(b33ZA)) deallocate(b33ZA)
  if (allocated(b8)) deallocate(b8)

  return
end subroutine deallocate_mf33_40
