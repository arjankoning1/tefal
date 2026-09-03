subroutine deallocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
  if (allocated(Ehist)) deallocate(Ehist)
  if (allocated(f0ex)) deallocate(f0ex)
  if (allocated(specexcl)) deallocate(specexcl)
!
  return
end subroutine deallocate_arrays
