subroutine deallocate_mf4
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF4
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(NL)) deallocate(NL)
  if (allocated(NL4r)) deallocate(NL4r)
  if (allocated(NP4r)) deallocate(NP4r)
  if (allocated(NR4r)) deallocate(NR4r)

  if (allocated(E4hr)) deallocate(E4hr)
  if (allocated(E4r)) deallocate(E4r)
  if (allocated(f4r)) deallocate(f4r)
  if (allocated(legr)) deallocate(legr)
  if (allocated(x4r)) deallocate(x4r)

  if (allocated(INTER4)) deallocate(INTER4)
  if (allocated(INTERh)) deallocate(INTERh)

  if (allocated(NBT4)) deallocate(NBT4)
  if (allocated(NBTh)) deallocate(NBTh)

  if (allocated(NP4)) deallocate(NP4)
  if (allocated(NR4)) deallocate(NR4)

  if (allocated(E4)) deallocate(E4)
  if (allocated(E4h)) deallocate(E4h)
  if (allocated(f4)) deallocate(f4)
  if (allocated(leg)) deallocate(leg)
  if (allocated(x4)) deallocate(x4)

  return
end subroutine deallocate_mf4
