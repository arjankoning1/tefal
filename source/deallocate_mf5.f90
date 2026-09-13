subroutine deallocate_mf5
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF5
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(INTER5)) deallocate(INTER5)
  if (allocated(INTER5e)) deallocate(INTER5e)
  if (allocated(INTER5e2)) deallocate(INTER5e2)

  if (allocated(LF)) deallocate(LF)

  if (allocated(NBT5)) deallocate(NBT5)
  if (allocated(NBT5e)) deallocate(NBT5e)
  if (allocated(NBT5e2)) deallocate(NBT5e2)

  if (allocated(NE5e)) deallocate(NE5e)
  if (allocated(NF)) deallocate(NF)
  if (allocated(NP5)) deallocate(NP5)
  if (allocated(NR5)) deallocate(NR5)
  if (allocated(NR5e)) deallocate(NR5e)
  if (allocated(NR5e2)) deallocate(NR5e2)

  if (allocated(E5)) deallocate(E5)
  if (allocated(E5p)) deallocate(E5p)
  if (allocated(gE5)) deallocate(gE5)
  if (allocated(pE)) deallocate(pE)
  if (allocated(TM5)) deallocate(TM5)
  if (allocated(U)) deallocate(U)

  return
end subroutine deallocate_mf5
