subroutine deallocate_mf8_10
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate temporary arrays for MF8-10
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(INTER10)) deallocate(INTER10)
  if (allocated(NBT10)) deallocate(NBT10)
  if (allocated(NP10)) deallocate(NP10)
  if (allocated(NR10)) deallocate(NR10)

  if (allocated(INTERZA)) deallocate(INTERZA)
  if (allocated(NBTZA)) deallocate(NBTZA)
  if (allocated(NPZA)) deallocate(NPZA)
  if (allocated(NRZA)) deallocate(NRZA)

  if (allocated(E10)) deallocate(E10)
  if (allocated(xsiso)) deallocate(xsiso)

  if (allocated(E10ZA)) deallocate(E10ZA)
  if (allocated(xsrpZA)) deallocate(xsrpZA)

  return
end subroutine deallocate_mf8_10
