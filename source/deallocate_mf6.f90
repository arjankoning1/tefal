subroutine deallocate_mf6
  use A0_tefal_mod
  implicit none

  if (allocated(NL)) deallocate(NL)
  if (allocated(leg)) deallocate(leg)

  if (allocated(INTER6ea)) deallocate(INTER6ea)
  if (allocated(INTER6y)) deallocate(INTER6y)

  if (allocated(LANG)) deallocate(LANG)
  if (allocated(LAW)) deallocate(LAW)
  if (allocated(LEP)) deallocate(LEP)
  if (allocated(LIP)) deallocate(LIP)

  if (allocated(NA)) deallocate(NA)
  if (allocated(NBT6ea)) deallocate(NBT6ea)
  if (allocated(NBT6y)) deallocate(NBT6y)
  if (allocated(ND)) deallocate(ND)
  if (allocated(NE6ea)) deallocate(NE6ea)
  if (allocated(NEP)) deallocate(NEP)
  if (allocated(NP6y)) deallocate(NP6y)
  if (allocated(NR6ea)) deallocate(NR6ea)
  if (allocated(NR6y)) deallocate(NR6y)
  if (allocated(NW)) deallocate(NW)

  if (allocated(AWP)) deallocate(AWP)
  if (allocated(b6)) deallocate(b6)
  if (allocated(E6)) deallocate(E6)
  if (allocated(Ey)) deallocate(Ey)
  if (allocated(Y)) deallocate(Y)
  if (allocated(ZAP)) deallocate(ZAP)

  if (allocated(flagrec)) deallocate(flagrec)
  if (allocated(b6gam)) deallocate(b6gam)
  if (allocated(b6rec)) deallocate(b6rec)

  return
end subroutine deallocate_mf6
