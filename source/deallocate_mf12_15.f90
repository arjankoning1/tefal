subroutine deallocate_mf12_15
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for MF12-15
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
  if (allocated(INTERg)) deallocate(INTERg)
  if (allocated(LG12)) deallocate(LG12)
  if (allocated(LO12)) deallocate(LO12)
  if (allocated(NBTg)) deallocate(NBTg)
  if (allocated(NPg)) deallocate(NPg)
  if (allocated(NRg)) deallocate(NRg)
  if (allocated(LP12)) deallocate(LP12)
  if (allocated(LPg)) deallocate(LPg)
  if (allocated(LFg)) deallocate(LFg)
  if (allocated(NS12)) deallocate(NS12)
  if (allocated(NT12)) deallocate(NT12)

  if (allocated(E12)) deallocate(E12)
  if (allocated(Eg)) deallocate(Eg)
  if (allocated(Egk)) deallocate(Egk)
  if (allocated(ES12)) deallocate(ES12)
  if (allocated(Esk)) deallocate(Esk)
  if (allocated(ESNS)) deallocate(ESNS)
  if (allocated(TP12)) deallocate(TP12)
  if (allocated(xsgtotyield)) deallocate(xsgtotyield)
  if (allocated(xsgyield)) deallocate(xsgyield)

  if (allocated(E13)) deallocate(E13)
  if (allocated(xsg)) deallocate(xsg)
  if (allocated(xsgtot)) deallocate(xsgtot)

  if (allocated(LI14)) deallocate(LI14)

  if (allocated(INTER15)) deallocate(INTER15)
  if (allocated(INTER15g)) deallocate(INTER15g)
  if (allocated(INTER15ge)) deallocate(INTER15ge)
  if (allocated(NBT15)) deallocate(NBT15)
  if (allocated(NBT15g)) deallocate(NBT15g)
  if (allocated(NBT15ge)) deallocate(NBT15ge)
  if (allocated(NE15g)) deallocate(NE15g)
  if (allocated(NP15)) deallocate(NP15)
  if (allocated(NP15ge)) deallocate(NP15ge)
  if (allocated(NR15)) deallocate(NR15)
  if (allocated(NR15g)) deallocate(NR15g)
  if (allocated(NR15ge)) deallocate(NR15ge)

  if (allocated(E15)) deallocate(E15)
  if (allocated(E15ge)) deallocate(E15ge)
  if (allocated(EPy)) deallocate(EPy)
  if (allocated(ge)) deallocate(ge)
  if (allocated(Pg)) deallocate(Pg)
!
  return
end subroutine deallocate_mf12_15
