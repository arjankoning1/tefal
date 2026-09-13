subroutine allocate_mf4
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF4
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
! Arrays for adopted/read MF4 data
!
  allocate(NL(numen4))
  allocate(NL4r(numen4))
  allocate(NP4r(numen4+1))
  allocate(NR4r(numen4+1))

  allocate(E4hr(numen4+1))
  allocate(E4r(numen4))
  allocate(f4r(numen4+1,numang+3))
  allocate(legr(numen4,0:numl))
  allocate(x4r(numen4+1,numang+3))
!
! Arrays for generated MF4 data
!
  allocate(INTER4(numen4,numint))
  allocate(INTERh(numint))

  allocate(NBT4(numen4,numint))
  allocate(NBTh(numint))

  allocate(NP4(numen4))
  allocate(NR4(numen4))

  allocate(E4(numen4))
  allocate(E4h(numen4+1))
  allocate(f4(numen4+1,numang+3))
  allocate(leg(numen4,0:numl))
  allocate(x4(numen4+1,numang+3))
!
! Initialization
!
  NL = 0
  NL4r = 0
  NP4r = 0
  NR4r = 0

  E4hr = 0.
  E4r = 0.
  f4r = 0.
  legr = 0.
  x4r = 0.

  INTER4 = 0
  INTERh = 0

  NBT4 = 0
  NBTh = 0

  NP4 = 0
  NR4 = 0

  E4 = 0.
  E4h = 0.
  f4 = 0.
  leg = 0.
  x4 = 0.
!
  return
end subroutine allocate_mf4
