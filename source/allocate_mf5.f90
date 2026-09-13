subroutine allocate_mf5
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF5
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: nen5
!
! One additional point may be added at EMAX.
!
  nen5 = numenin + 1
!
! Interpolation information
!
  allocate(INTER5(numsecea,numint))
  allocate(INTER5e(numsecea,numint))
  allocate(INTER5e2(numsecea,nen5,numint))

  allocate(LF(numsecea))

  allocate(NBT5(numsecea,numint))
  allocate(NBT5e(numsecea,numint))
  allocate(NBT5e2(numsecea,nen5,numint))

  allocate(NE5e(numsecea))
  allocate(NF(numsecea,nen5))
  allocate(NP5(numsecea))
  allocate(NR5(numsecea))
  allocate(NR5e(numsecea))
  allocate(NR5e2(numsecea,nen5))
!
! Energy and spectrum information
!
  allocate(E5(numsecea,nen5))

  allocate(E5p(numsecea,numenin))
  allocate(pE(numsecea,numenin))

  allocate(gE5(numsecea,nen5,10*numen2))

  allocate(TM5(numsecea,2*nen5))

  allocate(U(numsecea))
!
! Initialization
!
  INTER5 = 0
  INTER5e = 0
  INTER5e2 = 0

  LF = 0

  NBT5 = 0
  NBT5e = 0
  NBT5e2 = 0

  NE5e = 0
  NF = 0
  NP5 = 0
  NR5 = 0
  NR5e = 0
  NR5e2 = 0

  E5 = 0.
  E5p = 0.
  gE5 = 0.
  pE = 0.
  TM5 = 0.
  U = 0.

  EFL = 0.
  EFH = 0.

  return
end subroutine allocate_mf5
