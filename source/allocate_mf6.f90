subroutine allocate_mf6
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF6
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: nen6
!
! Maximum number of incident-energy points required for an MF6
! distribution. Extra points may be added at cutoff and EMAX.
!
  nen6 = max(Nenspec, Nenang) + 3
!
! Shared MF4/MF6 work arrays
!
  allocate(NL(nen6))
  allocate(leg(nen6,0:numl))
!
! MF6 arrays
!
  allocate(INTER6ea(numsec,numint))
  allocate(INTER6y(numsec,numint))

  allocate(LANG(numsec))
  allocate(LAW(numsec))
  allocate(LEP(numsec))
  allocate(LIP(numsec))

  allocate(NA(numsec,nen6))
  allocate(NBT6ea(numsec,numint))
  allocate(NBT6y(numsec,numint))
  allocate(ND(numsec,nen6))
  allocate(NE6ea(numsec))
  allocate(NEP(numsec,nen6))
  allocate(NP6y(numsec))
  allocate(NR6ea(numsec))
  allocate(NR6y(numsec))
  allocate(NW(numsec,nen6))

  allocate(AWP(numsec))
  allocate(E6(numsec,nen6))
!
! Yield grids can contain essentially the full incident-energy grid.
!
  allocate(Ey(numsec,numenin))
  allocate(Y(numsec,numenin))

  allocate(ZAP(numsec))
!
! Energy-angle arrays
!
  allocate(b6(numsecea,Nenspec+3,40*numen2))
  allocate(b6gam(Nenspec+3,40*numen2))
  allocate(b6rec(numsec,Nenspec+3,2*numenrec))

  allocate(flagrec(numsec))
!
! Initialization
!
  NL = 0
  leg = 0.

  INTER6ea = 0
  INTER6y = 0
  LANG = 0
  LAW = 0
  LEP = 0
  LIP = 0

  NA = 0
  NBT6ea = 0
  NBT6y = 0
  ND = 0
  NE6ea = 0
  NEP = 0
  NP6y = 0
  NR6ea = 0
  NR6y = 0
  NW = 0

  AWP = 0.
  b6 = 0.
  E6 = 0.
  Ey = 0.
  Y = 0.
  ZAP = 0.

  flagrec = .false.
  b6gam = 0.
  b6rec = 0.

  kpart = 0

  return
end subroutine allocate_mf6
