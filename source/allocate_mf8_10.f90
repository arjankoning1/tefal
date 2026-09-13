subroutine allocate_mf8_10
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate temporary arrays for MF8-10
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: nen10
!
! Maximum number of MF9/10 energy points.
!
! E10 may contain the TALYS grid plus:
! - one extra point at the high-energy cutoff
! - the shifted cutoff point
! - EMAX
!
  nen10 = numinc + 3
!
! MF9/10 interpolation information
!
  allocate(INTER10(nummt,numiso,numint))
  allocate(NBT10(nummt,numiso,numint))
  allocate(NP10(nummt,numiso))
  allocate(NR10(nummt,numiso))

  allocate(INTERZA(numsec,numint))
  allocate(NBTZA(numsec,numint))
  allocate(NPZA(numsec))
  allocate(NRZA(numsec))
!
! MF9/10 cross sections
!
  allocate(E10(nummt,numiso,nen10))
  allocate(xsiso(nummt,numiso,nen10))

  allocate(E10ZA(numsec,nen10))
  allocate(xsrpZA(numsec,nen10))
!
! Initialization
!
  INTER10 = 0
  NBT10 = 0
  NP10 = 0
  NR10 = 0

  INTERZA = 0
  NBTZA = 0
  NPZA = 0
  NRZA = 0

  E10 = 0.
  xsiso = 0.
  E10ZA = 0.
  xsrpZA = 0.
!
  return
end subroutine allocate_mf8_10
