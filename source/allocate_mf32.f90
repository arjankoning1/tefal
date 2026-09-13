subroutine allocate_mf32
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF32
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: ncov32
  integer :: nb32
!
! Maximum number of compact covariance records.
!
  ncov32 = numrescov * numrescov / 2
!
! Maximum size of conventional MF32 covariance array.
!
  nb32 = 6 * numrescov + &
 &       numrescov * numrespar * (numrescov * numrespar + 1)
!
! Compact covariance arrays.
!
  allocate(covdigit(ncov32,14))
  allocate(covix32(ncov32,2))
!
! MF32 control arrays.
!
  allocate(ISR(numres))
  allocate(LCOMP(numres))
  allocate(MLS(numres))
  allocate(NJS32(numres,numlres))
  allocate(NLS32(numres))
!
! MF32 resonance information.
!
  allocate(AJ32(numres,numlres,numjres))
  allocate(D32(numres,numlres,numjres))
  allocate(DAP(numres))
  allocate(GF32(numres,numlres,numjres))
  allocate(GG32(numres,numlres,numjres))
  allocate(GNO32(numres,numlres,numjres))
  allocate(GX32(numres,numlres,numjres))
!
! Covariance arrays.
!
  allocate(b32(nb32))
  allocate(b32URR(5*numjres*(numjres+1)))
!
! Initialization.
!
  covdigit = '    '
  covix32 = 0

  ISR = 0
  LCOMP = 0
  MLS = 0
  NJS32 = 0
  NLS32 = 0

  AJ32 = 0.
  D32 = 0.
  DAP = 0.
  GF32 = 0.
  GG32 = 0.
  GNO32 = 0.
  GX32 = 0.

  b32 = 0.
  b32URR = 0.

  return
end subroutine allocate_mf32
