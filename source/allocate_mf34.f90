subroutine allocate_mf34
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF34
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: i
  integer :: nE34max
  integer :: nT34max
!
! Determine number of covariance-energy points retained by make34.
! make34 uses covstep=2 and retains the first and last points.
!
  nE34max = 0
  do i = 1, Nchanleg
    if (mod(i,2) == 1 .and. i /= 1 .and. i /= Nchanleg) cycle
    nE34max = nE34max + 1
  enddo
!
! MF34 contains the energy grid followed by the upper triangular
! covariance matrix. This is equal to N*(N+1)/2 entries.
!
  nT34max = nE34max * (nE34max + 1) / 2
!
! Allocate MF34 arrays.
!
  allocate(NE34(Nleg34,Nleg34))
  allocate(NI34(Nleg34,Nleg34))
  allocate(NT34(Nleg34,Nleg34))
  allocate(b34(Nleg34,Nleg34,nT34max))
!
! Initialization.
!
  NE34 = 0
  NI34 = 0
  NT34 = 0
  b34 = 0.

  return
end subroutine allocate_mf34
