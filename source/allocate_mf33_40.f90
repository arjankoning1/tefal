subroutine allocate_mf33_40
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF33 and MF40
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: ncovtot
  integer :: nread
!
! Maximum number of covariance matrix entries.
!
  ncovtot = 1 + Nencov * Nencov
!
! External ENDF data are read in blocks of 6 values.
! Round the read arrays upward to a multiple of 6.
!
  nread = 6 * ((ncovtot + 5) / 6)
!
! Arrays for adopted covariance data.
!
  allocate(b33read(Nchancov,Nchancov,nread))
  allocate(b33MTread(Nchancov,nread))
  allocate(b8read(Nchancov,nread))
!
! Arrays for generated MF33/MF40 data.
!
  allocate(b33(Nchancov,Nchancov,ncovtot))
  allocate(b33MT(Nchancov,ncovtot))
  allocate(b33ZA(Ncovrp,ncovtot))
  allocate(b8(Nchancov,ncovtot))
!
! Initialization.
!
  b33read = 0.
  b33MTread = 0.
  b8read = 0.

  b33 = 0.
  b33MT = 0.
  b33ZA = 0.
  b8 = 0.

  return
end subroutine allocate_mf33_40
