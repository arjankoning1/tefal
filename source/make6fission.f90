subroutine make6fission
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF6 for partial cross sections
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! Variables for initialization of ENDF format
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
! Variables for ENDF format
!   NK            ! number of subsections
! Variables for MF1
!   EMAX          ! upper limit of energy range for evaluation
! Variables for MF3
!   E3            ! incident energy for MF3 (in ENDF - 6 format)
! Variables for MF4
!   LCT           ! LAB / CM flag
! Variables for MF6
!   AWP           ! product mass
!   Ey            ! incident energy for yields (in ENDF - 6 format)
!   INTER6y       ! interpolation scheme
!   LAW           ! flag for distribution function
!   LIP           ! product modifier flag
!   NBT6y         ! separation value for interpolation scheme
!   NP6y          ! number of incident energies for yields
!   NR6y          ! number of interpolation ranges for yields
!   Y             ! product yield (in ENDF - 6 format)
!   ZAP           ! product identifier
!
! *** Declaration of local data
!
  implicit none
  integer :: MT              ! MT-number
!
! ***************** Make MF6 for fission cross sections ****************
!
! write6 : subroutine to write MF6
!
  MT = 18
  mfexist(6) = .true.
  mtexist(6, MT) = .true.
  ZAP(1) = 1.
  AWP(1) = 1.
  LAW(1) = 3
  NK(6, MT) = 1
  LCT = 2
  LIP(1) = 0
  NP6y(1) = 2
  NR6y(1) = 1
  NBT6y(1, 1) = NP6y(1)
  INTER6y(1, 1) = 2
  Ey(1, 1) = E3(MT, 1)
  Ey(1, 2) = EMAX
  Y(1, 1) = 1.
  Y(1, 2) = 1.
  call write6(MT)
  return
end subroutine make6fission
! Copyright A.J. Koning 2021
