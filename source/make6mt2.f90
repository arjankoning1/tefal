subroutine make6mt2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF6 for MT2
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
! Constants
!   parA          ! mass number of particle
!   parspin       ! spin of particle
!   parZ          ! charge number of particle
! Variables for reaction initialization
!   Eangindex     ! enegy index for angular distributions
!   Nenang        ! number of incident energies for ang. dist.
!   relmass       ! mass relative to neutron mass
!   rmu           ! cosine of the angle
! Variables for info from TALYS
!   eninc         ! incident energy
!   k0            ! index of incident particle
!   numinc        ! number of incident energies
! Variables for initialization of ENDF format
!   AWR           ! standard mass parameter
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
!   ZA            ! standard charge parameter
! Variables for total cross sections in ENDF format
!   xselas        ! total elastic cross section
! Variables for angular distributions in ENDF format
!   fcpang        ! scattering angular distribution for charged particles
!   limang        ! smallest angle for charged - particle elastic scattering
! Variables for ENDF format
!   NK            ! number of subsections
! Variables for MF1
!   EMAX          ! upper limit of energy range for evaluation
! Variables for MF4
!   LCT           ! LAB / CM flag
!   NL            ! Legendre order or number of cosines
! Variables for MF6
!   AWP           ! product mass
!   b6            ! energy - angle values
!   E6            ! incident energy (in ENDF - 6 format) for distribution
!   Ey            ! incident energy for yields (in ENDF - 6 format)
!   INTER6ea      ! interpolation scheme
!   INTER6y       ! interpolation scheme
!   LAW           ! flag for distribution function
!   LIDP          ! indicates that particles are identical
!   LIP           ! product modifier flag
!   LTP           ! representation
!   NBT6ea        ! separation value for interpolation scheme
!   NBT6y         ! separation value for interpolation scheme
!   NE6ea         ! number of incident energies for distribution
!   NP6y          ! number of incident energies for yields
!   NR6ea         ! number of interpolation ranges for distribution
!   NR6y          ! number of interpolation ranges for yields
!   NW            ! number of words
!   SPIpar        ! particle spin
!   Y             ! product yield (in ENDF - 6 format)
!   ZAP           ! product identifier
!
! *** Declaration of local data
!
  implicit none
  integer :: iang              ! running variable for angle
  integer :: ib                ! counter
  integer :: iE                ! energy counter
  integer :: k                 ! counter
  integer :: MT                ! MT-number
  integer :: nen               ! energy counter
  integer :: nin               ! counter for incident energy
!
! ***** Charged-particle elastic scattering angular distributions ******
!
! LAW    : flag for distribution function
!
  MT = 2
  mtexist(6, MT) = .true.
  do k = 1, 2
    if (k == 1) then
      ZAP(k) = 1000. * parZ(k0) + parA(k0)
      AWP(k) = real(relmass(k0))
    else
      ZAP(k) = ZA
      AWP(k) = AWR
    endif
    LIP(k) = 0
!
! 1. Angular distributions
!
    iE = 0
    if (k == 1) then
      LAW(k) = 5
      SPIpar = parspin(k0)
      LIDP = 0
      do nen = 1, Nenang
        nin = Eangindex(nen)
        if (nin < numinc .and. xselas(nin) == 0.) cycle
        iE = iE + 1
        E6(k, iE) = eninc(nin) * 1.e6
        ib = 0
        do iang = 90, limang, -1
          ib = ib + 1
          b6(k, iE, ib) = rmu(iang)
          ib = ib + 1
          b6(k, iE, ib) = fcpang(nen, iang)
        enddo
        NL(6, MT, iE) = 91 - limang
        NW(k, iE) = 2 * NL(6, MT, iE)
      enddo
!
! Recoils
!
    else
      LAW(k) = 4
    endif
!
! 2. Yields
!
    NP6y(k) = 2
    NR6y(k) = 1
    NBT6y(k, 1) = NP6y(k)
    INTER6y(k, 1) = 2
    Ey(k, 1) = E6(1, 1)
    Y(k, 1) = 1.
    Ey(k, 2) = EMAX
    Y(k, 2) = 1.
!
! ENDF parameters
!
! write6  : subroutine to write MF6
!
    NE6ea(k) = iE
    NR6ea(k) = 1
    NBT6ea(k, 1) = NE6ea(k)
    INTER6ea(k, 1) = 2
  enddo
  LTP = 12
  LCT = 2
  NK(6, MT) = 2
  mfexist(6) = .true.
  call write6(MT)
  return
end subroutine make6mt2
! Copyright A.J. Koning 2021
