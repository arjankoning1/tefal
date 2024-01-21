subroutine make6discrete(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF6 for discrete levels
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
! Definition of single and double precision variables
!   sgl           ! single precision kind
! Variables for input of ENDF library type
!   flaghigh      ! flag for high energies ( > 20 MeV)
! Constants
!   parA          ! mass number of particle
!   parmass       ! mass of particle in a.m.u.
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Variables for reaction initialization
!   Eangindex     ! enegy index for angular distributions
!   nlevmax       ! number of included discrete levels
!   Nenang        ! number of incident energies for ang. dist.
!   numcut        ! number of energies before high - energy format
!   numcut4       ! number of energies before MF4 high - energy format
!   relmass       ! mass relative to neutron mass
! Variables for info from TALYS
!   Ainit         ! mass number of initial compound nucleus
!   eninc         ! incident energy
!   k0            ! index of incident particle
!   Zinit         ! charge number of initial compound nucleus
! Variables for initialization of ENDF format
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
! Variables for residual production cross sections in ENDF format
!   nucmass       ! mass of nucleus
! Variables for discrete state cross sections in ENDF format
!   cleg0         ! Legendre coefficients
!   ncleg         ! number of Legendre coefficients
!   xsdisc        ! discrete state cross section
! Variables for photon production in ENDF format
!   Egamdis       ! energy of gamma ray
!   Ngamdis       ! number of gamma ray lines per level
!   yieldg        ! total discrete gamma yield per level
!   yieldratio    ! yield ratio for level
! Variables for ENDF format
!   NK            ! number of subsections
! Variables for MF1
!   EMAX          ! upper limit of energy range for evaluation
! Variables for MF3
!   E3            ! incident energy for MF3 (in ENDF - 6 format)
!   EthMT         ! threshold energy
! Variables for MF4
!   LCT           ! LAB / CM flag
!   leg           ! Legendre coefficients (in ENDF - 6 format)
!   NL            ! Legendre order or number of cosines
! Variables for MF6
!   AWP           ! product mass
!   b6            ! energy - angle values
!   E6            ! incident energy (in ENDF - 6 format) for distribution
!   Ey            ! incident energy for yields (in ENDF - 6 format)
!   INTER6ea      ! interpolation scheme
!   INTER6y       ! interpolation scheme
!   LANG          ! flag for angular representation
!   LAW           ! flag for distribution function
!   LEP           ! interpolation scheme for secondary energy
!   LIP           ! product modifier flag
!   NA            ! number of angular parameters
!   NBT6ea        ! separation value for interpolation scheme
!   NBT6y         ! separation value for interpolation scheme
!   ND            ! number of discrete energies
!   NE6ea         ! number of incident energies for distribution
!   NEP           ! number of secondary energy points
!   NP6y          ! number of incident energies for yields
!   NR6ea         ! number of interpolation ranges for distribution
!   NR6y          ! number of interpolation ranges for yields
!   NW            ! number of words
!   Y             ! product yield (in ENDF - 6 format)
!   ZAP           ! product identifier
!
! *** Declaration of local data
!
  implicit none
  logical   :: goahead            ! flag to determine existence of reaction
  integer   :: A                  ! mass number of target nucleus
  integer   :: i                  ! counter
  integer   :: ib                 ! counter
  integer   :: iE                 ! energy counter
  integer   :: k                  ! counter
  integer   :: L                  ! counter for Legendre coefficients
  integer   :: MT                 ! MT-number
  integer   :: nen                ! energy counter
  integer   :: nex                ! discrete level
  integer   :: nin                ! counter for incident energy
  integer   :: Nix                ! neutron number index for residual nucleus
  integer   :: type               ! particle type
  integer   :: Z                  ! charge number of target nucleus
  integer   :: Zix                ! charge number index for residual nucleus
  real(sgl) :: Ein                ! incident energy
!
! **************************** Make MF6 ********************************
!
! LAW     : flag for distribution function
!
  goahead = .false.
  if (MT >= 50 .and. MT <= 50 + nlevmax) then
    nex = MT - 50
    type = 1
    goahead = .true.
  endif
  if (MT >= 600 .and. MT <= 600 + nlevmax) then
    nex = MT - 600
    type = 2
    goahead = .true.
  endif
  if (MT >= 650 .and. MT <= 650 + nlevmax) then
    nex = MT - 650
    type = 3
    goahead = .true.
  endif
  if (MT >= 700 .and. MT <= 700 + nlevmax) then
    nex = MT - 700
    type = 4
    goahead = .true.
  endif
  if (MT >= 750 .and. MT <= 750 + nlevmax) then
    nex = MT - 750
    type = 5
    goahead = .true.
  endif
  if (MT >= 800 .and. MT <= 800 + nlevmax) then
    nex = MT - 800
    type = 6
    goahead = .true.
  endif
  if ( .not. goahead) return
!
! Loop over subsections
!
  if (mod(MT, 50) == 0 .or. yieldg(type, nex) == 0.) then
    NK(6, MT) = 2
  else
    NK(6, MT) = 3
  endif
  do k = 1, NK(6, MT)
    iE = 0
    if (k == 1) then
      ZAP(k) = 1000. * parZ(type) + parA(type)
      AWP(k) = real(relmass(type))
      Y(k, 1) = 1.
      Y(k, 2) = 1.
    endif
    if (k == 2) then
      Zix = parZ(type)
      Nix = parN(type)
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      ZAP(k) = 1000. * Z + A
      AWP(k) = real(nucmass(Zix, Nix) / parmass(1))
      Y(k, 1) = 1.
      Y(k, 2) = 1.
    endif
    if (k == 3) then
      ZAP(k) = 0.
      AWP(k) = 0.
      Y(k, 1) = yieldg(type, nex)
      Y(k, 2) = yieldg(type, nex)
    endif
    Ey(k, 1) = E3(MT, 1)
    Ey(k, 2) = EMAX
    NP6y(k) = 2
    NR6y(k) = 1
    NBT6y(k, 1) = NP6y(k)
    INTER6y(k, 1) = 2
!
! k=1: Legendre coefficients
!
    if (k == 1) then
      LAW(k) = 2
      LIP(k) = 0
      LANG(k) = 0
      E6(k, 1) = E3(MT, 1)
      NL(6, MT, 1) = 2
      NW(k, 1) = NL(6, MT, 1)
      leg(1, 1) = 0.
      leg(1, 2) = 0.
!
! To avoid very big files, we do not give angular information at each incident energy, but instead skip various energy points.
! The number of skipped points may be different for neutrons and charged particles.
!
      iE = 1
      do nen = 1, Nenang
        nin = Eangindex(nen)
        if (nin > numcut) exit
        Ein = eninc(nin)
        if (eninc(nin) <= EthMT(MT)) cycle
        if (xsdisc(type, nex, nin) == 0.) cycle
        iE = iE + 1
        E6(k, iE) = Ein * 1.e6
        if (k0 == 0) then
          NL(6, MT, iE) = 2
        else
          NL(6, MT, iE) = ncleg(type, nex, nen)
        endif
        NW(k, iE) = NL(6, MT, iE)
        leg(iE, 1) = 0.
        leg(iE, 2) = 0.
        if (nin <= numcut4) then
          do L = 0, NL(6, MT, iE)
            leg(iE, L) = cleg0(type, nex, nen, L)
          enddo
        else
          NL(6, MT, iE) = NL(6, MT, iE - 1)
          NW(k, iE) = NW(k, iE - 1)
          do L = 0, NL(6, MT, iE)
            leg(iE, L) = leg(iE - 1, L)
          enddo
        endif
      enddo
      if (xsdisc(type, nex, numcut) == 0.) then
        iE = iE + 1
        E6(k, iE) = Ein * 1.e6
        NL(6, MT, iE) = NL(6, MT, iE - 1)
        NW(k, iE) = NW(k, iE - 1)
        leg(iE, 1) = 0.
        leg(iE, 2) = 0.
        do L = 0, NL(6, MT, iE)
          leg(iE, L) = leg(iE - 1, L)
        enddo
      endif
!
! High energies
!
      if (flaghigh) then
        iE = iE + 1
        E6(k, iE) = EMAX
        NL(6, MT, iE) = NL(6, MT, iE - 1)
        NW(k, iE) = NW(k, iE - 1)
        do L = 0, NL(6, MT, iE)
          leg(iE, L) = leg(iE - 1, L)
        enddo
      endif
    endif
!
! k=2: recoils
!
    if (k == 2) then
      LAW(k) = 4
      LIP(k) = 0
    endif
!
! k=3: photons
!
    if (k == 3) then
      LAW(k) = 1
      LIP(k) = 0
      LANG(k) = 1
      LEP(k) = 2
      E6(k, 1) = E3(MT, 1)
      E6(k, 2) = EMAX
      ND(k, 1) = Ngamdis(type, nex)
      ND(k, 2) = ND(k, 1)
      NA(k, 1) = 0
      NA(k, 2) = 0
      NEP(k, 1) = ND(k, 1)
      NEP(k, 2) = NEP(k, 1)
      NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
      NW(k, 2) = NEP(k, 2) * (NA(k, 2) + 2)
      ib = 0
      do i = 1, ND(k, 1)
        ib = ib + 1
        b6(k, 1, ib) = Egamdis(type, nex, i) * 1.e6
        b6(k, 2, ib) = b6(k, 1, ib)
        ib = ib + 1
        b6(k, 1, ib) = yieldratio(type, nex, i)
        b6(k, 2, ib) = b6(k, 1, ib)
      enddo
      iE = 2
    endif
!
! 2. ENDF-6 parameters
!
! write6  : subroutine to write MF6
!
    NE6ea(k) = iE
    NR6ea(k) = 1
    NBT6ea(k, 1) = NE6ea(k)
    INTER6ea(k, 1) = 2
  enddo
  LCT = 2
  mtexist(6, MT) = .true.
  mfexist(6) = .true.
  call write6(MT)
  return
end subroutine make6discrete
! Copyright A.J. Koning 2021
