subroutine make4discrete(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF4 for discrete levels
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
!   sgl          ! single precision kind
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   Eahigh       ! upper energy of MT values to be adopted
!   Ealow        ! lower energy of MT values to be adopted
! Variables for input of ENDF library type
!   flaghigh     ! flag for high energies ( > 20 MeV)
! Variables for reaction initialization
!   Eangindex    ! enegy index for angular distributions
!   eninccut     ! last incident energy before high - energy format
!   Nenang       ! number of incident energies for ang. dist.
!   nlevmax      ! number of included discrete levels
!   numcut4      ! number of energies before MF4 high - energy format
!   rmu          ! cosine of the angle
! Variables for info from TALYS
!   eninc        ! incident energy
!   k0           ! index of incident particle
!   numinc       ! number of incident energies
! Variables for initialization of ENDF format
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
! Variables for discrete state cross sections in ENDF format
!   cleg0        ! Legendre coefficients
!   ncleg        ! number of Legendre coefficients
!   xsdisc       ! discrete state cross section
! Variables for angular distributions in ENDF format
!   fang         ! scattering angular distribution
! Variables for EndF format
!   INTER        ! interpolation scheme
!   NBT          ! separation value for interpolation scheme
!   NR           ! number of interpolation ranges
! Variables for MF1
!   EMAX         ! upper limit of energy range for evaluation
! Variables for MF3
!   E3           ! incident energy for MF3 (in ENDF - 6 format)
!   EthMT        ! threshold energy
! Variables for MF4
!   E4           ! incident energy for MF4 (in ENDF - 6 format)
!   E4h          ! incident energy for MF4 (in ENDF - 6 format)
!   E4hr         ! incident energy for MF4 (in ENDF - 6 format)
!   E4r          ! incident energy for MF4 (in ENDF - 6 format)
!   f4           ! angular distribution
!   f4r          ! angular distribution
!   INTER4       ! interpolation scheme
!   INTERh       ! interpolation scheme
!   LCT          ! LAB / CM flag
!   leg          ! Legendre coefficients (in ENDF - 6 format)
!   legr         ! Legendre coefficients (in ENDF - 6 format)
!   LI4          ! isotropy flag
!   LTT          ! representation
!   LVT          ! specification of transformation matrix
!   NBT4         ! separation value for interpolation scheme
!   NBTh         ! separation value for interpolation scheme
!   NE           ! number of incident energies
!   NE4r         ! number of incident energies (MF4 only)
!   NEh          ! number of incident energies (MF4 only)
!   NEhr         ! number of incident energies (MF4 only)
!   NL           ! Legendre order or number of cosines
!   NL4r         ! number of Legendre coefficients
!   NP4          ! number of incident energies
!   NP4r         ! number of angles (MF4 only)
!   NR4          ! number of interpolation ranges
!   NRh          ! number of interpolation ranges
!   x4           ! cosine of the angle
!   x4r          ! cosine of the angle
!
! *** Declaration of local data
!
  implicit none
  logical   :: adoptread            ! logical for adopting data
  logical   :: goahead              ! flag to determine existence of reaction
  integer   :: i                    ! counter
  integer   :: ia                   ! counter for alpha particles
  integer   :: iang                 ! running variable for angle
  integer   :: iE                   ! energy counter
  integer   :: iEh                  ! counter for energy
  integer   :: L                    ! counter for Legendre coefficients
  integer   :: MF                   ! MF-number
  integer   :: MT                   ! MT-number
  integer   :: nen                  ! energy counter
  integer   :: nencut4              ! energy point at cut off
  integer   :: nex                  ! discrete level
  integer   :: nin                  ! counter for incident energy
  integer   :: type                 ! particle type
  real(sgl) :: Eev                  ! energy in eV
  real(sgl) :: Ein                  ! incident energy
!
! **************************** Make MF4 ********************************
!
  MF = 4
  goahead = .false.
  if (MT >= 51 .and. MT <= 50 + nlevmax) then
    nex = MT - 50
    type = k0
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
  if (adopt(MF, MT)) then
    iE = NE
  else
    E4(1) = E3(MT, 1)
    NL(4, MT, 1) = 2
    leg(1, 1) = 0.
    leg(1, 2) = 0.
    iE = 1
  endif
  nencut4 = Nenang
  if (adopt(MF, MT)) then
    adoptread = .true.
  else
    adoptread = .false.
  endif
  do nen = 1, Nenang
    nin = Eangindex(nen)
    if (nin > numcut4) then
      nencut4 = nen - 1
      exit
    endif
    Ein = eninc(nin)
    if (Ein <= EthMT(MT)) cycle
    if (xsdisc(type, nex, nin) == 0.) cycle
    Eev = Ein * 1.e6
    if (adoptread .and. Eev >= Ealow(MF, MT) .and. Eev <= Eahigh(MF, MT)) then
      do i = 1, NE4r
        if (E4r(i) > Eahigh(MF, MT)) exit
        iE = iE + 1
        E4(iE) = E4r(i)
        NL(4, MT, iE) = NL4r(i)
        leg(iE, 0) = cleg0(type, nex, nen, 0)
        do L = 1, NL(4, MT, i)
          leg(iE, L) = legr(i, L)
        enddo
      enddo
      adoptread = .false.
    else
      if (adopt(MF, MT) .and. Eev >= Ealow(MF, MT) .and. Eev <= Eahigh(MF, MT)) cycle
      iE = iE + 1
      E4(iE) = Ein * 1.e6
      NL(4, MT, iE) = ncleg(type, nex, nen)
      leg(iE, 1) = 0.
      leg(iE, 2) = 0.
      do L = 0, NL(4, MT, iE)
        leg(iE, L) = cleg0(type, nex, nen, L)
      enddo
    endif
  enddo
  if (xsdisc(type, nex, numcut4) == 0.) then
    iE = iE + 1
    E4(iE) = Ein * 1.e6
    NL(4, MT, iE) = NL(4, MT, iE - 1)
    leg(iE, 1) = 0.
    leg(iE, 2) = 0.
    do L = 0, NL(4, MT, iE)
      leg(iE, L) = leg(iE - 1, L)
    enddo
  endif
  LTT = 1
!
! Tabulated angular distributions
!
  if (adopt(MF, MT)) then
    adoptread = .true.
  else
    adoptread = .false.
  endif
  if (numcut4 < numinc) then
    iEh = 0
    do nen = 1, Nenang
      nin = Eangindex(nen)
      if (nen < nencut4) cycle
      Ein = eninc(nin)
      if (Ein > eninccut) goto 200
      Eev = Ein * 1.e6
      if (adoptread .and. Eev >= Ealow(MF, MT) .and. Eev <= Eahigh(MF, MT)) then
        do i = 1, NEhr
          if (E4hr(i) > Eahigh(MF, MT)) exit
          iEh = iEh + 1
          E4(iEh) = E4hr(i)
          NP4(iEh) = NP4r(i)
          NR4(iEh) = 1
          NBT4(iEh, 1) = NP4r(i)
          INTER4(iEh, 1) = 2
          do ia = 1, NP4(iEh)
            x4(iEh, ia) = x4r(i, ia)
            f4(iEh, ia) = f4r(i, ia)
          enddo
        enddo
        adoptread = .false.
      else
        if (adopt(MF, MT) .and. Eev >= Ealow(MF, MT) .and. Eev <= Eahigh(MF, MT)) cycle
        iEh = iEh + 1
        E4h(iEh) = Ein * 1.e6
        ia = 0
        do iang = 90, 0, -1
          ia = ia + 1
          x4(iEh, ia) = rmu(iang)
          f4(iEh, ia) = fang(type, nex, nen, iang)
        enddo
        NP4(iEh) = ia
        NR4(iEh) = 1
        NBT4(iEh, 1) = ia
        INTER4(iEh, 1) = 2
      endif
    enddo
    NEh = iEh
    NRh = 1
    NBTh(1) = iEh
    INTERh(1) = 2
    LTT = 3
  endif
!
! High energies
!
200 if (flaghigh) then
    iE = iE + 1
    E4(iE) = EMAX
    NL(4, MT, iE) = NL(4, MT, iE - 1)
    do L = 0, NL(4, MT, iE)
      leg(iE, L) = leg(iE - 1, L)
    enddo
  endif
!
! 2. ENDF-6 parameters
!
  LVT = 0
  LI4 = 0
  LCT = 2
  NE = iE
  NR(4, MT) = 1
  NBT(4, MT, 1) = iE
  INTER(4, MT, 1) = 2
  mtexist(4, MT) = .true.
  mfexist(4) = .true.
  return
end subroutine make4discrete
! Copyright A.J. Koning 2021
