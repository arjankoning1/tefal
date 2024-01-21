subroutine make13
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF13
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
!   sgl             ! single precision kind
! All global variables
!   nummt           ! number of MT numbers
! Variables for input of ENDF structure
!   flagcapt6       ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagpart6       ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
! Variables for input of specific ENDF data
!   adopt           ! logical for existence of MF information (per MT)
! Variables for input of ENDF library type
!   flaghigh        ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps          ! energy shift at MT5 cutoff energy (in eV)
! Variables for reaction initialization
!   Egamindex       ! enegy index for gamma cross sections
!   Nengam          ! number of incident energies for gamna c.s.
!   numcut          ! number of energies before high - energy format
! Variables for initialization of ENDF format
!   blank2          ! blank string
!   FEND            ! ENDF - 6 format
!   idnum           ! number of different exclusive cross sections
!   MAT             ! MAT number
!   mfexist         ! flag for existence of MF - number
!   mtexist         ! flag for existence of MT - number
!   MTid            ! channel identifier for MT - number
! Variables for info from TALYS
!   eninc           ! incident energy
!  Variables for partial cross sections in ENDF format
!   idchannel       ! identifier for channel
! Variables for photon production in ENDF format
!   Egamma          ! gamma energy
!   Estart          ! starting level
!   Ngam            ! number of gamma transitions for this nucleus
!   xsgam           ! exclusive discrete gamma - ray cross section
!   xsgamdisctot    ! total discrete gamma - ray cross section
! Variables for ENDF format
!   INTER           ! interpolation scheme
!   NK              ! number of subsections
!   NP              ! number of incident energies
!   NBT             ! separation value for interpolation scheme
!   NR              ! number of interpolation ranges
! Variables for MF1
!   EMAX            ! upper limit of energy range for evaluation
! Variables for MF3
!   EthMT           ! threshold energy
! Variables for MF12_15
!   E13             ! incident energy (in ENDF - 6 format)
!   Eg              ! gamma energy
!   Egk             ! gamma energy
!   Esk             ! starting level (in ENDF - 6 format)
!   INTERg          ! interpolation scheme
!   LFg             ! photo energy distribution law
!   LPg             ! primary photon flag
!   NBTg            ! separation value for interpolation scheme
!   NPg             ! number of incident energies
!   xsg             ! gamma - ray cross section (in ENDF - 6 format)
!   xsgtot          ! total discrete photon production cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                 ! counter
  integer   :: id                ! counter for deuterons
  integer   :: idc               ! help variable
  integer   :: iE                ! energy counter
  integer   :: iEg               ! index for gamma energy
  integer   :: igam              ! counter for gammas
  integer   :: MF                ! MF-number
  integer   :: MT                ! MT-number
  integer   :: nen               ! energy counter
  integer   :: nencut            ! number of energies before high-energy format
  integer   :: nin               ! counter for incident energy
  integer   :: NS                ! line number
  real(sgl) :: Ein               ! incident energy
  real(sgl) :: x                 ! help variable
!
! ***************************** Make MF13 ******************************
!
  MF = 13
  NS = 0
  open (unit = 2, file = 'MF13', status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(3, MT)) cycle
    if (MTid(MT) <=  -1) cycle
    if (adopt(12, MT)) cycle
    if ((MT >= 51 .and. MT <= 91) .or. MT >= 600) cycle
    if (MT == 102) then
      if (flagcapt6) cycle
    else
      if (flagpart6) cycle
    endif
    do id = 0, idnum
      if (idchannel(id) == MTid(MT)) then
        idc = id
        exit
      endif
    enddo
!
! 1. Total discrete gamma-ray production cross sections
!
    E13(idc, 1) = EthMT(MT) * 1.e6
    xsgtot(idc, 1) = 0.
    iE = 1
    nencut = Nengam
    do nen = 1, Nengam
      nin = Egamindex(nen)
      if (nin > numcut) then
        nencut = nen - 1
        exit
      endif
      Ein = eninc(nin)
      x = xsgamdisctot(idc, nen)
      if (x == 0.) cycle
      iE = iE + 1
      E13(idc, iE) = Ein * 1.e6
      xsgtot(idc, iE) = x * 1.e-3
    enddo
    if (iE == 1) cycle
    x = xsgamdisctot(idc, nencut)
    if (x == 0.) then
      iE = iE + 1
      E13(idc, iE) = Ein * 1.e6
      xsgtot(idc, iE) = 0.
    endif
!
! High energies
!
    if (flaghigh) then
      iE = iE + 1
      E13(idc, iE) = Ein * 1.e6 + cuteps
      xsgtot(idc, iE) = 0.
      iE = iE + 1
      E13(idc, iE) = EMAX
      xsgtot(idc, iE) = 0.
    endif
!
! 2. Gamma-ray production cross sections per discrete transition
!
    igam = 0
    do i = 1, Ngam(idc)
      if (Egamma(idc, i) >= 20.) cycle
      igam = igam + 1
      Esk(idc, igam) = Estart(idc, i) * 1.e6
      Egk(idc, igam) = Egamma(idc, i) * 1.e6
      Eg(idc, igam, 1) = EthMT(MT) * 1.e6
      xsg(idc, igam, 1) = 0.
      iEg = 1
      do nen = 1, Nengam
        nin = Egamindex(nen)
        if (nin > numcut) exit
        Ein = eninc(nin)
        x = xsgam(idc, i, nen)
        if (x == 0.) cycle
        iEg = iEg + 1
        Eg(idc, igam, iEg) = Ein * 1.e6
        xsg(idc, igam, iEg) = x * 1.e-3
      enddo
      if (iEg == 1) then
        igam = igam - 1
        cycle
      endif
      x = xsgam(idc, i, numcut)
      if (x == 0.) then
        iEg = iEg + 1
        Eg(idc, igam, iEg) = Ein * 1.e6
        xsg(idc, igam, iEg) = 0.
      endif
!
! High energies
!
      if (flaghigh) then
        iEg = iEg + 1
        Eg(idc, igam, iEg) = Ein * 1.e6 + cuteps
        xsg(idc, igam, iEg) = 0.
        iEg = iEg + 1
        Eg(idc, igam, iEg) = EMAX
        xsg(idc, igam, iEg) = 0.
      endif
      NPg(MT, igam) = iEg
      NBTg(MT, igam, 1) = iEg
      INTERg(MT, igam, 1) = 2
      LPg(MT, igam) = 0
      LFg(MT, igam) = 2
    enddo
    NK(MF, MT) = igam
    if (igam == 1) Eg(idc, igam, 1) = E13(idc, 1)
!
! ENDF-6 parameters
!
! write13: subroutine to write MF13
!
    NP(MF, MT) = iE
    NR(MF, MT) = 1
    NBT(MF, MT, 1) = iE
    INTER(MF, MT, 1) = 2
    if (mtexist(3, MT) .and. NK(MF, MT) >= 1) then
      call write13(MT)
      mtexist(MF, MT) = .true.
      mfexist(MF) = .true.
    endif
  enddo
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine make13
! Copyright A.J. Koning 2021
