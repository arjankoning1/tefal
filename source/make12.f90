subroutine make12
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF12
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
!   sgl            ! single precision kind
! All global variables
!   nummt          ! number of MT numbers
! Variables for input of ENDF structure
!   flagcapt6      ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagdisc6      ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
!   flaggam13      ! flag to use MF13 for gamma prod. instead of MF12 (if not in MF6)
!   flaggamdisc    ! flag to store gamma prod. per level (MT51..) instead of MT4
!   flagpart6      ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
! Variables for input of specific ENDF data
!   adopt          ! logical for existence of MF information (per MT)
! Variables for input of ENDF library type
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps         ! energy shift at MT5 cutoff energy (in eV)
! Variables for partial cross sections in ENDF format
!   idchannel      ! identifier for channel
! Variables for reaction initialization
!   Egamindex      ! enegy index for gamma cross sections
!   Nengam         ! number of incident energies for gamna c.s.
!   numcut         ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for initialization of ENDF format
!   blank2         ! blank string
!   FEND           ! ENDF - 6 format
!   idnum          ! number of different exclusive cross sections
!   MAT            ! MAT number
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
!   MTid           ! channel identifier for MT - number
! Variables for discrete state cross sections in ENDF format
!   Qdisc          ! Q - value
! Variables for photon production in ENDF format
!   branchlevel    ! level to which branching takes place
!   branchratio    ! branch ratio
!   Egamma         ! gamma energy
!   Nbranch        ! number of branches for level
!   Ngam           ! number of gamma transitions for this nucleus
!   yieldcont      ! continuum gamma - ray yield
!   yielddisc      ! discrete gamma - ray yield
!   yieldtot       ! total gamma - ray yield
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NK             ! number of subsections
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF1
!   EMAX           ! upper limit of energy range for evaluation
! Variables for MF3
!   EthMT          ! threshold energy
!   QI             ! Q - value (in ENDF - 6 format)
! Variables for MF12_15
!   E12            ! incident energy (in ENDF - 6 format)
!   Eg             ! gamma energy
!   Egk            ! gamma energy
!   ES12           ! energy of level
!   Esk            ! starting level (in ENDF - 6 format)
!   ESNS           ! energy of mother level
!   INTERg         ! interpolation scheme
!   LFg            ! photo energy distribution law
!   LG12           ! type setters
!   LO12           ! type setters
!   LP12           ! origin of photons
!   LPg            ! primary photon flag
!   NBTg           ! separation value for interpolation scheme
!   NPg            ! number of incident energies
!   NS12           ! number of levels below the present one
!   NT12           ! number of transitions for which data are given
!   TP12           ! probability of direct transition
!   xsgtotyield    ! total discrete photon multiplicity
!   xsgyield       ! gamma - ray multiplicity (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagpartmt            ! flag for partial MT numbers
  integer   :: i                     ! counter
  integer   :: id                    ! counter for deuterons
  integer   :: idc                   ! help variable
  integer   :: iE                    ! energy counter
  integer   :: iEg                   ! index for gamma energy
  integer   :: igam                  ! counter for gammas
  integer   :: ii                    ! counter
  integer   :: MF                    ! MF-number
  integer   :: MT                    ! MT-number
  integer   :: nen                   ! energy counter
  integer   :: nencut                ! number of energies before high-energy format
  integer   :: nex                   ! discrete level
  integer   :: nex0                  ! base number for discrete level
  integer   :: nin                   ! counter for incident energy
  integer   :: Nphoton               ! total number of photons (NJOY restriction)
  integer   :: NS                    ! line number
  integer   :: type                  ! particle type
  real(sgl) :: Ein                   ! incident energy
  real(sgl) :: Exd                   ! excitation energy
  real(sgl) :: x                     ! help variable
!
! ***************************** Make MF12 ******************************
!
! read12   : subroutine to read MF12 from existing ENDF-6 data library
!
  MF = 12
  NS = 0
  Nphoton = 0
  open (unit = 2, file = 'MF12', status = 'replace')
  do MT = 4, nummt
    if ( .not. mtexist(3, MT)) cycle
    if (adopt(13, MT)) cycle
    if (MT == 18 .or. MT == 19 .or. MT == 20 .or. MT == 21 .or. MT == 38) cycle
    if (MT == 600 .or. MT == 650 .or. MT == 700 .or. MT == 750 .or. MT == 800) cycle
    if (MT == 649 .or. MT == 699 .or. MT == 749 .or. MT == 799 .or. MT == 849) cycle
    flagpartmt = (MT <= 50 .or. (MT > 90 .and. MT < 600))
    if (flaggamdisc) then
      if (MT == 4 .or. (MT >= 103 .and. MT <= 107)) cycle
    else
      if ((MT >= 51 .and. MT <= 91) .or. MT >= 600) cycle
    endif
    if (MT == 102 .and. flagcapt6) cycle
    if (flagpartmt .and. flagpart6 .and. MT /= 102) cycle
    if (flagpartmt .and. flaggam13) cycle
    if ( .not. flagpartmt .and. flagdisc6) cycle
    iE = 0
    if (adopt(MF, MT)) then
      call read12(MT)
      goto 100
    endif
    do id = 0, idnum
      if (idchannel(id) == MTid(MT)) then
        idc = id
        exit
      endif
    enddo
!
! 1. Total gamma-ray yields
!
    if (flagpartmt) then
      nencut = Nengam
      do nen = 1, Nengam
        nin = Egamindex(nen)
        if (nin > numcut) then
          nencut = nen - 1
          exit
        endif
        Ein = eninc(nin)
        x = yieldtot(idc, nen)
        if (nin > 1 .and. x == 0.) cycle
        iE = iE + 1
        E12(idc, iE) = Ein * 1.e6
        xsgtotyield(idc, iE) = x
      enddo
      if (iE == 1) cycle
      x = yieldtot(idc, nencut)
      if (x == 0.) then
        iE = iE + 1
        E12(idc, iE) = Ein * 1.e6
        xsgtotyield(idc, iE) = 0.
      endif
!
! High energies
!
      if (flaghigh) then
        iE = iE + 1
        E12(idc, iE) = Ein * 1.e6 + cuteps
        xsgtotyield(idc, iE) = x
        iE = iE + 1
        E12(idc, iE) = EMAX
        xsgtotyield(idc, iE) = x
      endif
!
! 2. Gamma-ray yields per discrete transition
!
      igam = 0
      do i = 1, Ngam(idc) + 1
        if (Egamma(idc, i) >= 20.) cycle
        igam = igam + 1
        if (i == Ngam(idc) + 1) then
          Esk(idc, igam) = 0.
          Egk(idc, igam) = 0.
        else
          Esk(idc, igam) = 0.
          Egk(idc, igam) = Egamma(idc, i) * 1.e6
        endif
        Eg(idc, igam, 1) = EthMT(MT) * 1.e6
        iEg = 0
        do nen = 1, Nengam
          nin = Egamindex(nen)
          if (nin > numcut) exit
          Ein = eninc(nin)
          if (i == Ngam(idc) + 1) then
            x = yieldcont(idc, nen)
          else
            x = yielddisc(idc, i, nen)
          endif
          if (nin > 1 .and. x == 0.) cycle
          iEg = iEg + 1
          Eg(idc, igam, iEg) = Ein * 1.e6
          xsgyield(idc, igam, iEg) = x
        enddo
        if (iEg == 1) then
          igam = igam - 1
          cycle
        endif
        x = 0.
        if (i == Ngam(idc) + 1) then
          x = yieldcont(idc, nencut)
          LFg(MT, igam) = 1
        else
          x = yielddisc(idc, i, nencut)
          LFg(MT, igam) = 2
        endif
        if (x == 0.) then
          iEg = iEg + 1
          Eg(idc, igam, iEg) = Ein * 1.e6
          xsgyield(idc, igam, iEg) = 0.
        endif
!
! High energies
!
        if (flaghigh) then
          iEg = iEg + 1
          Eg(idc, igam, iEg) = Ein * 1.e6 + cuteps
          xsgyield(idc, igam, iEg) = 0.
          iEg = iEg + 1
          Eg(idc, igam, iEg) = EMAX
          xsgyield(idc, igam, iEg) = 0.
        endif
        NPg(MT, igam) = iEg
        NBTg(MT, igam, 1) = iEg
        INTERg(MT, igam, 1) = 2
        LPg(MT, igam) = 0
      enddo
      if (igam == 1) Eg(idc, igam, 1) = E12(idc, 1)
      NK(MF, MT) = igam
!
! ENDF-6 parameters
!
      LO12(MT) = 1
      LG12(MT) = 0
      NP(MF, MT) = iE
      NR(MF, MT) = 1
      NBT(MF, MT, 1) = iE
      INTER(MF, MT, 1) = 2
    else
!
! Transition probability arrays
!
! write12: subroutine to write MF12
!
      if (MT <= 90) then
        nex = MT - 50
        type = 1
      else
        nex = mod(MT, 50)
        type = (MT + 1 - 500) / 50
      endif
      LO12(MT) = 2
      LG12(MT) = 1
      LP12(MT) = 0
      NS12(MT) = nex
      if (type == k0 .and. nex <= Ltarget) then
        Exd = real(dble((Qdisc(type, 0) - Qdisc(type, nex - 1)) * 1.e6))
        nex0 = max(nex - 1, 1)
      else
        if (MT <= 90) then
          Exd = real(dble(abs(QI(MT))))
          if (Exd < 1.e4) Exd = Exd + 0.01
        else
          Exd = real(dble((Qdisc(type, 0) - Qdisc(type, nex)) * 1.e6))
        endif
        nex0 = nex
      endif
      ESNS(MT) = real(int(Exd))
      Nphoton = Nphoton + Nbranch(type, nex0)
      if (Nphoton > 100) cycle
      do i = 1, Nbranch(type, nex0)
        ii = branchlevel(type, nex0, i)
        Exd = real(dble((Qdisc(type, 0) - Qdisc(type, ii)) * 1.e6))
        ES12(MT, i) = real(int(Exd))
        TP12(MT, i) = 0.01 * branchratio(type, nex0, i)
      enddo
      NT12(MT) = Nbranch(type, nex0)
      NK(MF, MT) = NT12(MT)
    endif
  100   call write12(MT)
    mtexist(MF, MT) = .true.
    mfexist(MF) = .true.
  enddo
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine make12
! Copyright A.J. Koning 2021
