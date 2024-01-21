subroutine make3discrete
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3 for discrete levels and continuum
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
! All global variables
!   numenin      ! number of incident energies
! Variables for input of specific ENDF data
!   urrcomp      ! mode for competition in the URR, 0: none, 1:MT4, 2:all
! Variables for input of ENDF library type
!   flaghigh     ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps       ! energy shift at MT5 cutoff energy (in eV)
!   Eswitch      ! energy where ENDF - 6 representation is switched (in MeV)
!   NMTmax       ! maximum number of MT numbers
! Variables for input of ENDF structure
!   flagpara     ! flag to include partial cross sections for alphs
!   flagpard     ! flag to include partial cross sections for deuterons
!   flagparh     ! flag to include partial cross sections for helions
!   flagparp     ! flag to include partial cross sections for protons
!   flagpart     ! flag to include partial cross sections for tritons
! Variables for reaction initialization
!   EHres        ! upper energy in resonance range
!   eninccut     ! last incident energy before high - energy format
!   nlevmax      ! number of included discrete levels
!   numcut       ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc        ! incident energy
!   k0           ! index of incident particle
!   Ltarget      ! excited level of target
!   numinc       ! number of incident energies
! Variables for initialization of ENDF format
!   AWR          ! standard mass parameter
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
! Variables for discrete state cross sections in ENDF format
!   edisc        ! incident energy for discrete level cross section
!   Ethdisc      ! threshold energy
!   numendisc    ! number of incident energies including discrete energy grid
!   Qdisc        ! Q - value
!   xscont       ! continuum cross section
!   xsdisc       ! discrete state cross section
!   xsintdisc    ! interpolated discrete cross section
! Variables for ENDF format
!   INTER        ! interpolation scheme
!   NBT          ! separation value for interpolation scheme
!   NP           ! number of incident energies
!   NR           ! number of interpolation ranges
! Variables for MF1
!   EMAX         ! upper limit of energy range for evaluation
! Variables for MF3
!   E3           ! incident energy for MF3 (in ENDF - 6 format)
!   EthMT        ! threshold energy
!   LR3          ! breakup flag
!   NMT          ! total number of MT sections
!   QI           ! Q - value (in ENDF - 6 format)
!   QM           ! Q - value (in ENDF - 6 format)
!   xs           ! cross section
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagU                        ! flag for URR
  integer   :: i                            ! counter
  integer   :: iE                           ! energy counter
  integer   :: j                            ! counter
  integer   :: MF                           ! MF-number
  integer   :: MT                           ! MT-number
  integer   :: MT1                          ! MT-number for first discrete level
  integer   :: MTc                          ! MT-number for continuum
  integer   :: MTlast                       ! MT number for last discrete level
  integer   :: nen                          ! energy counter
  integer   :: nen1                         ! energy counter
  integer   :: nex                          ! discrete level
  integer   :: nex0                         ! base number for discrete level
  integer   :: nin                          ! counter for incident energy
  integer   :: Nlin                         ! number of linearized points
  integer   :: type                         ! particle type
  real(sgl) :: Eadd(10*numenin)             ! energy (to add)
  real(sgl) :: ee(0:numenin)                ! energy
  real(sgl) :: Ehigh                        ! help variable
  real(sgl) :: Qval                         ! Q-value
  real(sgl) :: x(numenin)                   ! help variable
  real(sgl) :: xsadd(10*numenin)            ! cross section (to add)
!
! **************************** Make MF3 ********************************
!
! Loop over ejectiles
!
  MF = 3
  do type = 1, 6
    if (type == 1) then
      if (k0 == 1) then
        MT1 = 51
      else
        MT1 = 50
      endif
      MTc = 91
    else
      MT1 = 500 + 50 * type
      MTc = MT1 + 49
    endif
!
! 1. Discrete cross sections
!
    MT = MT1 - 1
    MTlast = MT
    if ((MT >= 49 .and. MT <= 90) .and. .not. flagparn) cycle
    if ((MT >= 599 .and. MT <= 648) .and. .not. flagparp) cycle
    if ((MT >= 649 .and. MT <= 698) .and. .not. flagpard) cycle
    if ((MT >= 699 .and. MT <= 748) .and. .not. flagpart) cycle
    if ((MT >= 749 .and. MT <= 798) .and. .not. flagparh) cycle
    if ((MT >= 799 .and. MT <= 848) .and. .not. flagpara) cycle
    do nex = 0, nlevmax
      if (k0 /= 0 .and. type == 1 .and. nex == Ltarget) cycle
      MT = MT + 1
      Qval = Qdisc(type, 0)
      QM(MT) = Qval * 1.e6
      if (mod(MT, 50) == 0) then
        nex0 = 0
      else
        nex0 = nex
      endif
      QI(MT) = Qdisc(type, nex0) * 1.e6
      if (k0 == 0 .and. Qdisc(type, nex0) >= 0.) then
       iE = 0
      else
        EthMT(MT) = max(Ethdisc(type, nex0), EminMeV)
        if (nlevmax == 0 .and. k0 == type) QI(MT) = EthMT(MT) * 1.e6 * AWR / (AWR + 1.)
        E3(MT, 1) = EthMT(MT) * 1.e6
        if (eninc(1) == EminMeV) then
          xs(MT, 1) = xsdisc(type, nex, 1) * 1.e-3
        else
          xs(MT, 1) = 0.
        endif
        iE = 1
      endif
!
! Linearize exothermic partial cross sections
!
      Ehigh = 0.
      if (Qdisc(type, nex) > 0.) then
        if (EHres > 0.) then
          Ehigh = EHres
        else
          Ehigh = 1000.
        endif
        do nin = 1, numendisc(type)
          ee(nin) = edisc(type, nin)
        enddo
        call locate(ee, 1, numendisc(type), Ehigh * 1.e-6, nen)
        if (nen >= 1) then
          do i = 1, nen
            x(i) = xsintdisc(type, nex, i)
          enddo
          Nlin = 100
          nen1 = min(nen, numinc)
          call linear(ee, x, nen1, Eadd, xsadd, Nlin)
          if (nen == 1) Nlin = 1
          iE = 0
          do i = 1, Nlin
            iE = iE + 1
            E3(MT, iE) = Eadd(i) * 1.e6
            xs(MT, iE) = xsadd(i) * 1.e-3
          enddo
        endif
      endif
      do nin = 1, numendisc(type)
        if (edisc(type, nin) <= EthMT(MT)) cycle
        if (edisc(type, nin) <= EHres * 1.e-6) cycle
        if (edisc(type, nin) <= Ehigh * 1.e-6) cycle
        if (edisc(type, nin) > Eswitch) cycle
        if (edisc(type, nin) > EMAX * 1.e-6) cycle
        if (xsintdisc(type, nex, nin) == 0..and. xsintdisc(type, nex, nin + 1) == 0.) cycle
        iE = iE + 1
        E3(MT, iE) = edisc(type, nin) * 1.e6
        xs(MT, iE) = xsintdisc(type, nex, nin) * 1.e-3
      enddo
      if (iE > 1) then
        MTlast = MT
      else
        MT = MT - 1
        cycle
      endif
      if (xsdisc(type, nex, numcut) == 0.) then
        iE = iE + 1
        E3(MT, iE) = eninccut * 1.e6
        xs(MT, iE) = 0.
      endif
!
! High energies
!
      if (flaghigh) then
        iE = iE + 1
        E3(MT, iE) = eninccut * 1.e6 + cuteps
        xs(MT, iE) = 0.
        iE = iE + 1
        E3(MT, iE) = EMAX
        xs(MT, iE) = 0.
      endif
!
! Prevent overlap of any competing cross section with RR
!
      flagU = .true.
      if (urrcomp == 2) flagU = .false.
      if (urrcomp == 1 .and. (MT < 51 .or. MT > 91)) flagU = .false.
      if (E3(MT, 1) > EHres) flagU = .false.
      if (flagU) then
        nin = iE
        do i = 1, nin
          if (E3(MT, i) <= EHres) then
            xs(MT, i) = 0.
          else
            do j = nin + 1, i + 1, -1
              E3(MT, j) = E3(MT, j - 1)
              xs(MT, j) = xs(MT, j - 1)
            enddo
            E3(MT, i) = EHres
            xs(MT, i) = 0.
            iE = iE + 1
            exit
          endif
        enddo
      endif
!
! Remove MT number if all cross sections are zero
!
      do i = 1, iE
        if (xs(MT, i) > 0.) goto 100
      enddo
      mtexist(MF, MT) = .false.
      cycle
!
! ENDF-6 parameters
!
  100     NP(MF, MT) = iE
      NR(MF, MT) = 1
      NBT(MF, MT, 1) = iE
      INTER(MF, MT, 1) = 2
      if (xs(MT, iE) == 0.) then
        INTER(MF, MT, NR(MF, MT)) = 2
        NBT(MF, MT, NR(MF, MT)) = NP(MF, MT)
      endif
      LR3(MT) = 0
      NMT = NMT + 1
      if (NMT > NMTmax) cycle
      mtexist(MF, MT) = .true.
      mfexist(MF) = .true.
    enddo
!
! 2. Continuum cross sections
!
    MT = MTc
    QM(MT) = Qdisc(type, 0) * 1.e6
    if (MTlast == MT1 - 1) then
      QI(MT) = QM(MT)
      EthMT(MT) = max(Ethdisc(type, 0), EminMeV)
    else
      QI(MT) = QI(MTlast)
      EthMT(MT) = EthMT(MTlast)
    endif
    if (k0 == 1 .and. QI(MT) > 0.) then
      if (EHres > 0.) then
        Ehigh = EHres
      else
        Ehigh = 1000.
      endif
      call locate(eninc, 1, numinc, Ehigh * 1.e-6, nen)
      if (nen >= 1) then
        do i = 1, nen + 1
          x(i) = xscont(type, i)
        enddo
        Nlin = 100
        nen1 = min(nen + 1, numinc)
        call linear(eninc, x, nen1, Eadd, xsadd, Nlin)
        iE = 0
        do i = 1, Nlin
          if (Eadd(i) <= Ehigh * 1.e-6) then
            iE = iE + 1
            E3(MT, iE) = Eadd(i) * 1.e6
            xs(MT, iE) = xsadd(i) * 1.e-3
          endif
        enddo
      endif
    else
      E3(MT, 1) = EthMT(MT) * 1.e6
      if (EthMT(MT) == EminMeV .and. eninc(1) == EminMeV) then
        xs(MT, 1) = xscont(type, 1) * 1.e-3
      else
        xs(MT, 1) = 0.
      endif
      iE = 1
    endif
    if (k0 == 0 .and. Qdisc(type, 0) >= 0.) iE = 0
    do nin = 1, numcut
      if (eninc(nin) < EthMT(MT)) cycle
      if (k0 == 1 .and. QI(MT) > 0..and.eninc(nin) <= Ehigh * 1.e-7) cycle
      if (eninc(nin) == EminMeV) cycle
      if (xscont(type, nin) == 0..and.xscont(type, min(nin + 1, numenin)) == 0.) cycle
      iE = iE + 1
      E3(MT, iE) = eninc(nin) * 1.e6
      xs(MT, iE) = xscont(type, nin) * 1.e-3
    enddo
    if (iE == 1) cycle
    if (xscont(type, numcut) == 0.) then
      iE = iE + 1
      E3(MT, iE) = eninccut * 1.e6
      xs(MT, iE) = 0.
    endif
!
! High energies
!
    if (flaghigh) then
      iE = iE + 1
      E3(MT, iE) = eninccut * 1.e6 + cuteps
      xs(MT, iE) = 0.
      iE = iE + 1
      E3(MT, iE) = EMAX
      xs(MT, iE) = 0.
    endif
!
! Remove MT number if all cross sections are zero
!
    do i = 1, iE
      if (xs(MT, i) > 0.) goto 300
    enddo
    mtexist(MF, MT) = .false.
    cycle
!
! ENDF-6 parameters
!
  300   NP(MF, MT) = iE
    NR(MF, MT) = 1
    NBT(MF, MT, 1) = iE
    INTER(MF, MT, 1) = 2
    LR3(MT) = 0
    NMT = NMT + 1
    if (NMT > NMTmax) cycle
    mtexist(MF, MT) = .true.
    mfexist(MF) = .true.
  enddo
  return
end subroutine make3discrete
! Copyright A.J. Koning 2021
