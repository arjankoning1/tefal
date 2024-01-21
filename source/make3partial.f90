subroutine make3partial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3 for partial cross sections
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
!   numenin        ! number of incident energies
!   nummt          ! number of MT numbers
! Variables for input of specific ENDF data
!   urrcomp        ! mode for competition in the URR, 0: none, 1:MT4, 2:all
! Variables for input of ENDF library type
!   flageaf        ! flag for EAF - formatted activation library
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps         ! energy shift at MT5 cutoff energy (in eV)
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
!   NMTmax         ! maximum number of MT numbers
! Variables for input of ENDF structure
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   flagmulti      ! flag to include multi - chance fission
!   flagpara       ! flag to include partial cross sections for alphs
!   flagpard       ! flag to include partial cross sections for deuterons
!   flagparp       ! flag to include partial cross sections for protons
!   flagpart       ! flag to include partial cross sections for tritons
! Variables for reaction initialization
!   EHres          ! upper energy in resonance range
!   eninccut       ! last incident energy before high - energy format
!   numcut         ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   numinc         ! number of incident energies
! Variables for initialization of ENDF format
!   idnum          ! number of different exclusive cross sections
!   LIS            ! state number of target nucleus
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
!   MTid           ! channel identifier for MT - number
!   MTinel         ! MT - number for inelastic scattering
! Variables for partial cross sections in ENDF format
!   Ethexcl        ! threshold energy
!   flagfission    ! flag for fission
!   idchannel      ! identifier for channel
!   Qexcl          ! Q - value
!   reacstring     ! string with reaction information
!   xsexcl         ! exclusive cross section
! Variables for discrete state cross sections in ENDF format
!   edisc          ! incident energy for discrete level cross section
!   Ethdisc        ! threshold energy
!   numendisc      ! number of incident energies including discrete energy grid
!   Qdisc          ! Q - value
!   xsbin          ! binary cross section
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF1
!   EMAX           ! upper limit of energy range for evaluation
!   ERN            ! total energy minus energy of neutrinos
! Variables for MF3
!   E3             ! incident energy for MF3 (in ENDF - 6 format)
!   E3adopt        ! incident energy adopted from other library
!   E3res          ! energy in resonance range
!   EthMT          ! threshold energy
!   LFS3           ! isomeric state number (EAF only)
!   LR3            ! breakup flag
!   mtstring       ! string with reaction information
!   NE3adopt       ! number of adopted incident energies
!   NE3res         ! number of energies in resonance range
!   NMT            ! total number of MT sections
!   QI             ! Q - value (in ENDF - 6 format)
!   QM             ! Q - value (in ENDF - 6 format)
!   xs             ! cross section
!   xs3adopt       ! cross section adopted from other library
!   xs3res         ! cross section in resonance range
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagfismt                    ! flag for fission
  logical   :: flagU                        ! flag for URR
  integer   :: i                            ! counter
  integer   :: idc                          ! help variable
  integer   :: iE                           ! energy counter
  integer   :: j                            ! counter
  integer   :: MF                           ! MF-number
  integer   :: MT                           ! MT-number
  integer   :: nen                          ! energy counter
  integer   :: nen1                         ! energy counter
  integer   :: nin                          ! counter for incident energy
  integer   :: Nlin                         ! number of linearized points
  integer   :: NN                           ! neutron number of residual nucleus
  integer   :: type                         ! particle type
  real(sgl) :: Eadd(10*numenin)             ! energy (to add)
  real(sgl) :: Eev                          ! energy in eV
  real(sgl) :: Ehigh                        ! help variable
  real(sgl) :: Eprev                        ! previous energy
  real(sgl) :: Qval                         ! Q-value
  real(sgl) :: x(numenin)                   ! help variable
  real(sgl) :: xsadd(10*numenin)            ! cross section (to add)
  real(sgl) :: xsint                        ! interpolated cross section
!
! ****************** Make MF3 for partial cross sections ***************
!
!
  MF = 3
Loop1:  do MT = 1, nummt
    if (MTid(MT) ==  -1) cycle
    if ( .not. (flagendfdet .or. flageaf) .and. MT /= 18) cycle
    if (Eswitch == 0..and.MT /= 18) cycle
    if ( .not. flaggpf .and. MT == MTinel) cycle
    if (MT >= 600) cycle
    flagfismt = (MT == 18 .or. MT == 19 .or. MT == 20 .or. MT == 21 .or. MT == 38)
    if ( .not. flagfission .and. flagfismt) cycle
    if (flagfis10 .and. flagfismt) cycle
    if (flagfismt .and. k0 > 1 .and. MT /= 18) cycle
    if (MT == 91) cycle
    do idc = -5, idnum
      if (idc <=  -1 .and. ( .not. flagfission)) cycle
      if (idc <=  -2 .and. ( .not. flagmulti .or. flageaf)) cycle
      if (idchannel(idc) == MTid(MT)) then
        iE = 0
!
! Q-values for fission
!
        if (idc <=  -1) then
          EthMT(MT) = EminMeV
          QM(MT) = ERN
          QI(MT) = QM(MT)
        else
!
! Q-values for other channels
!
          EthMT(MT) = max(Ethexcl(idc), EminMeV)
          if (MT == MTinel) EthMT(MT) = max(Ethdisc(k0, 1), EminMeV)
          if (EthMT(MT) >= eninccut) cycle Loop1
          Qval = Qexcl(idc)
          QM(MT) = Qval * 1.e6
          if (MT == MTinel .and. LIS == 0) then
            QI(MT) = Qdisc(k0, 1) * 1.e6
          else
            QI(MT) = QM(MT)
          endif
        endif
!
! Adopt background cross sections for resonance range or set them to 0
!
        if (k0 == 1 .and. .not. flageaf .and. (flagfismt .or. MT == 102)) then
          if (NE3res(MT) /= 0) then
            do i = 1, NE3res(MT)
              E3(MT, i) = E3res(MT, i)
              xs(MT, i) = xs3res(MT, i)
            enddo
            iE = NE3res(MT)
          else
            if (EHres > 0.) then
              E3(MT, 1) = EthMT(MT) * 1.e6
              xs(MT, 1) = 0.
              E3(MT, 2) = EHres
              xs(MT, 2) = 0.
              iE = 2
            else
              E3(MT, 1) = EthMT(MT) * 1.e6
              xs(MT, 1) = xsexcl(idc, 1) * 1.e-3
            endif
          endif
        else
!
! Set first energy for other channels
!
          E3(MT, 1) = EthMT(MT) * 1.e6
          if (eninc(1) == EminMeV) then
            xs(MT, 1) = xsexcl(idc, 1) * 1.e-3
          else
            xs(MT, 1) = 0.
          endif
          iE = 1
        endif
!
! Linearize exothermic partial cross sections
!
        Ehigh = 0.
        if (k0 == 1 .and. idc >= 0) then
          if (Qexcl(idc) > 0..and..not.MT == 102) then
            if (EHres > 0.) then
              Ehigh = EHres
            else
              Ehigh = 1000.
            endif
            call locate(eninc, 1, numinc, Ehigh * 1.e-6, nen)
            if (nen >= 1) then
              do  i = 1, nen + 1
                x(i) = xsexcl(idc, i)
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
          endif
        endif
!       if (k0 == 0) iE = 0
!
! MT4 and MT103-107
!
        if ( .not. flageaf .and. k0 /= 0 .and.  &
 &        ((flagparn .and. MT == 4) .or. (flagparp .and. MT == 103) .or. &
 &        (flagpard .and. MT == 104) .or. (flagpart .and. MT == 105) .or. &
 &        (flagpart .and. MT == 106) .or. (flagpara .and. MT == 107))) then
          if (MT == 4) then
            type = 1
          else
            type = MT - 101
          endif
          if (eninc(1) == EminMeV) xs(MT, 1) = xsbin(type, 1) * 1.e-3
          do nin = 1, numendisc(type)
            if (edisc(type, nin) > Eswitch) cycle
            if (edisc(type, nin) > EMAX * 1.e-6) cycle
            Eev = edisc(type, nin) * 1.e6
            if (Eev < EHres) cycle
            if (Eev < Ehigh) cycle
            if (edisc(type, nin) <= EthMT(MT)) cycle
            if (xsbin(type, nin) == 0..and.xsbin(type, nin + 1) == 0.) cycle
            iE = iE + 1
            E3(MT, iE) = Eev
            xs(MT, iE) = xsbin(type, nin) * 1.e-3
          enddo
        else
!
! Other MT numbers
!
          Eprev = 0.
          if (flagfismt) then
            NN = numinc
          else
            NN = numcut
          endif
          do nin = 1, NN
            Eev = eninc(nin) * 1.e6
            if (nin > 1 .and. Eprev < EHres .and. Eev > EHres) then
              xsint = xsexcl(idc, nin - 1) + (EHres - Eprev) / (Eev - Eprev) * &
                (xsexcl(idc, nin) - xsexcl(idc, nin - 1))
              if (xsint > 0.) then
                iE = iE + 1
                E3(MT, iE) = EHres
                xs(MT, iE) = xsint * 1.e-3
              endif
            endif
            Eprev = Eev
            if (Eev < EHres) cycle
            if (Eev < Ehigh) cycle
            if (eninc(nin) <= EthMT(MT)) cycle
            if (nin < NN .and. xsexcl(idc, nin) == 0..and. xsexcl(idc, min(nin + 1, numenin)) == 0.) cycle
            iE = iE + 1
            E3(MT, iE) = Eev
            xs(MT, iE) = xsexcl(idc, nin) * 1.e-3
          enddo
        endif
        if (iE <= 1) cycle
        if (xsexcl(idc, numcut) == 0.) then
          iE = iE + 1
          E3(MT, iE) = eninccut * 1.e6
          xs(MT, iE) = 0.
        endif
        if (flageaf) then
          LFS3(MT) = 99
          mtstring(MT) = reacstring(idc, 0)
        endif
        if (NE3adopt(MT) > 0) then
          do i = 1, NE3adopt(MT)
            E3(MT, i) = E3adopt(MT, i)
            xs(MT, i) = xs3adopt(MT, i)
          enddo
          iE = NE3adopt(MT)
        endif
        goto 100
      endif
    enddo
    cycle
!
! High energies
!
  100   if (flaghigh .and. .not. flageaf) then
      if ( .not. flagfismt) then
        iE = iE + 1
        E3(MT, iE) = eninccut * 1.e6 + cuteps
        xs(MT, iE) = 0.
        iE = iE + 1
        E3(MT, iE) = EMAX
        xs(MT, iE) = 0.
      endif
    endif
!
! Prevent overlap of any competing cross section with RR
!
    flagU = .true.
    if (urrcomp == 2) flagU = .false.
    if (urrcomp == 1 .and. MT /= 4) flagU = .false.
    if (flagfismt .or. MT == 102) flagU = .false.
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
! Remove MT number if all or only 1 cross section is zero
!
    j = 0
    do i = 1, iE
      if (xs(MT, i) > 0.) j = j + 1
      if (j > 1) goto 200
    enddo
    mtexist(MF, MT) = .false.
    cycle
!
! ENDF-6 parameters
!
  200   NP(MF, MT) = iE
    NR(MF, MT) = 1
    NBT(MF, MT, 1) = iE
    INTER(MF, MT, 1) = 2
    LR3(MT) = 0
    NMT = NMT + 1
    if ( .not. flageaf .and. NMT > NMTmax) return
    mtexist(MF, MT) = .true.
    mfexist(MF) = .true.
  enddo Loop1
  return
end subroutine make3partial
! Copyright A.J. Koning 2021
