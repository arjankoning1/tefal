subroutine make3mt5
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3 for MT5
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
! All global variables
!   numenin        ! number of incident energies
! Variables for input of ENDF library type
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps         ! energy shift at MT5 cutoff energy (in eV)
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
!   NMTmax         ! maximum number of MT numbers
! Variables for reaction initialization
!   eninccut       ! last incident energy before high - energy format
!   numcut         ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   numinc         ! number of incident energies
! Variables for initialization of ENDF format
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
! Variables for total cross sections in ENDF format
!   xsany          ! (x, anything) cross section (MT5)
!   xsnonel        ! nonelastic cross section
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   xsexcl         ! exclusive cross section
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF3
!   E3             ! incident energy for MF3 (in ENDF - 6 format)
!   E3res          ! energy in resonance range
!   EthMT          ! threshold energy
!   LR3            ! breakup flag
!   NE3res         ! number of energies in resonance range
!   NMT            ! total number of MT sections
!   QI             ! Q - value (in ENDF - 6 format)
!   QM             ! Q - value (in ENDF - 6 format)
!   xs             ! cross section
!
! *** Declaration of local data
!
  implicit none
  integer :: iE               ! energy counter
  integer :: MF               ! MF-number
  integer :: MT               ! MT-number
  integer :: nin              ! counter for incident energy
!
! ************************* Make MF3 for MT5 ***************************
!
  MF = 3
  MT = 5
  QM(MT) = 0.
  QI(MT) = 0.
  iE = 0
!
! A. Detailed ENDF-6 representation
!
  if (flagendfdet) then
    if (k0 == 0) then
      iE = 0
      EthMT(MT) = EminMeV
    endif
    if (k0 == 1) then
      iE = 1
      EthMT(MT) = EminMeV
      E3(MT, iE) = EthMT(MT) * 1.e6
      xs(MT, iE) = 0.
    endif
    if (k0 > 1 .and. Eswitch > 0.) then
      iE = 1
      EthMT(MT) = eninc(1)
      E3(MT, iE) = EthMT(MT) * 1.e6
      xs(MT, iE) = 0.
    endif
!
! 1. In the absence of any alternative within the ENDF-6 format, we store the (n,gn), (n,gp)...,(n,galpha) cross sections in MT5.
!    These generally have rather small values.
!
    if (flaggpf) then
      if (NE3res(1) /= 0) then
        iE = iE + 1
        E3(MT, iE) = E3res(1, NE3res(1))
        xs(MT, iE) = 0.
      endif
      do nin = 1, numcut
        if (nin /= numcut .and. xsany(nin) == 0..and. xsany(min(nin + 1, numenin)) == 0.) cycle
        if (eninc(nin) * 1.e6 < E3res(1, NE3res(1))) cycle
        if (eninc(nin) == EminMeV .and. iE == 1) cycle
        iE = iE + 1
        E3(MT, iE) = eninc(nin) * 1.e6
        xs(MT, iE) = real(xsany(nin)) * 1.e-3
      enddo
    else
      if (flaghigh) then
        iE = iE + 1
        E3(MT, iE) = eninccut * 1.e6
        xs(MT, iE) = 0.
      endif
    endif
!
! 2. Reaction cross sections for high energies
!
    if (flaghigh .and. numcut > 1) then
      iE = iE + 1
      E3(MT, iE) = eninccut * 1.e6 + cuteps
      xs(MT, iE) = real(xsnonel(numcut)) * 1.e-3
      if (flagfission .and. .not. flagfis10) xs(MT, iE) = xs(MT, iE) - real(dble(xsexcl(-1, nin)) * 1.e-3)
    endif
    do nin = numcut + 1, numinc
      if (nin < numinc) then
        if (xsany(nin) == 0..and.xsany(nin + 1) == 0.) cycle
      endif
      iE = iE + 1
      E3(MT, iE) = eninc(nin) * 1.e6
      xs(MT, iE) = real(xsany(nin)) * 1.e-3
    enddo
  else
    do nin = 1, numinc
      if (nin < numinc) then
        if (xsany(nin) == 0..and.xsany(nin + 1) == 0.) cycle
      endif
      iE = iE + 1
      E3(MT, iE) = eninc(nin) * 1.e6
      xs(MT, iE) = real(xsany(nin)) * 1.e-3
    enddo
  endif
!
! Final check for non-zero values
!
  NMT = NMT + 1
  if (NMT > NMTmax) return
  do nin = 1, numinc
    if (xs(MT, nin) > 0.) then
      mtexist(MF, MT) = .true.
      mfexist(MF) = .true.
      exit
    endif
  enddo
!
! 3. ENDF-6 parameters
!
  NP(MF, MT) = iE
  NR(MF, MT) = 1
  NBT(MF, MT, 1) = iE
  INTER(MF, MT, 1) = 2
  LR3(MT) = 0
  return
end subroutine make3mt5
! Copyright A.J. Koning 2021
