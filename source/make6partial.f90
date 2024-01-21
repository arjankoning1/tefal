subroutine make6partial(MT)
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
! Definition of single and double precision variables
!   sgl             ! single precision kind
! Variables for input of ENDF library type
!   flagbreakup     ! breakup flag
!   flagclean       ! flag to clean up double points
!   flaghigh        ! flag for high energies ( > 20 MeV)
! Variables for input of ENDF structure
!   flagpara        ! flag to include partial cross sections for alphs
!   flagpard        ! flag to include partial cross sections for deuterons
!   flagparh        ! flag to include partial cross sections for helions
!   flagparp        ! flag to include partial cross sections for protons
!   flagpart        ! flag to include partial cross sections for tritons
!   flagpart6       ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
!   flagrecoil      ! flag to include recoil information
! Variables for reaction initialization
!   Egamindex       ! enegy index for gamma cross sections
!   Especindex      ! enegy index for spectra
!   Nengam          ! number of incident energies for gamna c.s.
!   Nenspec         ! number of incident energies for spectra
!   numcut          ! number of energies before high - energy format
!   relmass         ! mass relative to neutron mass
! Variables for info from TALYS
!   Ainit           ! mass number of initial compound nucleus
!   eninc           ! incident energy
!   k0              ! index of incident particle
!   Zinit           ! charge number of initial compound nucleus
! Constants
!   parA            ! mass number of particle
!   parmass         ! mass of particle in a.m.u.
!   parN            ! neutron number of particle
!   parZ            ! charge number of particle
!   xsepslow        ! lower limit for cross sections in millibarns
! Variables for initialization of ENDF format
!   idnum           ! number of different exclusive cross sections
!   mfexist         ! flag for existence of MF - number
!   mtexist         ! flag for existence of MT - number
!   MTid            ! channel identifier for MT - number
! Variables for residual production cross sections in ENDF format
!   nucmass         ! mass of nucleus
! Variables for partial cross sections in ENDF format
!   Ehist           ! histogram emission energy
!   Ehistrec        ! histogram recoil energy
!   f0ex            ! energy distribution for exclusive channel
!   f0exrec         ! energy distribution for recoil
!   idchannel       ! identifier for channel
!   nbeg            ! first outgoing energy
!   nbegrec         ! first outgoing energy
!   nend            ! last outgoing energy
!   nendrec         ! last outgoing energy
!   nout            ! number of emission energies
!   noutrec         ! number of recoil energies
!   Qexcl           ! Q - value
!   xsexcl          ! exclusive cross section
! Variables for discrete state cross sections in ENDF format
!   xscont          ! continuum cross section
! Variables for photon production in ENDF format
!   Egamma          ! gamma energy
!   Ngam            ! number of gamma transitions for this nucleus
!   yielddisc       ! discrete gamma - ray yield
!   yieldgam        ! gamma yield
! Variables for spectra in ENDF format
!   buratio         ! break - up ratio
!   preeqratio      ! pre - equilibrium ratio
! Variables for ENDF format
!   NK              ! number of subsections
! Variables for MF1
!   EMAX            ! upper limit of energy range for evaluation
! Variables for MF3
!   E3              ! incident energy for MF3 (in ENDF - 6 format)
!   EthMT           ! threshold energy
! Variables for MF4
!   LCT             ! LAB / CM flag
! Variables for MF6
!   AWP             ! product mass
!   b6              ! energy - angle values
!   E6              ! incident energy (in ENDF - 6 format) for distribution
!   Ey              ! incident energy for yields (in ENDF - 6 format)
!   INTER6ea        ! interpolation scheme
!   INTER6y         ! interpolation scheme
!   LANG            ! flag for angular representation
!   LAW             ! flag for distribution function
!   LEP             ! interpolation scheme for secondary energy
!   LIP             ! product modifier flag
!   NA              ! number of angular parameters
!   NBT6ea          ! separation value for interpolation scheme
!   NBT6y           ! separation value for interpolation scheme
!   ND              ! number of discrete energies
!   NE6ea           ! number of incident energies for distribution
!   NEP             ! number of secondary energy points
!   NP6y            ! number of incident energies for yields
!   NR6ea           ! number of interpolation ranges for distribution
!   NR6y            ! number of interpolation ranges for yields
!   NW              ! number of words
!   Y               ! product yield (in ENDF - 6 format)
!   ZAP             ! product identifier
!
! *** Declaration of local data
!
  implicit none
  integer   :: A                  ! mass number of target nucleus
  logical   :: flaggam            ! logical for existence of gammas
  integer   :: i                  ! counter
  integer   :: ib                 ! counter
  integer   :: id                 ! counter for deuterons
  integer   :: idc                ! help variable
  integer   :: iE                 ! energy counter
  integer   :: iyield             ! particle yield
  integer   :: k                  ! counter
  integer   :: MT                 ! MT-number
  integer   :: ncont              ! number of continua
  integer   :: nen                ! energy counter
  integer   :: nencut             ! number of energies before high-energy format
  integer   :: nenout             ! counter for outgoing energy
  integer   :: nin                ! counter for incident energy
  integer   :: Nix                ! neutron number index for residual nucleus
  integer   :: type               ! particle type
  integer   :: Z                  ! charge number of target nucleus
  integer   :: Zix                ! charge number index for residual nucleus
  real(sgl) :: yield              ! logical to use MF9 or MF10
  real(sgl) :: xsex               ! help variable
  real(sgl) :: Ein                ! incident energy
!
! ***************** Make MF6 for partial cross sections ****************
!
  if (MTid(MT) == -1) return
  if (k0 == 1 .and. MT == 103 .and. flagparp) return
  if (k0 == 1 .and. MT == 104 .and. flagpard) return
  if (k0 == 1 .and. MT == 105 .and. flagpart) return
  if (k0 == 1 .and. MT == 106 .and. flagparh) return
  if (k0 == 1 .and. MT == 107 .and. flagpara) return
  if (MT == 18 .or. MT == 19 .or. MT == 20 .or. MT == 21 .or. MT == 38) return
  mfexist(6) = .true.
  idc = -99
  do id = 0, idnum
    if (idchannel(id) == MTid(MT)) then
      idc = id
      exit
    endif
  enddo
  if (idc == -99) return
  mtexist(6, MT) = .true.
  k = 0
!
! A. Particles
!
! 1. Product yields
!
! To avoid very big files, we do not give secondary information at each incident energy, but instead skip various energy points.
!
  Zix = 0
  Nix = 0
  do type = 1, 6
    iyield = mod(MTid(MT), 10 **(7 - type)) / (10 **(6 - type))
    if (iyield == 0) cycle
    Zix = Zix + iyield * parZ(type)
    Nix = Nix + iyield * parN(type)
    k = k + 1
    ZAP(k) = 1000. * parZ(type) + parA(type)
    AWP(k) = real(relmass(type))
    yield = real(iyield)
    Y(k, 1) = yield
    Y(k, 2) = yield
    Ey(k, 1) = E3(MT, 1)
    Ey(k, 2) = EMAX
    NP6y(k) = 2
    NR6y(k) = 1
    NBT6y(k, 1) = NP6y(k)
    INTER6y(k, 1) = 2
    LAW(k) = 1
    LIP(k) = 0
    if (flagbreakup .and. (type == 1 .or. type == 2)) then
      LANG(k) = 3
    else
      LANG(k) = 2
    endif
    LEP(k) = 1
!
! 2. Energy-angle distributions
!
! The first energy is the threshold energy, with zero energy-angle distribution
!
    E6(k, 1) = E3(MT, 1)
    ND(k, 1) = 0
    if (flagbreakup .and. (type == 1 .or. type == 2)) then
      NA(k, 1) = 2
    else
      NA(k, 1) = 1
    endif
    NEP(k, 1) = 2
    NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
    b6(k, 1, 1) = 0.
    b6(k, 1, 2) = 1.
    b6(k, 1, 3) = 0.
    b6(k, 1, 4) = 1.
    b6(k, 1, 5) = 0.
    b6(k, 1, 6) = 0.
!
! Loop over incident energies
!
    iE = 1
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (nin > numcut) cycle
      Ein = eninc(nin)
      if (Ein <= EthMT(MT)) cycle
      xsex = xsexcl(idc, nin)
      if (MT == 91) xsex = xscont(1, nin)
      if (MT == 649) xsex = xscont(2, nin)
      if (MT == 699) xsex = xscont(3, nin)
      if (MT == 749) xsex = xscont(4, nin)
      if (MT == 799) xsex = xscont(5, nin)
      if (MT == 849) xsex = xscont(6, nin)
      if (xsex <= xsepslow .and. nin /= numcut) cycle
      if (nout(idc, nen) == 0 .and. nin /= numcut) cycle
      iE = iE + 1
      E6(k, iE) = Ein * 1.e6
      ND(k, iE) = 0
      if (flagbreakup .and. (type == 1 .or. type == 2)) then
        NA(k, iE) = 2
      else
        NA(k, iE) = 1
      endif
!
! Store energy-angle distributions
!
      if (nbeg(idc, type, nen) == nend(idc, type, nen)) then
        NEP(k, iE) = 2
        b6(k, iE, 1) = 0.
        b6(k, iE, 2) = 1.
        b6(k, iE, 3) = 0.
        b6(k, iE, 4) = 1.
        b6(k, iE, 5) = 0.
        b6(k, iE, 6) = 0.
      else
        ib = 0
        do nenout = nbeg(idc, type, nen), nend(idc, type, nen)
          ib = ib + 1
          b6(k, iE, ib) = Ehist(idc, nen, type, nenout) * 1.e6
          ib = ib + 1
          b6(k, iE, ib) = f0ex(idc, nen, type, nenout) * 1.e-6
          ib = ib + 1
          if (f0ex(idc, nen, type, nenout) == 0.) then
            b6(k, iE, ib) = 0.
          else
            b6(k, iE, ib) = preeqratio(type, nen, nenout)
          endif
          if (flagbreakup .and. (type == 1 .or. type == 2)) then
            ib = ib + 1
            if (f0ex(idc, nen, type, nenout) == 0.) then
              b6(k, iE, ib) = 0.
            else
              b6(k, iE, ib) = buratio(type, nen, nenout)
            endif
          endif
        enddo
        NEP(k, iE) = nend(idc, type, nen) - nbeg(idc, type, nen) + 1
      endif
      NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
    enddo
!
! High energies
!
    if (flaghigh) then
      iE = iE + 1
      E6(k, iE) = EMAX
      NEP(k, iE) = NEP(k, iE - 1)
      ND(k, iE) = ND(k, iE - 1)
      NA(k, iE) = NA(k, iE - 1)
      NW(k, iE) = NW(k, iE - 1)
      do ib = 1, NW(k, iE)
        b6(k, iE, ib) = b6(k, iE - 1, ib)
      enddo
    endif
!
! ENDF-6 parameters
!
    NE6ea(k) = iE
    NR6ea(k) = 1
    NBT6ea(k, 1) = NE6ea(k)
    INTER6ea(k, 1) = 2
  enddo
!
! B. Recoils
!
! 1. Yields
!
  if (flagrecoil .and. MT /= 102) then
    k = k + 1
    Z = Zinit - Zix
    A = Ainit - Zix - Nix
    ZAP(k) = 1000. * Z + A
    AWP(k) = real(nucmass(Zix, Nix) / parmass(1))
    Ey(k, 1) = E3(MT, 1)
    Ey(k, 2) = EMAX
    Y(k, 1) = 1.
    Y(k, 2) = 1.
    NP6y(k) = 2
    NR6y(k) = 1
    NBT6y(k, 1) = NP6y(k)
    INTER6y(k, 1) = 2
    LAW(k) = 1
    LIP(k) = 0
    LANG(k) = 1
    LEP(k) = 1
!
! 2. Recoil energy distributions
!
    E6(k, 1) = E3(MT, 1)
    NA(k, 1) = 0
    ND(k, 1) = 0
    NEP(k, 1) = 2
    NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
    b6(k, 1, 1) = 0.
    b6(k, 1, 2) = 1.
    b6(k, 1, 3) = 1.
    b6(k, 1, 4) = 0.
    iE = 1
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (nin > numcut) cycle
      Ein = eninc(nin)
      if (Ein <= EthMT(MT)) cycle
      xsex = xsexcl(idc, nin)
      if (MT == 91) xsex = xscont(1, nin)
      if (MT == 649) xsex = xscont(2, nin)
      if (MT == 699) xsex = xscont(3, nin)
      if (MT == 749) xsex = xscont(4, nin)
      if (MT == 799) xsex = xscont(5, nin)
      if (MT == 849) xsex = xscont(6, nin)
      if (xsex <= xsepslow .and. nin /= numcut) cycle
      if (noutrec(idc, nen) == 0 .and. nin /= numcut) cycle
      iE = iE + 1
      E6(k, iE) = Ein * 1.e6
!
! Store recoil energy distributions
!
      if (nbegrec(idc, nen) == nendrec(idc, nen)) then
        NEP(k, iE) = 2
        b6(k, iE, 1) = 0.
        b6(k, iE, 2) = 1.
        b6(k, iE, 3) = 1.
        b6(k, iE, 4) = 0.
      else
        ib = 0
        do nenout = nbegrec(idc, nen), nendrec(idc, nen)
          ib = ib + 1
          b6(k, iE, ib) = Ehistrec(idc, nen, nenout) * 1.e6
          ib = ib + 1
          b6(k, iE, ib) = f0exrec(idc, nen, nenout) * 1.e-6
        enddo
        NEP(k, iE) = nendrec(idc, nen) - nbegrec(idc, nen) + 1
      endif
      NA(k, iE) = 0
      ND(k, iE) = 0
      NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
    enddo
!
! High energies
!
    if (flaghigh) then
      iE = iE + 1
      E6(k, iE) = EMAX
      NEP(k, iE) = NEP(k, iE - 1)
      ND(k, iE) = ND(k, iE - 1)
      NA(k, iE) = NA(k, iE - 1)
      NW(k, iE) = NW(k, iE - 1)
      do ib = 1, NW(k, iE)
        b6(k, iE, ib) = b6(k, iE - 1, ib)
      enddo
    endif
!
! ENDF-6 parameters
!
    NE6ea(k) = iE
    NR6ea(k) = 1
    NBT6ea(k, 1) = NE6ea(k)
    INTER6ea(k, 1) = 2
  endif
!
! 1. Photon yields
!
  if (flagpart6) then
    k = k + 1
    ZAP(k) = 0.
    AWP(k) = 0.
    Ey(k, 1) = E3(MT, 1)
    if (Qexcl(idc) <= 0..or.yieldgam(idc, 1) == 0.) then
      iE = 1
    else
      iE = 0
    endif
    yield = 0.
    nencut = Nengam
    do nen = 1, Nengam
      nin = Egamindex(nen)
      if (nin > numcut) then
        nencut = nen - 1
        exit
      endif
      Ein = eninc(nin)
      if (Ein <= EthMT(MT) .and. Qexcl(idc) <= 0.) cycle
      xsex = xsexcl(idc, nin)
      if (MT == 91) xsex = xscont(1, nin)
      if (MT == 649) xsex = xscont(2, nin)
      if (MT == 699) xsex = xscont(3, nin)
      if (MT == 749) xsex = xscont(4, nin)
      if (MT == 799) xsex = xscont(5, nin)
      if (MT == 849) xsex = xscont(6, nin)
      if (xsex <= xsepslow .and. nin /= numcut) cycle
      yield = yieldgam(idc, nen)
      if (yield == 0.) cycle
      iE = iE + 1
      if (iE == 1 .and. Qexcl(idc) > 0.) then
        Ey(k, iE) = E3(MT, 1)
      else
        Ey(k, iE) = Ein * 1.e6
      endif
      Y(k, iE) = yield
    enddo
    if (iE > 1 .and. yieldgam(idc, nencut) == 0.) then
      iE = iE + 1
      Ey(k, iE) = Ein * 1.e6
      Y(k, iE) = yield
    endif
    if (flaghigh) then
      iE = iE + 1
      Ey(k, iE) = EMAX
      Y(k, iE) = yield
    endif
    if (iE >= 2) Y(k, 1) = Y(k, 2)
    NP6y(k) = iE
    NR6y(k) = 1
    NBT6y(k, 1) = NP6y(k)
    INTER6y(k, 1) = 2
    LAW(k) = 1
    LIP(k) = 0
    LANG(k) = 1
!
! 2. Photon energy distributions
!
! For non-threshold reactions (e.g. MT102), the discrete gamma yields for the first energy are not zero.
!
    if (Qexcl(idc) <= 0..or.yieldgam(idc, 1) == 0.) then
      E6(k, 1) = E3(MT, 1)
      NA(k, 1) = 0
!
! Discrete photons
!
      ND(k, 1) = Ngam(idc)
      ib = 0
      do i = 1, Ngam(idc)
        ib = ib + 1
        b6(k, 1, ib) = Egamma(idc, i) * 1.e6
        ib = ib + 1
        b6(k, 1, ib) = 0.
      enddo
      ib = 2 * ND(k, 1) + 1
      b6(k, 1, ib) = 0.
      ib = ib + 1
      b6(k, 1, ib) = 1.
      ib = ib + 1
      b6(k, 1, ib) = 1.
      ib = ib + 1
      b6(k, 1, ib) = 0.
      NEP(k, 1) = 2 + ND(k, 1)
      NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
      iE = 1
    else
      iE = 0
    endif
    do nen = 1, Nengam
      nin = Egamindex(nen)
      if (nin > numcut) exit
      Ein = eninc(nin)
      if (Ein <= EthMT(MT) .and. Qexcl(idc) <= 0.) cycle
      xsex = xsexcl(idc, nin)
      if (MT == 91) xsex = xscont(1, nin)
      if (MT == 649) xsex = xscont(2, nin)
      if (MT == 699) xsex = xscont(3, nin)
      if (MT == 749) xsex = xscont(4, nin)
      if (MT == 799) xsex = xscont(5, nin)
      if (MT == 849) xsex = xscont(6, nin)
      if (xsex <= xsepslow .and. nin /= numcut) cycle
      if (yieldgam(idc, nen) == 0.) cycle
      iE = iE + 1
      if (iE == 1 .and. Qexcl(idc) > 0.) then
        E6(k, iE) = E3(MT, 1)
      else
        E6(k, iE) = Ein * 1.e6
      endif
      NA(k, iE) = 0
      ND(k, iE) = Ngam(idc)
      ib = 0
      flaggam = .false.
      do i = 1, Ngam(idc)
        ib = ib + 1
        b6(k, iE, ib) = Egamma(idc, i) * 1.e6
        ib = ib + 1
        b6(k, iE, ib) = yielddisc(idc, i, nen)
        if (yielddisc(idc, i, nen) > 0.) flaggam = .true.
      enddo
!
! The first outgoing continuum energy is zero.
!
      ib = 2 * ND(k, iE)
      ncont = nend(idc, 0, nen) - nbeg(idc, 0, nen)
      if (ncont <= 0) then
        if ( .not. flaggam) then
          ib = ib + 1
          b6(k, iE, ib) = 0.
          ib = ib + 1
          b6(k, iE, ib) = 1.
          ib = ib + 1
          b6(k, iE, ib) = 1.
          ib = ib + 1
          b6(k, iE, ib) = 0.
          ncont = 1
        else
          ib = ib + 1
          b6(k, iE, ib) = 0.
          ib = ib + 1
          b6(k, iE, ib) = 0.
          ib = ib + 1
          b6(k, iE, ib) = 1.
          ib = ib + 1
          b6(k, iE, ib) = 0.
          ncont = 1
        endif
      else
!
! Store energy distributions
!
        do nenout = nbeg(idc, 0, nen), nend(idc, 0, nen)
          ib = ib + 1
          b6(k, iE, ib) = Ehist(idc, nen, 0, nenout) * 1.e6
          ib = ib + 1
          b6(k, iE, ib) = f0ex(idc, nen, 0, nenout) * 1.e-6
        enddo
      endif
      NEP(k, iE) = ncont + 1 + ND(k, iE)
      NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
    enddo
    if (iE >= 1 .and. yieldgam(idc, nencut) == 0.) then
      iE = iE + 1
      E6(k, iE) = Ein * 1.e6
      NEP(k, iE) = NEP(k, iE - 1)
      ND(k, iE) = ND(k, iE - 1)
      NA(k, iE) = NA(k, iE - 1)
      NW(k, iE) = NW(k, iE - 1)
      do ib = 1, NW(k, iE)
        b6(k, iE, ib) = b6(k, iE - 1, ib)
      enddo
    endif
!
! High energies
!
    if (iE >= 1 .and. flaghigh) then
      iE = iE + 1
      E6(k, iE) = EMAX
      NEP(k, iE) = NEP(k, iE - 1)
      ND(k, iE) = ND(k, iE - 1)
      NA(k, iE) = NA(k, iE - 1)
      NW(k, iE) = NW(k, iE - 1)
      do ib = 1, NW(k, iE)
        b6(k, iE, ib) = b6(k, iE - 1, ib)
      enddo
    endif
!
! ENDF-6 parameters
!
! make6clean: subroutine to clean up MF6
! write6    : subroutine to write MF6
!
    NE6ea(k) = iE
    NR6ea(k) = 1
    NBT6ea(k, 1) = NE6ea(k)
    INTER6ea(k, 1) = 2
  endif
  NK(6, MT) = k
  if (k == 0) return
  if (iE <= 1 .or. NP6y(k) <= 1) NK(6, MT) = k - 1
  if (flagrecoil .or. MT == 102) then
    LCT = 3
  else
    LCT = 2
  endif
  LEP(k) = 1
  if (flagclean) call make6clean(MT)
  call write6(MT)
  return
end subroutine make6partial
! Copyright A.J. Koning 2021
