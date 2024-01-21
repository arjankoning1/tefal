subroutine processphoton
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process gamma data from TALYS
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
!   numchan         ! maximum number of exclusive channels
!   numen2          ! number of emission energies
!   numengam        ! number of incident energies for gamma cross sections
!   numenspec       ! number of incident energies for spectra
!   numgam          ! number of gamma lines
!   numlevels       ! maximum number of discrete levels
! Variables for input of ENDF structure
!   flagcapt6       ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flaggamspec     ! flag to store gamma prod. only as a spectrum and not per level
!   flagpara        ! flag to include partial cross sections for alphs
!   flagpard        ! flag to include partial cross sections for deuterons
!   flagparh        ! flag to include partial cross sections for helions
!   flagparp        ! flag to include partial cross sections for protons
!   flagpart        ! flag to include partial cross sections for tritons
!   flagpart6       ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
!   flagrenorm      ! flag for renormalization of spectra
! Constants
!   parA            ! mass number of particle
!   xsepshigh       ! upper limit for cross sections in millibarns
!   xsepslow        ! lower limit for cross sections in millibarns
! Variables for reaction initialization
!   Egamindex       ! enegy index for gamma cross sections
!   massN           ! mass of nucleus in neutron units
!   Nengam          ! number of incident energies for gamna c.s.
!   nlevmax         ! number of included discrete levels
! Variables for info from TALYS
!   eninc           ! incident energy
!   k0              ! index of incident particle
! Variables for initialization of ENDF format
!   idnum           ! number of different exclusive cross sections
! Variables for partial cross sections in ENDF format
!   Eparticles      ! total energy carried away by particles
!   Erecav          ! average recoil energy
!   Estartdis       ! starting level
!   f0ex            ! energy distribution for exclusive channel
!   Egammadis       ! gamma energy
!   Ehist           ! histogram emission energy
!   idchannel       ! identifier for channel
!   nbeg            ! first outgoing energy
!   nend            ! last outgoing energy
!   nout            ! number of emission energies
!   Qexcl           ! Q - value
!   specexcl        ! exclusive spectra
!   xsexcl          ! exclusive cross section
!   xsgamdis        ! exclusive discrete gamma - ray cross section
!   xsgamexcl       ! exclusive gamma cross section
! Variables for discrete state cross sections in ENDF format
!   xscont          ! continuum cross section
!   xsdisc          ! discrete state cross section
! Variables for photon production in ENDF format
!   Egamma          ! gamma energy
!   Estart          ! starting level
!   Ngam            ! number of gamma transitions for this nucleus
!   xsgam           ! exclusive discrete gamma - ray cross section
!   xsgamcont       ! gamma production cross section per residual nucleus
!   xsgamdisctot    ! total discrete gamma - ray cross section
!   xsgamtot        ! total (discrete + continuum) gamma - ray cross section
!   yieldcont       ! continuum gamma - ray yield
!   yielddisc       ! discrete gamma - ray yield
!   yieldg          ! total discrete gamma yield per level
!   yieldgam        ! gamma yield
!   yieldtot        ! total gamma - ray yield
!
! *** Declaration of local data
!
  implicit none
  logical   :: gamexist(0:numlevels, 0:numlevels)            ! logical for existence of gamma ray
  logical   :: limit                                         ! help variable
  integer   :: i1                                            ! value
  integer   :: i2                                            ! value
  integer   :: id                                            ! counter for deuterons
  integer   :: idc                                           ! help variable
  integer   :: igam                                          ! counter for gammas
  integer   :: igam2                                         ! counter
  integer   :: indexgam(0:numlevels, 0:numlevels)            ! new, one-dimensional index for gamma ray
  integer   :: nbeg0                                         ! first outgoing energy
  integer   :: nen                                           ! energy counter
  integer   :: nend0                                         ! last outgoing energy
  integer   :: nenout                                        ! counter for outgoing energy
  integer   :: nex                                           ! discrete level
  integer   :: Ng                                            ! number of discrete gamma-rays
  integer   :: nin                                           ! counter for incident energy
  integer   :: type                                          ! particle type
  real(sgl) :: Eavailable                                    ! available emission energy
  real(sgl) :: Eaverage                                      ! average emission energy
  real(sgl) :: Egammatmp                                     ! gamma energy
  real(sgl) :: eninccm                                       ! incident energy in C.M. frame
  real(sgl) :: Estarttmp                                     ! starting level
  real(sgl) :: Esum                                          ! total photon emission energy
  real(sgl) :: factor                                        ! multiplication factor
  real(sgl) :: gspec                                         ! gamma spectrum
  real(sgl) :: specmass                                      ! specific mass
  real(sgl) :: tot                                           ! total
  real(sgl) :: totC                                          ! total
  real(sgl) :: totE                                          ! total
  real(sgl) :: x                                             ! help variable
  real(sgl) :: xsex                                          ! help variable
  real(sgl) :: xsgamtmp(numengam)                            ! gamma cross section
  real(sgl) :: xsy                                           ! help variable
!
! *********************** Re-indexing of gamma rays ********************
!
  do idc = 0, numchan
    do i1 = 1, numlevels
      do i2 = 0, i1
        gamexist(i1, i2) = .false.
        indexgam(i1, i2) = 0
      enddo
    enddo
    Ng = 0
    do i1 = 1, numlevels
      do i2 = 0, i1
        do nen = 1, Nengam
          if (xsgamdis(idc, nen, i1, i2) > 0..and. Egammadis(idc, i1, i2) > 0.001 .and. Ng < numgam) then
            if ( .not. gamexist(i1, i2)) then
              gamexist(i1, i2) = .true.
              Ng = Ng + 1
              indexgam(i1, i2) = Ng
            endif
            igam = indexgam(i1, i2)
            xsgam(idc, igam, nen) = xsgamdis(idc, nen, i1, i2)
            Egamma(idc, igam) = Egammadis(idc, i1, i2)
            Estart(idc, igam) = Estartdis(idc, i1, i2)
          endif
        enddo
      enddo
    enddo
    Ngam(idc) = Ng
  enddo
!
! ********** Sort gamma ray energies into descending order *************
!
  do idc = 0, numchan
    do igam = 1, Ngam(idc)
      do igam2 = 1, igam
        if (Egamma(idc, igam) > Egamma(idc, igam2)) then
          Egammatmp = Egamma(idc, igam2)
          Estarttmp = Estart(idc, igam2)
          do nen = 1, Nengam
            xsgamtmp(nen) = xsgam(idc, igam2, nen)
          enddo
          Egamma(idc, igam2) = Egamma(idc, igam)
          Estart(idc, igam2) = Estart(idc, igam)
          do nen = 1, Nengam
            xsgam(idc, igam2, nen) = xsgam(idc, igam, nen)
          enddo
          Egamma(idc, igam) = Egammatmp
          Estart(idc, igam) = Estarttmp
          do nen = 1, Nengam
            xsgam(idc, igam, nen) = xsgamtmp(nen)
          enddo
        endif
      enddo
    enddo
!
! ************************** Determine limits **************************
!
    do igam = 1, Ngam(idc)
      limit = .false.
      do nen = 1, Nengam
        if (xsgam(idc, igam, nen) < xsepslow) xsgam(idc, igam, nen) = 0.
        if (xsgam(idc, igam, nen) >= xsepshigh) limit = .true.
      enddo
      if ( .not. limit) then
        do nen = 1, Nengam
          xsgam(idc, igam, nen) = 0.
        enddo
      endif
    enddo
!
! **** Construct total discrete gamma-ray cross section per nucleus ****
!
    do nen = 1, Nengam
      nin = Egamindex(nen)
      xsgamdisctot(idc, nen) = 0.
      do igam = 1, Ngam(idc)
        if (Egamma(idc, igam) >= 20.) cycle
        xsgamdisctot(idc, nen) = xsgamdisctot(idc, nen) + xsgam(idc, igam, nen)
      enddo
!
! Construct total (discrete+continuum) gamma-ray cross section per nucleus
!
      xsgamtot(idc, nen) = xsgamexcl(idc, nin)
      xsgamcont(idc, nen) = xsgamtot(idc, nen) - xsgamdisctot(idc, nen)
    enddo
!
! ************************** Determine limits **************************
!
    limit = .false.
    do nen = 1, Nengam
      if (xsgamdisctot(idc, nen) < xsepslow) xsgamdisctot(idc, nen) = 0.
      if (xsgamdisctot(idc, nen) >= xsepshigh) limit = .true.
    enddo
    if ( .not. limit) then
      do nen = 1, Nengam
        xsgamdisctot(idc, nen) = 0.
      enddo
    endif
    limit = .false.
    do nen = 1, Nengam
      if (xsgamcont(idc, nen) < xsepslow) xsgamcont(idc, nen) = 0.
      if (xsgamcont(idc, nen) >= xsepshigh) limit = .true.
    enddo
    if ( .not. limit) then
      do nen = 1, Nengam
        xsgamcont(idc, nen) = 0.
      enddo
    endif
    limit = .false.
    do nen = Nengam, 1, -1
      if (xsgamtot(idc, nen) < xsepslow) xsgamtot(idc, nen) = 0.
      if (xsgamtot(idc, nen) >= xsepshigh) limit = .true.
      if (idchannel(idc) == 0 .and. xsgamtot(idc, nen) == 0..and. &
        nen <= Nengam - 1) xsgamtot(idc, nen) = xsgamtot(idc, nen + 1)
    enddo
    if ( .not. limit) then
      do nen = 1, Nengam
        xsgamtot(idc, nen) = 0.
      enddo
    endif
  enddo
!
! **************** Construct exclusive gamma-ray yields ****************
!
  do idc = 0, idnum
    id = idchannel(idc)
    do nen = 1, Nengam
      nin = Egamindex(nen)
      xsex = xsexcl(idc, nin)
      if (flagparn .and. id == 100000) xsex = xscont(1, nin)
      if (flagparp .and. id == 10000) xsex = xscont(2, nin)
      if (flagpard .and. id == 1000) xsex = xscont(3, nin)
      if (flagpart .and. id == 100) xsex = xscont(4, nin)
      if (flagparh .and. id == 10) xsex = xscont(5, nin)
      if (flagpara .and. id == 1) xsex = xscont(6, nin)
      yieldgam(idc, nen) = 0.
      if (xsex == 0.) cycle
      if ((flagpart6 .and. id /= 0) .or. (flagcapt6 .and. id == 0)) then
        x = xsgamtot(idc, nen)
        if (id == 1 .or. id == 10 .or. id == 100 .or. id == 1000 .or. id == 10000 .or. id == 100000) then
          type = 6 - int(log10(id + 1.))
          do nex = 1, nlevmax
            x = x - yieldg(type, nex) * xsdisc(type, nex, nin)
          enddo
          x = max(x, 0.)
        endif
      else
        x = xsgamcont(idc, nen)
      endif
      if (xsex <= 1.e-9) then
        yieldgam(idc, nen) = 1.
      else
        yieldgam(idc, nen) = x / xsex
      endif
    enddo
  enddo
!
! ************************* Energy spectra *****************************
!
  do idc = 0, idnum
    do nen = Nengam, 1, -1
      nin = Egamindex(nen)
      xsex = xsexcl(idc, nin)
      id = idchannel(idc)
      if (id == 100000) xsex = xscont(1, nin)
      if (id == 10000) xsex = xscont(2, nin)
      if (id == 1000) xsex = xscont(3, nin)
      if (id == 100) xsex = xscont(4, nin)
      if (id == 10) xsex = xscont(5, nin)
      if (id == 1) xsex = xscont(6, nin)
      if (xsex == 0.) cycle
!
! Exclusive photon spectra
!
      nbeg0 = 0
      nend0 = 0
      if (xsgamcont(idc, nen) > 0.) then
!
! Determine first and last outgoing energy
!
        nbeg0 = 1
        nend0 = nout(idc, nen)
        do nenout = 2, nout(idc, nen)
          if (specexcl(idc, nen, 0, nenout) > 0.) then
            nbeg0 = nenout - 1
            exit
          endif
        enddo
        do nenout = nout(idc, nen), 1, -1
          if (specexcl(idc, nen, 0, nenout) > 0.) then
            nend0 = nenout + 1
            exit
          endif
        enddo
      endif
      nend0 = min(nend0, nout(idc, nen))
      nbeg(idc, 0, nen) = nbeg0
      nend(idc, 0, nen) = nend0
!
! Normalization of spectra
!
      if ((flagpart6 .and. id /= 0) .or. (flagcapt6 .and. id == 0)) then
        xsy = xsgamdisctot(idc, nen)
      else
        xsy = 0.
      endif
      nend0 = nend0 + 1
      nend(idc, 0, nen) = nend0
      f0ex(idc, nen, 0, nend0) = 0.
      do nenout = max(nbeg0 - 1, 0), nend0 - 1
        xsy = xsy + (Ehist(idc, nen, 0, nenout + 1) - Ehist(idc, nen, 0, nenout)) * specexcl(idc, nen, 0, nenout)
      enddo
      if (xsy == 0.) then
        xsgamdisctot(idc, nen) = 0.
        xsgamcont(idc, nen) = 0.
        nbeg(idc, 0, nen) = 0
        nend(idc, 0, nen) = 0
        cycle
      endif
      if ( .not. flagrenorm) xsy = xsgamdisctot(idc, nen) + xsgamcont(idc, nen)
      Ng = Ngam(idc)
      do nenout = nbeg0, nend0
        f0ex(idc, nen, 0, nenout) = max(specexcl(idc, nen, 0, nenout) / xsy, 1.e-13)
        if (flaggamspec) then
          do igam = 1, Ng
            if (Egamma(idc, igam) >= Ehist(idc, nen, 0, nenout) .and. &
              Egamma(idc, igam) < Ehist(idc, nen, 0, nenout + 1)) then
              gspec = xsgam(idc, igam, nen) / (Ehist(idc, nen, 0, nenout + 1) - Ehist(idc, nen, 0, nenout))
              f0ex(idc, nen, 0, nenout) = f0ex(idc, nen, 0, nenout) + gspec / xsy
            endif
          enddo
        endif
      enddo
      if (flaggamspec) Ng = 0
      if ((flagpart6 .and. id /= 0) .or. (flagcapt6 .and. id == 0)) then
        yieldtot(idc, nen) = 0.
        yieldcont(idc, nen) = 0.
        do igam = 1, Ng
          yielddisc(idc, igam, nen) = 0.
          if (xsy > 0.) yielddisc(idc, igam, nen) = xsgam(idc, igam, nen) / xsy
        enddo
      endif
      if ( .not. flagcapt6 .and. id == 0) then
        if (xsex /= 0.) then
          yieldtot(idc, nen) = xsgamtot(idc, nen) / xsex
          yieldcont(idc, nen) = xsgamcont(idc, nen) / xsex
          do igam = 1, Ng
            yielddisc(idc, igam, nen) = xsgam(idc, igam, nen) / xsex
          enddo
          if (xsgamcont(idc, nen) /= 0.) then
            do nenout = nbeg0, nend0
              f0ex(idc, nen, 0, nenout) = max(specexcl(idc, nen, 0, nenout) / xsgamcont(idc, nen), 1.e-13)
            enddo
          endif
        endif
      endif
!
! Extrapolate gamma production for E < 1 keV
!
      if (eninc(nin) < 0.001) then
        xsgamdisctot(idc, nen) = xsgamdisctot(idc, nen + 1)
        xsgamcont(idc, nen) = xsgamcont(idc, nen + 1)
        yieldtot(idc, nen) = yieldtot(idc, nen + 1)
        yieldcont(idc, nen) = yieldcont(idc, nen + 1)
        do igam = 1, Ng
          yielddisc(idc, igam, nen) = yielddisc(idc, igam, nen + 1)
        enddo
        nbeg(idc, 0, nen) = nbeg(idc, 0, nen + 1)
        nend(idc, 0, nen) = nend(idc, 0, nen + 1)
        do nenout = nbeg(idc, 0, nen), nend(idc, 0, nen)
          Ehist(idc, nen, 0, nenout) = Ehist(idc, nen + 1, 0, nenout)
          f0ex(idc, nen, 0, nenout) = f0ex(idc, nen + 1, 0, nenout)
        enddo
      endif
!
! Renormalize energy balance by adjusting total gamma multiplicity
!
      if (flagrenorm .and. yieldgam(idc, nen) > 0.) then
        tot = 0.
        totE = 0.
        do igam = 1, Ng
          tot = tot + yielddisc(idc, igam, nen)
          totE = totE + Egamma(idc, igam) * yielddisc(idc, igam, nen)
        enddo
        totC = 0.
        do nenout = max(nbeg0 - 1, 0), nend0 - 1
          totC = totC + (Ehist(idc, nen, 0, nenout + 1) - Ehist(idc, nen, 0, nenout)) * f0ex(idc, nen, 0, nenout)
          totE = totE + 0.5 * (Ehist(idc, nen, 0, nenout) + Ehist(idc, nen, 0, nenout + 1)) * &
 &          (Ehist(idc, nen, 0, nenout + 1) - Ehist(idc, nen, 0, nenout)) * f0ex(idc, nen, 0, nenout)
        enddo
        if ( .not. flagcapt6 .and. id == 0 .and. totC > 0.) then
          do nenout = max(nbeg0 - 1, 0), nend0 - 1
            f0ex(idc, nen, 0, nenout) = f0ex(idc, nen, 0, nenout) / totC
          enddo
        endif
        tot = tot + totC
        if (totE > 0..and.tot > 0.) then
          specmass = massN / (massN + parA(k0))
          eninccm = eninc(nin) * specmass
          Eavailable = max(eninccm + Qexcl(idc) - Eparticles(idc, nen), 0.)
          if (Erecav(idc, nen) < Eavailable) Eavailable = Eavailable - Erecav(idc, nen)
          Eaverage = totE / tot
          Esum = Eaverage * yieldgam(idc, nen)
          factor = Eavailable / Esum
          yieldgam(idc, nen) = yieldgam(idc, nen) * factor
          tot = Eaverage * yieldgam(idc, nen)
        endif
      endif
    enddo
    if (flaggamspec) Ngam(idc) = 0
  enddo
!
! Set some histogram energies
!
  do idc = 0, numchan
    do nen = numenspec - 1, 1, -1
      do type = 0, 6
        do nenout = 0, numen2
          if (Ehist(idc, nen, type, nenout) == 0.) Ehist(idc, nen, type, nenout) = Ehist(idc, nen + 1, type, nenout)
        enddo
      enddo
      nend0 = nend(idc, 0, nen)
      if (nend0 > 1 .and. Ehist(idc, nen, 0, nend0) == 0.) nend(idc, 0, nen) = nend(idc, 0, nen) - 1
    enddo
  enddo
  return
end subroutine processphoton
! Copyright A.J. Koning 2021
