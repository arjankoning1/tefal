subroutine make2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF2
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
!   numjres        ! number of j - values
!   numlres        ! number of l - values
!   numnrs         ! number of resonances
!   numres         ! number of resonance sections
! Variables for input of specific ENDF data
!   adopt          ! logical for existence of MF information (per MT)
!   lssfinp        ! 0: URR cross section from MF2, 1: URR cross section
!   urrenergy      ! upper energy of the URR in MeV
!   urrmode        ! 0: no URR, 1: URR from TALYS, 2: URR from data library
! Variables for input of ENDF structure
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   flagsubfis     ! flag to include subactinide fission
! Variables for reaction initialization
!   EHres          ! upper energy in resonance range
! Variables for info from TALYS
!   targetspin     ! spin of target
! Variables for initialization of ENDF format
!   AWR            ! standard mass parameter
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
!   ZA             ! standard charge parameter
! Variables for total cross sections in ENDF format
!   Rprime         ! potential scattering radius
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
! Variables for URR in ENDF format
!   Djlurr         ! average level spacing for resonances for j, l value
!   Eurr           ! incident energy with URR data
!   GFjlurr        ! average fission width for j, l value
!   GGjlurr        ! average radiative width for j, l value
!   GNjlurr        ! average reduced neutron width for j, l value
!   GXjlurr        ! average fission width for j, l value
!   Jurr           ! spin value
!   NJSurr         ! number of j - values
!   NLSurr         ! number of l - values
!   Nurr           ! number of incident energies with URR data
! Variables for MF1
!   EMAX           ! upper limit of energy range for evaluation
! Variables for MF2
!   ABN            ! abundance
!   AJU            ! spin
!   AMUF           ! number of degrees of freedom for fission width distribution
!   AMUG           ! number of degrees of freedom for gamma width distribution
!   AMUN           ! number of degrees of freedom for neutron width distribution
!   AMUX           ! number of degrees of freedom for competitive width distribution
!   AP             ! scattering radius
!   AWRI           ! ratio of isotope mass to neutron
!   D              ! average level spacing for resonances with spin J
!   EH             ! boundary for resonance range
!   EL             ! boundary for resonance range
!   Es             ! energy of energy - dependent width
!   GF             ! fission width of the resonance
!   GFA            ! first partial fission width
!   GFB            ! second partial fission width
!   GFu            ! average fission width
!   GG             ! gamma width of the resonance
!   GGu            ! average radiation width
!   GN0            ! average reduced neutron width
!   GX             ! average competitive reaction width
!   kINT           ! interpolation scheme
!   LFW            ! flag for average fission width
!   Lres           ! l - value
!   LRF            ! representation indicator
!   LRU            ! flag for resolved / unresolved
!   LSSF           ! flag for interpretation of MF3 cross sections
!   NAPS           ! flag for channel radius and scattering radius
!   NER            ! number of resonance energy ranges
!   NEu            ! number of energy points
!   NIS            ! number of isotopes in the material
!   NJS            ! number of J - states for a particular l - value
!   NLS            ! number of l - values
!   NRO            ! flag for energy dependence of scattering radius
!   SPI            ! target spin
!   ZAI            ! (Z, A) designation
!
! *** Declaration of local data
!
  implicit none
  integer   :: j                 ! counter
  integer   :: l                 ! counter
  integer   :: MF                ! MF-number
  integer   :: MT                ! MT-number
  integer   :: n                 ! counter
  integer   :: nen               ! energy counter
  integer   :: nenurr            ! energy counter for URR
  real(sgl) :: Ebound            ! boundary energy
  real(sgl) :: Efac              ! help variable
  real(sgl) :: Eu                ! upper energy of the URR in MeV
  real(sgl) :: Euprev            ! upper energy of the URR in MeV
  real(sgl) :: rad               ! scattering radius
!
! ******** Use resonance parameters from existing data library *********
!
! read2: subroutine to read MF2 from existing ENDF-6 data library
!
  MF = 2
  MT = 151
  if (adopt(MF, MT)) then
    call read2
    EHres = EH(NER)
!
! ************* If urrmode 1: Adopt URR parameters from TALYS **********
!
    if (lssfinp /=  -1) then
      do n = 1, NER
        LSSF(n) = lssfinp
      enddo
    endif
    if (urrmode == 1) then
      if (LRU(NER) /= 2) NER = NER + 1
      if (lssfinp /=  -1) LSSF(NER) = lssfinp
      LRU(NER) = 2
      LRF(NER) = 2
      NRO(NER) = 0
      NAPS(NER) = 0
      SPI(NER) = targetspin
      AP(NER) = Rprime * 0.1
      LSSF(NER) = 1
      NLS(NER) = NLSurr
      EL(NER) = EH(NER - 1)
      Eu = urrenergy * 1.e6
      if (EL(NER) >= 0.8 * Eu) urrenergy = 1.5 * EL(NER) * 1.e-6
      EH(NER) = urrenergy * 1.e6
      nenurr = 0
      do nen = 1, Nurr
        if (Eurr(nen) > urrenergy) then
          nenurr = nen
          exit
        endif
      enddo
      do l = 1, NLS(NER)
        AWRI(NER, l) = AWR
        Lres(NER, l) = l - 1
        NJS(NER, l) = NJSurr(l)
 Loop1: do j = 1, NJS(NER, l)
          AJU(NER, l, j) = Jurr(l, j)
          kINT(NER, l, j) = 5
          AMUX(NER, l, j) = 0.
          AMUN(NER, l, j) = 1.
          AMUG(NER, l, j) = 0.
          if (flagfission) then
            AMUF(NER, l, j) = 1.
          else
            AMUF(NER, l, j) = 0.
          endif
          n = 0
          Euprev = 1.e-10
          do nen = 1, nenurr
            Eu = Eurr(nen) * 1.e6
            if ((Eu >= EL(NER) .and. Euprev < EL(NER)) .or. nen == nenurr) then
              if (nen == nenurr) then
                Ebound = EH(NER)
                Efac = (Ebound - Euprev) / (Eu - Euprev)
              else
                Ebound = EL(NER)
                Efac = 1.
              endif
              n = n + 1
              Es(NER, l, j, n) = Ebound
              D(NER, l, j, n) = Djlurr(nen - 1, l, j) + Efac * (Djlurr(nen, l, j) - Djlurr(nen - 1, l, j))
              GX(NER, l, j, n) = GXjlurr(nen - 1, l, j) + Efac * (GXjlurr(nen, l, j) - GXjlurr(nen - 1, l, j))
              GN0(NER, l, j, n) = GNjlurr(nen - 1, l, j) + Efac * (GNjlurr(nen, l, j) - GNjlurr(nen - 1, l, j))
              GGu(NER, l, j, n) = GGjlurr(nen - 1, l, j) + Efac * (GGjlurr(nen, l, j) - GGjlurr(nen - 1, l, j))
              GFu(NER, l, j, n) = GFjlurr(nen - 1, l, j) + Efac * (GFjlurr(nen, l, j) - GFjlurr(nen - 1, l, j))
            endif
            if (Eu > EL(NER) .and. nen < nenurr) then
              n = n + 1
              Es(NER, l, j, n) = Eu
              D(NER, l, j, n) = Djlurr(nen, l, j)
              GX(NER, l, j, n) = GXjlurr(nen, l, j)
              GN0(NER, l, j, n) = GNjlurr(nen, l, j)
              GGu(NER, l, j, n) = GGjlurr(nen, l, j)
              GFu(NER, l, j, n) = GFjlurr(nen, l, j)
            endif
            if (n == 0) cycle Loop1
            GX(NER, l, j, n) = abs(GX(NER, l, j, n))
            GN0(NER, l, j, n) = abs(GN0(NER, l, j, n))
            GGu(NER, l, j, n) = abs(GGu(NER, l, j, n))
            GFu(NER, l, j, n) = abs(GFu(NER, l, j, n))
            if (n == 40) then
              Es(NER, l, j, n) = EH(NER)
              exit
            endif
            Euprev = Eu
          enddo
          if (GX(NER, l, j, n) > 0.) then
            if (GX(NER, l, j, n - 1) > 0.) then
              AMUX(NER, l, j) = 1.
            else
              GX(NER, l, j, n) = 0.
            endif
          endif
          NEu(NER, l, j) = n
        enddo Loop1
      enddo
    endif
    if (urrmode >= 1 .and. LSSF(NER) == 1) then
      EHres = EL(NER)
    else
      EHres = EH(NER)
    endif
    if (urrmode == 0 .and. NER > 1) NER = NER - 1
  else
!
! *********** Potential scattering radius provided by TALYS ************
!
! If nothing is available, we include a trivial MF2. We use r=1.35*(A**1/3) for the scattering radius.
!
    LFW = 0
    NER = 1
    EL(1) = EmineV
    EH(1) = EMAX
    LRU(1) = 0
    LRF(1) = 0
    NRO(1) = 0
    NAPS(1) = 0
    SPI(1) = targetspin
    rad = Rprime * 0.1
    AP(1) = rad
    NLS(1) = 0
  endif
  NIS = 1
  ZAI = ZA
  ABN = 1.
  mfexist(MF) = .true.
  mtexist(MF, MT) = .true.
!
! Tiny fission widths for subactinide fission
!
  if (flagfission .and. flagsubfis .and. .not. flagfis10) then
    LFW = 1
    do n = 1, numres
      do l = 1, numlres
        do j = 1, numjres
          do nen = 1, numnrs
            if (GFu(n, l, j, nen) == 0.) GFu(n, l, j, nen) = GGu(n, l, j, nen) * 1.e-10
          enddo
        enddo
        do nen = 1, numnrs
          if (GF(n, l, nen) == 0.) GF(n, l, nen) = GG(n, l, nen) * 1.e-10
          if (GFA(n, l, nen) == 0.) GFA(n, l, nen) = GG(n, l, nen) * 1.e-10
          if (GFB(n, l, nen) == 0.) GFB(n, l, nen) = GG(n, l, nen) * 1.e-10
        enddo
      enddo
    enddo
  endif
  return
end subroutine make2
! Copyright A.J. Koning 2021
