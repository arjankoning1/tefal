subroutine write2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF2
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
! All global variables
!   numenin       ! number of incident energies
!   numint        ! number of interpolation sections
! Variables for initialization of ENDF format
!   AWR           ! standard mass parameter
!   blank2        ! blank string
!   FEND          ! ENDF - 6 format
!   MAT           ! MAT number
!   SEND          ! ENDF - 6 format
!   ZA            ! standard charge parameter
! Variables for ENDF format
!   INTER         ! interpolation scheme
!   NBT           ! separation value for interpolation scheme
!   NP            ! number of incident energies
!   NR            ! number of interpolation ranges
! Variables for MF2
!   ABN           ! abundance
!   AJ            ! spin of the resonance
!   AJ7           ! spin
!   AJU           ! spin
!   AMUF          ! number of degrees of freedom for fission width distribution
!   AMUG          ! number of degrees of freedom for gamma width distribution
!   AMUN          ! number of degrees of freedom for neutron width distribution
!   AMUX          ! number of degrees of freedom for competitive width distribution
!   AP            ! scattering radius
!   APE           ! scattering radius
!   APE7          ! effective channel radius
!   APL           ! l - dependent scattering radius
!   APT7          ! true channel radius
!   AWRI          ! ratio of isotope mass to neutron
!   BND7          ! boundary condition for this channel
!   D             ! average level spacing for resonances with spin J
!   E2            ! incident energy for MF2 (in ENDF - 6 format)
!   EH            ! boundary for resonance range
!   EL            ! boundary for resonance range
!   Er            ! resonance energy in LAB system
!   ER7           ! energy of resonance in eV
!   Es            ! energy of energy - dependent width
!   GAM7          ! channel width or reduced width
!   GF            ! fission width of the resonance
!   GFA           ! first partial fission width
!   GFB           ! second partial fission width
!   GFu           ! average fission width
!   GG            ! gamma width of the resonance
!   GGu           ! average radiation width
!   GN            ! neutron width of the resonance
!   GN0           ! average reduced neutron width
!   GT            ! total width of the resonance
!   GX            ! average competitive reaction width
!   IA7           ! spin of first particle in pair
!   IB7           ! spin of second particle in pair
!   IFG7          ! flag for gamma width
!   INTERU        ! URR interpolation scheme
!   KBK7          ! R - matrix background parameter
!   kINT          ! interpolation scheme
!   KPS7          ! non - hard - sphere phase shift parameter
!   KRL7          ! flag for relativistic kinematics
!   KRM7          ! flag of formula to be used (MLBW, RM, etc.) for R - matrix
!   L7            ! orbital angular momentum
!   LAD           ! flag to indicate computation of angular distributions
!   LFW           ! flag for average fission width
!   Lres          ! l - value
!   LRF           ! representation indicator
!   LRU           ! flag for resolved / unresolved
!   LRX           ! flag to indicate competitive width
!   LSSF          ! flag for interpretation of MF3 cross sections
!   MA7           ! mass of first particle in pair
!   MB7           ! mass of second particle in pair
!   MTP7          ! reaction type for this particle pair
!   NAPS          ! flag for channel radius and scattering radius
!   NBTU          ! separation value for URR interpolation scheme
!   NCH7          ! number of channels per J
!   NER           ! number of resonance energy ranges
!   NEu           ! number of energy points
!   NIS           ! number of isotopes in the material
!   NJS           ! number of J - states for a particular l - value
!   NJS7          ! number of J - states for a particular l - value
!   NLS           ! number of l - values
!   NLSC          ! number of l - values for convergence of angular distributions
!   NPU           ! number of URR incident energies
!   NPP7          ! total number of particle pairs
!   NRO           ! flag for energy dependence of scattering radius
!   NRS           ! number of resolved resonances per l - value
!   NRS7          ! number of resonances
!   NRU           ! number of URR interpolation ranges
!   NX7           ! number of lines required for all resonances
!   PA7           ! parity of first particle in pair
!   PB7           ! parity of second particle in pair
!   PJ7           ! parity
!   PNT7          ! flag for penetrability
!   PPI7          ! particle - pair number for this channel
!   QP7           ! Q - value for this particle pair
!   QX            ! Q - value to be added to C.M. incident energy
!   SCH7          ! channel spin
!   SHF7          ! flag for shift factor
!   SPI           ! target spin
!   ZA7           ! charge of first particle in pair
!   ZAI           ! (Z, A) designation
!   ZB7           ! charge of second particle in pair
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: ie                        ! counter
  integer   :: ii                        ! counter
  integer   :: j                         ! counter
  integer   :: jj                        ! counter
  integer   :: k                         ! counter
  integer   :: ll                        ! angular momentum
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: n                         ! counter
  integer   :: NN                        ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(2*numenin)              ! help variable
!
! ***************************** Write MF2 ******************************
!
! hrwrite: subroutine to write header with real values
!
! ********* Use resonance parameters from existing data library ********
!
  MF = 2
  MT = 151
  NS = 0
  open (unit = 2, file = 'MF2', status = 'replace')
  call hrwrite(ZA, AWR, 0, 0, NIS, 0, MAT, MF, MT, NS)
  call hrwrite(ZAI, ABN, 0, LFW, NER, 0, MAT, MF, MT, NS)
  do n = 1, NER
    call hrwrite(EL(n), EH(n), LRU(n), LRF(n), NRO(n), NAPS(n), MAT, MF, MT, NS)
!
! A. Potential scattering radius
!
    if (LRU(n) == 0) call hrwrite(SPI(n), AP(n), 0, 0, NLS(n), 0, MAT, MF, MT, NS)
!
! B. Resolved resonance parameters
!
    if (LRU(n) == 1) then
!
! 1. SLBW and MLBW
!
      if (LRF(n) <= 2) then
        if (NRO(n) == 1) then
!
! Energy dependent scattering radius
!
! Write interpolation ranges
!
! kwrite: subroutine to write integer value block
! xwrite: subroutine to write real value block
!
          call hrwrite(0., 0., 0, 0, NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
          do i = 1, NR(MF, MT)
            ii = 2 * i - 1
            Nval(ii) = NBT(MF, MT, i)
            Nval(ii + 1) = INTER(MF, MT, i)
          enddo
          NN = 2 * NR(MF, MT)
          call kwrite(NN, Nval, MAT, MF, MT, NS)
!
! Write values
!
          do i = 1, NP(MF, MT)
            ii = 2 * i - 1
            x(ii) = E2(n, i)
            x(ii + 1) = APE(n, i)
          enddo
          NN = 2 * NP(MF, MT)
          call xwrite(NN, x, MAT, MF, MT, NS)
        endif
        call hrwrite(SPI(n), AP(n), 0, 0, NLS(n), 0, MAT, MF, MT, NS)
        do ll = 1, NLS(n)
          call hrwrite(AWRI(n, ll), QX(n, ll), Lres(n, ll), LRX(n, ll), 6 * NRS(n, ll), NRS(n, ll), MAT, MF, MT, NS)
          do ie = 1, NRS(n, ll)
            call rwrite(Er(n, ll, ie), AJ(n, ll, ie), GT(n, ll, ie), GN(n, ll, ie), GG(n, ll, ie), GF(n, ll, ie), MAT, MF, MT, NS)
          enddo
        enddo
      endif
!
! 2. Reich-Moore
!
      if (LRF(n) == 3) then
        if (NRO(n) == 1) then
          call hrwrite(0., 0., 0, 0, NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
          do i = 1, NR(MF, MT)
            ii = 2 * i - 1
            Nval(ii) = NBT(MF, MT, i)
            Nval(ii + 1) = INTER(MF, MT, i)
          enddo
          NN = 2 * NR(MF, MT)
          call kwrite(NN, Nval, MAT, MF, MT, NS)
!
! Write values
!
          do i = 1, NP(MF, MT)
            ii = 2 * i - 1
            x(ii) = E2(n, i)
            x(ii + 1) = APE(n, i)
          enddo
          NN = 2 * NP(MF, MT)
          call xwrite(NN, x, MAT, MF, MT, NS)
        endif
        call hrwrite(SPI(n), AP(n), LAD(n), 0, NLS(n), NLSC(n), MAT, MF, MT, NS)
        do ll = 1, NLS(n)
          call hrwrite(AWRI(n, ll), APL(n, ll), Lres(n, ll), 0, 6 * NRS(n, ll), NRS(n, ll), MAT, MF, MT, NS)
          do ie = 1, NRS(n, ll)
            call rwrite(Er(n, ll, ie), AJ(n, ll, ie), GN(n, ll, ie), GG(n, ll, ie), GFA(n, ll, ie), GFB(n, ll, ie), MAT, MF, MT, NS)
          enddo
        enddo
      endif
!
! 3. R-matrix limited format
!
      if (LRF(n) == 7) then
        call hrwrite(0., 0., IFG7, KRM7, NJS7, KRL7, MAT, MF, MT, NS)
        call hrwrite(0., 0., NPP7, 0, 12 * NPP7, 2 * NPP7, MAT, MF, MT, NS)
        do i = 1, NPP7
          call rwrite(MA7(i), MB7(i), ZA7(i), ZB7(i), IA7(i), IB7(i), MAT, MF, MT, NS)
          call rwrite(QP7(i), PNT7(i), SHF7(i), MTP7(i), PA7(i), PB7(i), MAT, MF, MT, NS)
        enddo
        do j = 1, NJS7
          call hrwrite(AJ7(j), PJ7(j), KBK7(j), KPS7(j), 6 * NCH7(j), NCH7(j), MAT, MF, MT, NS)
          do i = 1, NCH7(j)
             call rwrite(PPI7(j, i), L7(j, i), SCH7(j, i), BND7(j, i), APE7(j, i), APT7(j, i), MAT, MF, MT, NS)
          enddo
          call hrwrite(0., 0., 0, NRS7(j), 6 * NX7(j), NX7(j), MAT, MF, MT, NS)
          NN = 1 + NCH7(j)
          do ie = 1, NRS7(j)
            x(1) = ER7(j, ie)
            do k = 1, NCH7(j)
              x(k + 1) = GAM7(j, ie, k)
            enddo
            call xwrite(NN, x, MAT, MF, MT, NS)
          enddo
        enddo
      endif
    endif
!
! C. Unresolved resonance parameters
!
    if (LRU(n) == 2) then
!
! 1. Energy independent parameters
!
      if (LRF(n) == 1) then
        call hrwrite(SPI(n), AP(n), LSSF(n), 0, NLS(n), 0, MAT, MF, MT, NS)
        do ll = 1, NLS(n)
          call hrwrite(AWRI(n, ll), 0., Lres(n, ll), 0, 6 * NJS(n, ll), NJS(n, ll), MAT, MF, MT, NS)
          do jj = 1, NJS(n, ll)
            call rwrite(D(n, ll, jj, 1), AJU(n, ll, jj), AMUN(n, ll, jj), GN0(n, ll, jj, 1), GGu(n, ll, jj, 1), 0., MAT, MF, MT, NS)
          enddo
        enddo
      endif
!
! 2. Energy dependent parameters
!
      if (LRF(n) == 2) then
        if (NRO(n) == 1) then
          call hrwrite(0., 0., 0, 0, NRU, NPU, MAT, MF, MT, NS)
          do i = 1, NRU
            ii = 2 * i - 1
            Nval(ii) = NBTU
            Nval(ii + 1) = INTERU
          enddo
          NN = 2 * NRU
          call kwrite(NN, Nval, MAT, MF, MT, NS)
!
! Write values
!
          do i = 1, NPU
            ii = 2 * i - 1
            x(ii) = E2(n, i)
            x(ii + 1) = APE(n, i)
          enddo
          NN = 2 * NPU
          call xwrite(NN, x, MAT, MF, MT, NS)
        endif
        call hrwrite(SPI(n), AP(n), LSSF(n), 0, NLS(n), 0, MAT, MF, MT, NS)
        do ll = 1, NLS(n)
          call hrwrite(AWRI(n, ll), 0., Lres(n, ll), 0, NJS(n, ll), 0, MAT, MF, MT, NS)
          do jj = 1, NJS(n, ll)
            call hrwrite(AJU(n, ll, jj), 0., kINT(n, ll, jj), 0, 6 * NEu(n, ll, jj) + 6, NEu(n, ll, jj), MAT, MF, MT, NS)
            call rwrite(0., 0., AMUX(n, ll, jj), AMUN(n, ll, jj), AMUG(n, ll, jj), AMUF(n, ll, jj), MAT, MF, MT, NS)
            do ie = 1, NEu(n, ll, jj)
              call rwrite(Es(n, ll, jj, ie), D(n, ll, jj, ie), &
 &              GX(n, ll, jj, ie), GN0(n, ll, jj, ie), GGu(n, ll, jj, ie), GFu(n, ll, jj, ie), MAT, MF, MT, NS)
            enddo
          enddo
        enddo
      endif
    endif
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write2
! Copyright A.J. Koning 2021
