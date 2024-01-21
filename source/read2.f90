subroutine read2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF2 from existing ENDF-6 data library
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! Variables for input of specific ENDF data
!   adoptfile     ! name of library for MF information (per MT)
!   urrmode       ! 0: no URR, 1: URR from TALYS, 2: URR from data library
! Variables for reaction initialization
!   includeres    ! flag to include resonance parameters
! Variables for TALYS info
!   Atarget       ! mass number of nucleus
!   Ztarget       ! charge number of nucleus
! Variables for ENDF format
!   INTER         ! interpolation scheme
!   NBT           ! separation value for interpolation scheme
!   NP            ! number of incident energies
!   NR            ! number of interpolation ranges
! Variables for MF2
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
!   NPP7          ! total number of particle pairs
!   NPU           ! number of URR incident energies
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
!   ZB7           ! charge of second particle in pair
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)   :: MTstr     ! string for MT number
  character(len=80)  :: string    ! line with parameter value
  character(len=132) :: afile    ! name of library for MF information (per MT)
  integer            :: Ar        ! mass number
  integer            :: i         ! counter
  integer            :: ie        ! counter
  integer            :: is        ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: istat     ! error code
  integer            :: j         ! counter
  integer            :: jj        ! counter
  integer            :: k         ! counter
  integer            :: ll        ! angular momentum
  integer            :: MF        ! MF-number
  integer            :: MT        ! MT-number
  integer            :: n         ! counter
  integer            :: nlin      ! number of lines
  integer            :: NN        ! neutron number of residual nucleus
  integer            :: nrest     ! help variable
  integer            :: Zr        ! charge number index for residual nucleus
  real(sgl)          :: rZAI      ! (Z,A) designation
!
! ****************** Read resonance parameters from MF2 ****************
!
  MF = 2
  MT = 151
  afile = adoptfile(MF, MT)
  open (unit = 3, file = afile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(afile, istat)
!
! Read until MF2 is found
!
  MTstr = '     '
  write(MTstr(1:2), '(i2)') MF
  write(MTstr(3:5), '(i3)') MT
  do
    read(3, '(a80)', iostat = istat) string
    if (istat ==  -1) then
      close(unit = 3)
      return
    endif
    if (string(71:75) == MTstr) exit
  enddo
  read(string(45:55), '(i11)') NIS
!
! Loop over isotopes
!
  LRU(1) = 0
  do is = 1, NIS
    read(3, '(e11.6, 22x, 2i11)', iostat = istat) rZAI, LFW, NER
    if (istat /= 0) call read_error(afile, istat)
    do n = 1, NER
      LSSF(n) = 0
      read(3, '(2e11.6, 4i11)', iostat = istat) EL(n), EH(n), LRU(n), LRF(n), NRO(n), NAPS(n)
      if (istat /= 0) call read_error(afile, istat)
!
! A. Potential scattering radius
!
      if (LRU(n) == 0) read(3, '(2e11.6, 22x, i11)', iostat = istat) SPI(n), AP(n), NLS(n)
      if (istat /= 0) call read_error(afile, istat)
!
! B. Resolved resonance parameters:
!
      if (LRU(n) == 1) then
!
! 1. SLBW and MLBW
!
        if (LRF(n) <= 2) then
!
! NRO=1: Energy dependent scattering radius
!
          if (NRO(n) == 1) then
            read(3, '(44x, 2i11)', iostat = istat) NR(MF, MT), NP(MF, MT)
            if (istat /= 0) call read_error(afile, istat)
            NN = NR(MF, MT)
            do i = 1, NN - 2, 3
              read(3, '(6i11)', iostat = istat) (NBT(MF, MT, j), INTER(MF, MT, j), j = i, i + 2)
            if (istat /= 0) call read_error(afile, istat)
            enddo
            nlin = NN / 3
            nrest = NN - nlin * 3
            if (nrest == 1) read(3, '(6i11)', iostat = istat) NBT(MF, MT, nlin*3+1), INTER(MF, MT, nlin * 3 + 1)
            if (istat /= 0) call read_error(afile, istat)
            if (nrest == 2) read(3, '(6i11)', iostat = istat) (NBT(MF, MT, nlin*3+j), INTER(MF, MT, nlin * 3 + j), j = 1, 2)
            if (istat /= 0) call read_error(afile, istat)
            NN = NP(MF, MT)
            do i = 1, NN - 2, 3
              read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = i, i+2)
              if (istat /= 0) call read_error(afile, istat)
            enddo
            nlin = NN / 3
            nrest = NN - nlin * 3
            if (nrest /= 0) read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = 3 * nlin + 1, 3 * nlin + 3)
            if (istat /= 0) call read_error(afile, istat)
          endif
          read(3, '(2e11.6, 22x, i11)', iostat = istat) SPI(n), AP(n), NLS(n)
          if (istat /= 0) call read_error(afile, istat)
          do ll = 1, NLS(n)
            read(3, '(2e11.6, 2i11, 11x, i11)', iostat = istat) AWRI(n, ll), QX(n, ll), Lres(n, ll), LRX(n, ll), NRS(n, ll)
            if (istat /= 0) call read_error(afile, istat)
            do ie = 1, NRS(n, ll)
              read(3, '(6e11.6)', iostat = istat) Er(n, ll, ie), AJ(n, ll, ie), &
                GT(n, ll, ie), GN(n, ll, ie), GG(n, ll, ie), GF(n, ll, ie)
              if (istat /= 0) call read_error(afile, istat)
              GT(n, ll, ie) = abs(GT(n, ll, ie))
              GN(n, ll, ie) = abs(GN(n, ll, ie))
              GG(n, ll, ie) = abs(GG(n, ll, ie))
              GF(n, ll, ie) = abs(GF(n, ll, ie))
            enddo
          enddo
        endif
!
! 2. Reich-Moore
!
        if (LRF(n) == 3) then
          if (NRO(n) == 1) then
            read(3, '(44x, 2i11)', iostat = istat) NR(MF, MT), NP(MF, MT)
            if (istat /= 0) call read_error(afile, istat)
            NN = NR(MF, MT)
            do i = 1, NN - 2, 3
              read(3, '(6i11)', iostat = istat) (NBT(MF, MT, j), INTER(MF, MT, j), j = i, i + 2)
              if (istat /= 0) call read_error(afile, istat)
            enddo
            nlin = NN / 3
            nrest = NN - nlin * 3
            if (nrest == 1) read(3, '(6i11)', iostat = istat) NBT(MF, MT, nlin*3+1), INTER(MF, MT, nlin * 3 + 1)
            if (istat /= 0) call read_error(afile, istat)
            if (nrest == 2) read(3, '(6i11)', iostat = istat) (NBT(MF, MT, nlin*3+j), INTER(MF, MT, nlin * 3 + j), j = 1, 2)
            if (istat /= 0) call read_error(afile, istat)
            NN = NP(MF, MT)
            do i = 1, NN - 2, 3
              read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = i, i+2)
              if (istat /= 0) call read_error(afile, istat)
            enddo
            nlin = NN / 3
            nrest = NN - nlin * 3
            if (nrest /= 0) read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = 3 * nlin + 1, 3 * nlin + 3)
            if (istat /= 0) call read_error(afile, istat)
          endif
          read(3, '(2e11.6, i11, 11x, 2i11)', iostat = istat) SPI(n), AP(n), LAD(n), NLS(n), NLSC(n)
          if (istat /= 0) call read_error(afile, istat)
          do ll = 1, NLS(n)
            read(3, '(2e11.6, i11, 22x, i11)', iostat = istat) AWRI(n, ll), APL(n, ll), Lres(n, ll), NRS(n, ll)
            if (istat /= 0) call read_error(afile, istat)
            do ie = 1, NRS(n, ll)
              read(3, '(6e11.6)', iostat = istat) Er(n, ll, ie), AJ(n, ll, ie), &
                GN(n, ll, ie), GG(n, ll, ie), GFA(n, ll, ie), GFB(n, ll, ie)
              if (istat /= 0) call read_error(afile, istat)
              GT(n, ll, ie) = abs(GT(n, ll, ie))
              GN(n, ll, ie) = abs(GN(n, ll, ie))
              GG(n, ll, ie) = abs(GG(n, ll, ie))
            enddo
          enddo
        endif
!
! 3. R-matrix limited format
!
        if (LRF(n) == 7) then
          read(3, '(22x, 4i11)', iostat = istat) IFG7, KRM7, NJS7, KRL7
          if (istat /= 0) call read_error(afile, istat)
          read(3, '(22x, i11)', iostat = istat) NPP7
          if (istat /= 0) call read_error(afile, istat)
          do i = 1, NPP7
            read(3, '(6e11.6)', iostat = istat) MA7(i), MB7(i), ZA7(i), ZB7(i), IA7(i), IB7(i)
            if (istat /= 0) call read_error(afile, istat)
            read(3, '(6e11.6)', iostat = istat) QP7(i), PNT7(i), SHF7(i), MTP7(i), PA7(i), PB7(i)
            if (istat /= 0) call read_error(afile, istat)
            if (ABS(QP7(i)) <= 1.e-20) QP7(i) = 0.
          enddo
          do j = 1, NJS7
            read(3, '(2e11.6, 2i11, 11x, i11)', iostat = istat) AJ7(j), PJ7(j), KBK7(j), KPS7(j), NCH7(j)
            if (istat /= 0) call read_error(afile, istat)
            do i = 1, NCH7(j)
              read(3, '(6e11.6)', iostat = istat) PPI7(j, i), L7(j, i), SCH7(j, i), BND7(j, i), APE7(j, i), APT7(j, i)
              if (istat /= 0) call read_error(afile, istat)
            enddo
            read(3, '(33x, i11, 11x, i11)', iostat = istat) NRS7(j), NX7(j)
            if (istat /= 0) call read_error(afile, istat)
            NN = NCH7(j) / 6
            do ie = 1, NRS7(j)
              read(3, '(6e11.6)', iostat = istat) ER7(j, ie), (GAM7(j, ie, k), k = 1, 5)
              if (istat /= 0) call read_error(afile, istat)
              do jj = 1, NN
                read(3, '(6e11.6)', iostat = istat) (GAM7(j, ie, 6*jj-1+k), k = 1, 6)
                if (istat /= 0) call read_error(afile, istat)
              enddo
            enddo
          enddo
        endif
      endif
!
! C. Unresolved resonance parameters:
!
      if (LRU(n) == 2) then
        if (urrmode == 2) then
!
! 1. Energy independent parameters
!
          if (LRF(n) == 1) then
            read(3, '(2e11.6, i11, 11x, i11)', iostat = istat) SPI(n), AP(n), LSSF(n), NLS(n)
            if (istat /= 0) call read_error(afile, istat)
            do ll = 1, NLS(n)
              read(3, '(e11.6, 11x, i11, 22x, i11)', iostat = istat) AWRI(n, ll), Lres(n, ll), NJS(n, ll)
              if (istat /= 0) call read_error(afile, istat)
              do jj = 1, NJS(n, ll)
                read(3, '(5e11.6)', iostat = istat) D(n, ll, jj, 1), AJU(n, ll, jj), &
                  AMUN(n, ll, jj), GN0(n, ll, jj, 1), GGu(n, ll, jj, 1)
                if (istat /= 0) call read_error(afile, istat)
              enddo
            enddo
          endif
!
! 2. Energy dependent parameters
!
          if (LRF(n) == 2) then
            if (NRO(n) == 1) then
              read(3, '(44x, 2i11)', iostat = istat) NRU, NPU
              if (istat /= 0) call read_error(afile, istat)
              NN = NRU
              read(3, '(2i11)') NBTU, INTERU
              NN = NPU
              do i = 1, NN - 2, 3
                read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = i, i+2)
                if (istat /= 0) call read_error(afile, istat)
              enddo
              nlin = NN / 3
              nrest = NN - nlin * 3
              if (nrest /= 0) read(3, '(6e11.6)', iostat = istat) (E2(n, j), APE(n, j), j = 3 * nlin + 1, 3 * nlin + 3)
              if (istat /= 0) call read_error(afile, istat)
            endif
            read(3, '(2e11.6, i11, 11x, i11)', iostat = istat) SPI(n), AP(n), LSSF(n), NLS(n)
            if (istat /= 0) call read_error(afile, istat)
            do ll = 1, NLS(n)
              read(3, '(e11.6, 11x, i11, 11x, i11)', iostat = istat) AWRI(n, ll), Lres(n, ll), NJS(n, ll)
              if (istat /= 0) call read_error(afile, istat)
              do jj = 1, NJS(n, ll)
                read(3, '(e11.6, 11x, i11, 22x, i11)', iostat = istat) AJU(n, ll, jj), kINT(n, ll, jj), NEu(n, ll, jj)
                if (istat /= 0) call read_error(afile, istat)
                read(3, '(22x, 4e11.6)', iostat = istat) AMUX(n, ll, jj), AMUN(n, ll, jj), AMUG(n, ll, jj), AMUF(n, ll, jj)
                if (istat /= 0) call read_error(afile, istat)
                do ie = 1, NEu(n, ll, jj)
                  read(3, '(6e11.6)', iostat = istat) Es(n, ll, jj, ie), D(n, ll, jj, ie), &
                    GX(n, ll, jj, ie), GN0(n, ll, jj, ie), GGu(n, ll, jj, ie), GFu(n, ll, jj, ie)
                  if (istat /= 0) call read_error(afile, istat)
                enddo
              enddo
            enddo
          endif
        endif
!
! Set urrmode=0 if no URR is given by either TALYS or other data library
!
      else
        if (urrmode == 2 .and. n == NER) urrmode = 0
      endif
    enddo
!
! Jump out of the loop if this was the correct isotope
!
    Zr = int(rZAI) / 1000
    Ar = mod(int(rZAI), 1000)
    if (Zr == Ztarget .and. Ar == Atarget) exit
  enddo
  close (unit = 3)
  if (LRU(1) /= 0) includeres = .true.
  return
end subroutine read2
! Copyright A.J. Koning 2021
