subroutine make33(MF)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF33 and MF40
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
!   numchan        ! maximum number of exclusive channels
!   numchancov     ! number of channels with inter - channel correlations
!   numencov       ! number of incident energies for covariances
!   numencovtot    ! number of energies for covariances
!   numenin        ! number of incident energies
!   nummf          ! number of MF numbers
!   nummt          ! number of MT numbers
! Variables for input of ENDF structure
!   flagsubfis     ! flag to include subactinide fission
! Variables for input of specific ENDF data
!   adopt          ! logical for existence of MF information (per MT)
! Variables for input of ENDF library type
!   flaggpf        ! flag for general purpose library
!   flageaf        ! flag for EAF - formatted activation library
! Variables for ENDF limits, switches and tolerances
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
! Variables for ENDF covariance input
!   covdiscrete    ! number of disc. inelastic levels with covariances
!   flagparcov     ! flag to include covariances for MT600 - 849
! Variables for info from TALYS
!   Atarget        ! mass number of nucleus
!   k0             ! index of incident particle
!   Ztarget        ! charge number of nucleus
! Variables for reaction initialization
!   EHres          ! upper energy in resonance range
! Variables for initialization of ENDF format
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
! Variables for covariances in ENDF format
!   Ecov           ! energy grid for covariances
!   flagMTint      ! flag for channel with inter - MT covariance data
!   MTindex        ! index for MT number
!   MTindexiso     ! index for isomer of MT number
!   MTintindex     ! index for MT number with inter - MT covariance data
!   Nchancov       ! number of channels with covariance data
!   Ncovrp         ! neutron number of residual product
!   Nencov         ! number of energies with covariance data
!   Nisocov        ! number of isomers with covariances per MT number
!   Rcov           ! relative covariance matrix
!   relerr         ! relative cross section uncertainty
!   Rmt            ! relative covariance matrix within same MT number
!   Rrp            ! covariance element for residual cross section
!   xserr          ! cross section uncertainty
! Variables for MF2
!   LRU            ! flag for resolved / unresolved
!   NER            ! number of resonance energy ranges
! Variables for MF3
!   EthMT          ! threshold energy
! Variables for MF8_10
!   EthZA          ! threshold energy
!   NSt            ! number of final states
! Variables for MF31_40
!   b33            ! covariance matrix element
!   b33MT          ! covariance matrix element
!   b33MTread      ! covariance matrix element
!   b33read      ! covariance matrix element
!   b33ZA          ! covariance matrix element
!   b8             ! covariance matrix element
!   b8read         ! covariance matrix element
!   Emincov        ! minimum energy for covariance information
!   IZAP1          ! second IZAP - number
!   LB             ! flag for meaning of numbers
!   LBread       ! flag for meaning of numbers
!   LB8            ! flag for meaning of numbers
!   LB8read        ! flag for meaning of numbers
!   LBZA           ! flag for meaning of numbers
!   LS             ! symmetry flag
!   LSread       ! symmetry flag
!   LSZA           ! symmetry flag
!   MAT1           ! second MAT - number
!   MAT1read     ! second MAT-number
!   MAT1ZA         ! second MAT - number
!   MT33           ! second MT - number
!   MT33read     ! second MT-number
!   MT33ZA         ! second MT - number
!   MTL            ! lumped reaction identifier
!   MTLread      ! lumped reaction identifier
!   NC33read     ! number of NC-type sub-subsections
!   NC33           ! number of NC - type sub - subsections
!   NC33ZA         ! number of NC - type sub - subsections
!   NE8            ! number of entries for LB = 8 section
!   NE8read        ! number of entries for LB=8 section
!   NE33           ! number of energies in energy array
!   NE33read       ! number of energies in energy array
!   NE33ZA         ! number of energies in energy array
!   NI33           ! number of NI - type sub - subsections
!   NI33read     ! number of NI-type sub-subsections
!   NI33ZA         ! number of NI - type sub - subsections
!   NL33           ! number of subsections
!   NL33read     ! number of subsections
!   NT8read      ! total number of entries
!   NT33           ! total number of entries
!   NT33read     ! total number of entries
!   NT33ZA         ! total number of entries
!   NT8            ! total number of entries
!   XLFS1          ! second state discrete level number
!   XMF1           ! second MF - number
!   XLFS1read    ! second state discrete level number
!   XMF1read     ! second MF-number
!
! *** Declaration of local data
!
  implicit none
  logical   :: flag33(nummf, nummt)                                          ! flag for existence of covariance information per MT
  logical   :: intercov(numchan, numchan)                                    ! flag for inter-channel covariances
  integer   :: covstep                                                       ! energy step for covariances
  integer   :: i                                                             ! counter
  integer   :: i1                                                            ! value
  integer   :: i2                                                            ! value
  integer   :: ib                                                            ! counter
  integer   :: ichan                                                         ! counter for channels
  integer   :: ichan2                                                        ! counter for channels
  integer   :: icov(numenin)                                                 ! index for covariance
  integer   :: iE                                                            ! energy counter
  integer   :: iE5                                                           ! energy counter
  integer   :: iE5b                                                          ! energy counter
  integer   :: iE6                                                           ! energy counter
  integer   :: iE6b                                                          ! energy counter
  integer   :: iE8                                                           ! energy counter
  integer   :: iEb                                                           ! energy counter
  integer   :: iza                                                           ! counter for Z,A combinations
  integer   :: j                                                             ! counter
  integer   :: MF                                                            ! MF-number
  integer   :: MT                                                            ! MT-number
  integer   :: MT2                                                           ! MT number
  integer   :: N33read                                                       ! number of energies in energy array
  integer   :: NE2cov(numchan)                                               ! number of energies
  integer   :: nen                                                           ! energy counter
  integer   :: Nsec                                                          ! number of sections
  integer   :: Nthresh(nummt)                                                ! energy point at threshold
  integer   :: NthZA                                                         ! counter
  real(sgl) :: E33read                                                       ! energy of covariance grid
  real(sgl) :: ELB5(0:numencov+1)                                            ! energy of covariance grid
  real(sgl) :: ELB6(0:numencov+1)                                            ! energy of covariance grid
  real(sgl) :: ELB6b(numchan, 0:numencov+1)                                  ! energy of covariance grid
  real(sgl) :: ELB8(0:numencov+1)                                            ! energy of covariance grid
  real(sgl) :: err                                                           ! error
  real(sgl) :: EZA(0:numencov+1)                                             ! energy of covariance grid
  real(sgl) :: Rfinal(numchancov, numencov, numchancov, numencov)            ! relative covariance matrix
  real(sgl) :: Rmt8(numchan, numencovtot)                                    ! relative covariance matrix within same MT number
  real(sgl) :: Rmtfinal(numchan, numencov, numencov)                         ! relative covariance matrix within same MT number
  real(sgl) :: RZA(numchan, numencov, numencov)                              ! realtive covariance matrix for residual cross section
!
! ***************************** Make MF33 ******************************
!
! read35     : subroutine to read MF35 from existing ENDF-6 data library
! locate     : subroutine to find value in ordered table
!
! 1. Covariance data
!
! Determine minimum energy from which covariance information is given
!
  covstep = 2
  intercov = .false.
  Emincov = 0.
  if (flaggpf .and. k0 == 1) then
    if (LRU(NER) > 0) Emincov = EHres * 1.e-6
  endif
  NL33read = 0
  NC33read = 0
  NI33read = 0
  NT33read = 0
  NE33read = 0
  NE2cov = 0
  LB8read = 0
  NE8read = 0
  NT8read = 0
  MT33read = 0
  MTLread = 0
  MAT1read = 0
  LSread = 0
  LBread = 0
  b33read = 0.
  b33MTread = 0.
  RZA = 0.
  Rmtfinal = 0.
  Rmt8 = 0.
  Rfinal = 0.
  b8read = 0.
  XMF1read = 0.
  XLFS1read = 0.
  b33 = 0.
  b33MT = 0.
  b8 = 0.
  b33ZA = 0.
  EZA = 0.
  ELB5 = 0.
  ELB6 = 0.
  ELB8 = 0.
  icov = 1
Loop1:  do MT = 1, nummt
    Nthresh(MT) = 1000
    if ( .not. mtexist(MF - 30, MT)) cycle
!
! Determine threshold energy
!
    do i = 2, Nencov
      if (Ecov(i - 1) < EthMT(MT) .and. Ecov(i) >= EthMT(MT)) then
        Nthresh(MT) = i
        cycle Loop1
      endif
    enddo
  enddo Loop1
  Ecov(Nencov + 1) = Ecov(Nencov)
  call locate(Ecov, 1, Nencov + 1, Emincov, nen)
  Ecov(nen + 1) = Emincov
  if (Ecov(nen + 1) == Ecov(nen) .and. nen >= 1) Ecov(nen + 1) = 1.001 * Ecov(nen)
  do MT = 1, nummt
    flag33(MF, MT) = .false.
    if ( .not. mtexist(MF - 30, MT)) cycle
    if (MF == 40) then
      if ( .not. mtexist(10, MT)) cycle
      if (MT <= 3) cycle
      if (MT >= 51 .and. MT <= 91) cycle
      if (MT >= 201) cycle
!
! ********** Covariances for residual production cross sections ********
!
      if (MT == 5 .and. Ncovrp > 0) then
        do iza = 1, Ncovrp
          NthZA = 1000
          do i = 2, Nencov
            if (Ecov(i - 1) < EthZA(iza) .and. Ecov(i) >= EthZA(iza)) then
              NthZA = i
              exit
            endif
          enddo
          iE = 1
          EZA(1) = Ecov(1)
          do i = 2, Nencov
            if (Ecov(i) < EthZA(iza)) cycle
!
! Energy grid
!
            if (i == NthZA) then
              iE = iE + 1
              EZA(iE) = EthZA(iza)
            else
              if (mod(i, covstep) == 1 .and. i /= Nencov) cycle
              iE = iE + 1
              EZA(iE) = Ecov(i)
            endif
!
! Intra-MT covariances
!
            iEb = iE - 1
            do j = i, nencov
              if (mod(j, covstep) == 1 .and. j /= i .and. j /= Nencov) cycle
              iEb = iEb + 1
              RZA(ichan, iE, iEb) = Rrp(iza, i, j)
            enddo
          enddo
          NE33ZA(iza, 1) = iE
!
! Covariances in ENDF format
!
! Energy grid
!
          do iE = 1, NE33ZA(iza, 1)
            b33ZA(iza, iE) = EZA(iE) * 1.e6
          enddo
!
! Intra-MT correlations
!
          ib = NE33ZA(iza, 1)
          do iE = 2, NE33ZA(iza, 1)
            do iEb = iE, NE33ZA(iza, 1)
              ib = ib + 1
              b33ZA(iza, ib) = RZA(iza, iE, iEb)
            enddo
          enddo
          b33ZA(iza, 1) = 1.e-5
          MT33ZA(iza, 1) = 5
          NC33ZA(iza, 1) = 0
          NT33ZA(iza, 1) = (NE33ZA(iza, 1) * (NE33ZA(iza, 1) + 1)) / 2
          NI33ZA(iza, 1) = 1
          MAT1ZA(iza, 1) = 0
          LSZA(iza, 1) = 1
          LBZA(iza, 1) = 5
        enddo
        flag33(MF, MT) = .true.
        NL33(MF, MT) = 1
        goto 400
      endif
    endif
!
! ********** Covariances for all other channels ************************
!
    if (MT >= 51 .and. MT <= 91 .and. covdiscrete == 0) cycle
    if (MT >= 51 .and. MT <= 90 .and. MT > 50 + covdiscrete) cycle
    if (MT >= 600 .and. .not. flagparcov) cycle
    do ichan = 1, Nchancov
      if (MT /= MTindex(ichan)) cycle
!
! ******************* Adoption of external covariance data *************
!
      E33read = 0.
      N33read = 0
      if (adopt(MF, MT) .and. MTindexiso(ichan) ==  -1) then
        call read33(MT)
        N33read = NE33read(ichan, 1)
        do iE5 = 1, N33read
          ELB5(iE5) = b33MTread(ichan, iE5) * 1.e-6
        enddo
        E33read = ELB5(N33read)
        ib = N33read
        do iE5 = 2, N33read
          do iE5b = iE5, N33read
            ib = ib + 1
            Rmtfinal(ichan, iE5, iE5b) = b33MTread(ichan, ib)
          enddo
        enddo
        iE5 = N33read
        do iE8 = 1, NE8read(ichan)
          ELB8(iE8) = b8read(ichan, 2 * iE8 - 1) * 1.e-6
          Rmt8(ichan, iE8) = b8read(ichan, 2 * iE8)
        enddo
        iE8 = NE8read(ichan)
      else
        ELB5(1) = Ecov(1)
        iE5 = 1
        ELB8(1) = Ecov(1)
        iE8 = 1
      endif
      ELB6(1) = Ecov(1)
      iE6 = 1
!
! *********************** Covariances for fast range *******************
!
      do i = 2, Nencov
        if (Ecov(i) < EthMT(MT)) cycle
        if (Ecov(i) <= 1.001 * E33read) cycle
        if (Ecov(i) > Eswitch .and. (MT >= 6 .or. MT == 4) .and. MT /= 18) cycle
!
! Energy grid
!
        if (i == Nthresh(MT)) then
          iE5 = iE5 + 1
          ELB5(iE5) = EthMT(MT)
          iE6 = iE6 + 1
          ELB6(iE6) = EthMT(MT)
          iE8 = iE8 + 1
          ELB8(iE8) = EthMT(MT)
        else
          if (mod(i, covstep) == 1 .and. i /= Nencov) cycle
          iE5 = iE5 + 1
          ELB5(iE5) = Ecov(i)
          iE6 = iE6 + 1
          ELB6(iE6) = Ecov(i)
          iE8 = iE8 + 1
          ELB8(iE8) = Ecov(i)
        endif
!
! Variances
!
        if (flageaf) then
          err = relerr(ichan, i) **2
          Rmt8(ichan, iE8) = err
        else
          err = (xserr(ichan, i) * 1.e-3) **2
          Rmt8(ichan, iE8) = 1.e-3 * err
          if (i == 2) then
            err = (xserr(ichan, 1) * 1.e-3) **2
            Rmt8(ichan, 1) = 1.e-3 * err
          endif
        endif
!
! Intra-MT covariances
!
        iE5b = iE5 - 1
        do j = i, nencov
          if (Ecov(j) <= 1.001 * E33read) cycle
          if (mod(j, covstep) == 1 .and. j /= i .and. j /= Nencov) cycle
          iE5b = iE5b + 1
          Rmtfinal(ichan, iE5, iE5b) = Rmt(ichan, i, j)
        enddo
!
! Inter-MT correlations
!
        NE2cov = 0
        if (MF == 33 .and. flagMTint(ichan)) then
          i1 = MTintindex(ichan)
          do MT2 = MT + 1, nummt
            do ichan2 = 1, Nchancov
              if (MT2 /= MTindex(ichan2)) cycle
              if ( .not. flagMTint(ichan2)) cycle
              if ( .not. mtexist(MF - 30, MT2)) cycle
              if (MT2 >= 51 .and. MT2 <= 91 .and. covdiscrete == 0) cycle
              if (MT2 >= 51 .and. MT2 <= 90 .and. MT2 > 50 + covdiscrete) cycle
              if (flagsubfis .and. MT2 == 18) cycle
              i2 = MTintindex(ichan2)
              intercov(ichan, ichan2) = .true.
              ELB6b(ichan2, 1) = Ecov(1)
              iE6b = 1
              do j = 2, Nencov
                if (Ecov(j) < EthMT(MT2)) cycle
                if (Ecov(j) > Eswitch .and. (MT2 >= 6 .or. MT2 == 4) .and. MT /= 18) cycle
                if (j == Nthresh(MT2)) then
                  iE6b = iE6b + 1
                  ELB6b(ichan2, iE6b) = Ecov(j)
                else
                  if (mod(j, covstep) == 1 .and. j /= Nencov) cycle
                  iE6b = iE6b + 1
                  ELB6b(ichan2, iE6b) = Ecov(j)
                endif
                if (i > 1 .and. j > 1) then
!   if (Ecov(i - 1) >= Emincov .and. Ecov(j - 1) >= Emincov)
                  Rfinal(i1, iE6, i2, iE6b) = Rcov(i1, i, i2, j)
                  if (Ecov(j - 1) < EthMT(MT2) .and. Ecov(j) >= EthMT(MT2)) Rfinal(i1, iE6, i2, iE6b) = 0.
                endif
              enddo
              NE2cov(ichan2) = iE6b
            enddo
          enddo
        endif
      enddo
      NE33(ichan, 1) = iE5
      NE33(ichan, 2) = iE6
      NE8(ichan) = iE8
!
! *********************** Covariances in ENDF format *******************
!
! Energy grid
!
      do iE5 = 1, NE33(ichan, 1)
        b33MT(ichan, iE5) = ELB5(iE5) * 1.e6
      enddo
!
! Variances
!
      do iE8 = 1, NE8(ichan)
        b8(ichan, 2 * iE8 - 1) = ELB8(iE8) * 1.e6
        b8(ichan, 2 * iE8) = max(Rmt8(ichan, iE8), 1.e-20)
      enddo
!
! Set last element to zero (Caleb Mattoon)
!
      b8(ichan, 2 * NE8(ichan)) = 0.
!
! Intra-MT correlations
!
      ib = NE33(ichan, 1)
      if ( .not. flageaf) then
        do iE5 = 2, NE33(ichan, 1)
          do iE5b = iE5, NE33(ichan, 1)
            ib = ib + 1
            b33MT(ichan, ib) = Rmtfinal(ichan, iE5, iE5b)
          enddo
        enddo
      endif
      b33MT(ichan, 1) = 1.e-5
      if ((b33MT(ichan, 1) == b33MT(ichan, 2)) .and. (b33MT(ichan, 2) /= b33MT(ichan, 3))) &
 &      b33MT(ichan, 2) = 0.5 * (b33MT(ichan, 1) + b33MT(ichan, 3))
!
! ENDF parameters
!
      if (MF == 33) then
        flag33(MF, MT) = .true.
      else
        if (MTindexiso(ichan) /=  -1) flag33(MF, MT) = .true.
        if (Nisocov(MT) /= NSt(MT)) flag33(MF, MT) = .false.
      endif
      if (mtexist(10, 18)) flag33(40, MT) = .true.
      MT33(ichan, 1) = MTindex(ichan)
      NC33(ichan, 1) = 0
      NT33(ichan, 1) = (NE33(ichan, 1) * (NE33(ichan, 1) + 1)) / 2
      NI33(ichan, 1) = 1
      XMF1(ichan, 1) = 0.
      if (flagsubfis .and. MT == 18) XMF1(ichan, 1) = 10
      if (MF == 40) XMF1(ichan, 1) = 10
      XLFS1(ichan, 1) = 0.
      if (flageaf) then
        LB8 = 1
      else
        if (MF == 33) then
          NI33(ichan, 1) = 2
          LB8 = 8
        endif
      endif
      MAT1(ichan, 1) = 0
      IZAP1(ichan, 1) = 0
      LS(ichan, 1) = 1
      LB(ichan, 1) = 5
      NT8(ichan) = 2 * NE8(ichan)
!
! Inter-MT correlations
!
      Nsec = 1
      if (MF == 33 .and. flagMTint(ichan)) then
        i1 = MTintindex(ichan)
        do MT2 = MT + 1, nummt
          do ichan2 = 1, Nchancov
            if (MT2 /= MTindex(ichan2)) cycle
            if ( .not. flagMTint(ichan2)) cycle
            if ( .not. mtexist(MF - 30, MT2)) cycle
            if (MT2 >= 51 .and. MT2 <= 91 .and. covdiscrete == 0) cycle
            if (MT2 >= 51 .and. MT2 <= 90 .and. MT2 > 50 + covdiscrete) cycle
            if (flagsubfis .and. MT2 == 18) cycle
            i2 = MTintindex(ichan2)
            Nsec = Nsec + 1
            NE33(ichan, Nsec) = NE33(ichan, 2)
            do iE6 = 1, NE33(ichan, Nsec)
              b33(i1, Nsec, iE6) = ELB6(iE6) * 1.e6
            enddo
            ib = NE33(ichan, Nsec)
            do iE6b = 1, NE2cov(ichan2)
              ib = ib + 1
              b33(i1, Nsec, ib) = ELB6b(ichan2, iE6b) * 1.e6
            enddo
            do iE6 = 2, NE33(ichan, Nsec)
              do iE6b = 2, NE2cov(ichan2)
                ib = ib + 1
                b33(i1, Nsec, ib) = Rfinal(i1, iE6, i2, iE6b)
              enddo
            enddo
            MT33(ichan, Nsec) = MTindex(ichan2)
            NC33(ichan, Nsec) = 0
            NT33(ichan, Nsec) = 1 + NE33(ichan, Nsec) * NE2cov(ichan2)
            NI33(ichan, Nsec) = 1
            MAT1(ichan, Nsec) = 0
            IZAP1(ichan, Nsec) = 1000 * Ztarget + Atarget
            XMF1(ichan, Nsec) = 0
            LS(ichan, Nsec) = 0
            LB(ichan, Nsec) = 6
          enddo
        enddo
      endif
      NL33(MF, MT) = Nsec
!
! 2. ENDF-6 parameters
!
      MTL(MT) = 0
    enddo
  400   if (flag33(MF, MT)) then
      mtexist(MF, MT) = .true.
      mfexist(MF) = .true.
    endif
  enddo
  return
end subroutine make33
! Copyright A.J. Koning 2021
