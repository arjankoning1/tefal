subroutine write1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF1
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
!   numen6         ! number of incident energies
!   numint         ! number of interpolation sections
!   nummf          ! number of MF numbers
!   nummt          ! number of MT numbers
! Variables for info from TALYS
!   k0             ! index of incident particle
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
! Variables for initialization of ENDF format
!   AWR            ! standard mass parameter
!   blank1         ! blank string
!   blank2         ! blank string
!   CONT           ! ENDF - 6 format
!   FEND           ! ENDF - 6 format
!   LIS            ! state number of target nucleus
!   LISO           ! isomeric state number
!   MAT            ! MAT number
!   mtexist        ! flag for existence of MT - number
!   SEND           ! ENDF - 6 format
!   TEXT           ! ENDF - 6 format
!   VALUE          ! ENDF - 6 format
!   ZA             ! standard charge parameter
! Variables for MF0
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF1
!   ALAB           ! Lab
!   AUTH           ! author(s)
!   AWI            ! mass of projectile in neutron units
!   Cnubar         ! coefficients of the polynomial expansion
!   DDATE          ! date of distribution
!   DEBDEL         ! uncertainty of total energy of delayed beta's
!   DEFR           ! uncertainty of kinetic energy of fragments
!   DEGD           ! uncertainty of total energy of delayed gamma rays
!   DEGP           ! uncertainty of total energy of prompt gamma rays
!   DEND           ! uncertainty of kinetic energy of delayed fission neutrons
!   DENP           ! uncertainty of kinetic energy of prompt fission neutrons
!   DENU           ! uncertainty of energy of neutrinos
!   DERN           ! uncertainty of total energy minus energy of neutrinos
!   DET            ! uncertainty of sum of all partial energies
!   E1             ! incident energy for MF1 (in ENDF - 6 format)
!   EBDEL          ! total energy of delayed beta's
!   EDATE          ! date of evaluationon
!   EFR            ! kinetic energy of fragments
!   EGD            ! total energy of delayed gamma rays
!   EGP            ! total energy of prompt gamma rays
!   ELIS           ! excitation energy of target nucleus
!   EMAX           ! upper limit of energy range for evaluation
!   END1           ! kinetic energy of delayed fission neutrons
!   ENDATE         ! date string
!   ENP            ! kinetic energy of prompt fission neutrons
!   ENU            ! energy of neutrinos
!   ERN            ! total energy minus energy of neutrinos
!   ET             ! sum of all partial energies
!   HSUB1          ! string
!   HSUB2          ! string
!   HSUB3          ! string
!   lambda         ! decay constant for precursor
!   LDRV           ! special evaluation flag
!   LFI            ! flag for fission
!   LNU            ! polynomial(LNU = 1) or tabulated (LNU = 2) representation
!   LREL           ! release number
!   LRP            ! flag that indicates presence of file 2 (resonances)
!   MODN           ! modification indicator
!   NC             ! number of lines for MF / MT number
!   NCnubar        ! number of terms in polynomial expansion
!   NFOR           ! library format
!   NLIB           ! library identifier
!   NMOD           ! modification number of evaluation
!   NNF            ! number of precursor families
!   NSUB           ! sub - library number
!   ntxt           ! number of text lines
!   nubar          ! average number of neutrons per fission
!   nuspont        ! nu for spontaneous fission
!   NVER           ! library version number
!   NWD            ! number of text lines
!   NXC            ! number of directory lines
!   RDATE          ! date of revision
!   REF            ! reference
!   STA            ! target stability flag
!   TEMP           ! target temperature
!   txt            ! text line
!   ZSYMAM         ! string
!
! *** Declaration of local data
!
  implicit none
  character(len=11) :: endf              ! function for a number in ENDF-6 format
  character(len=11) :: nus               ! nu for spontaneous fission
  integer           :: i                 ! counter
  integer           :: ii                ! counter
  integer           :: IMF               ! counter for MF-number
  integer           :: IMT               ! counter for MT-number
  integer           :: j                 ! counter
  integer           :: MF                ! MF-number
  integer           :: MT                ! MT-number
  integer           :: N                 ! neutron number of residual nucleus
  integer           :: NS                ! line number
  integer           :: nutype            ! help variable
  integer           :: Nval(2*numint)    ! value
  real(sgl)         :: x(2*numen6)       ! help variable
!
! ***************************** Write MF1 ******************************
!
! hrwrite  : subroutine to write header with real values
!
  MF = 1
  MT = 451
  NS = 0
  open (unit = 2, file = 'MF1', status = 'replace')
  call hrwrite(ZA, AWR, LRP, LFI, NLIB, NMOD, MAT, MF, MT, NS)
  call hrwrite(ELIS, STA, LIS, LISO, 0, NFOR, MAT, MF, MT, NS)
  call hrwrite(AWI, EMAX, LREL, 0, NSUB, NVER, MAT, MF, MT, NS)
  call hrwrite(TEMP, 0., LDRV, 0, NWD, NXC, MAT, MF, MT, NS)
  write(2, '(2a11, a10, 1x, a33, i4, i2, i3, i5)') ZSYMAM, ALAB, EDATE, AUTH, MAT, MF, MT, NS
  write(2, '(1x, a21, a10, 1x, a10, 12x, a6, 5x, i4, i2, i3, i5)') REF, DDATE, RDATE, ENDATE, MAT, MF, MT, NS
  write(2, fmt = TEXT) HSUB1, MAT, MF, MT, NS
  write(2, fmt = TEXT) HSUB2, MAT, MF, MT, NS
  write(2, fmt = TEXT) HSUB3, MAT, MF, MT, NS
  do i = 1, ntxt
    write(2, fmt = TEXT) txt(i), MAT, MF, MT, NS
  enddo
  do IMF = 1, nummf
    do IMT = 1, nummt
      if (mtexist(IMF, IMT)) then
        write(2, fmt = CONT) blank1, blank1, IMF, IMT, NC(IMF, IMT), MODN(IMF, IMT), MAT, MF, MT, NS
      endif
    enddo
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
!
! ********************** Fission neutrons ******************************
!
! rwrite     : subroutine to write line with real values
!
  if (k0 <= 1 .and. flagfission) then
    do MT = 452, 456
      if ( .not. mtexist(MF, MT)) cycle
      if (MT == 452) nutype = 1
      if (MT == 455) nutype = 2
      if (MT == 456) nutype = 3
      call hrwrite(ZA, AWR, 0, LNU(MT), 0, 0, MAT, MF, MT, NS)
      if (LNU(MT) == 1) then
        if (MT == 452) then
          N = NCnubar
          call hrwrite(0., 0., 0, 0, N, 0, MAT, MF, MT, NS)
          do i = 1, N
            x(i) = Cnubar(i)
          enddo
          call xwrite(N, x, MAT, MF, MT, NS)
        endif
        if (MT == 455) then
          N = NNF
          call hrwrite(0., 0., 0, 0, N, 0, MAT, MF, MT, NS)
          do i = 1, N
            x(i) = lambda(i)
          enddo
          call xwrite(N, x, MAT, MF, MT, NS)
        endif
        if (MT == 455 .or. MT == 456) then
          call hrwrite(0., 0., 0, 0, 1, 0, MAT, MF, MT, NS)
          nus = endf(nuspont)
          write(2, fmt = VALUE) nus, (blank1, j = 1, 5), MAT, MF, MT, NS
        endif
      else
        if (MT == 455) then
          N = NNF
          call hrwrite(0., 0., 0, 0, N, 0, MAT, MF, MT, NS)
          do i = 1, N
            x(i) = lambda(i)
          enddo
          call xwrite(N, x, MAT, MF, MT, NS)
        endif
        call hrwrite(0., 0., 0, 0, NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
! kwrite: subroutine to write integer value block
!
        do i = 1, NR(MF, MT)
          ii = 2 * i - 1
          Nval(ii) = NBT(MF, MT, i)
          Nval(ii + 1) = INTER(MF, MT, i)
        enddo
        N = 2 * NR(MF, MT)
        call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write tabulated values
!
! xwrite: subroutine to write real value block
!
        do i = 1, NP(MF, MT)
          ii = 2 * i - 1
          x(ii) = E1(nutype, i)
          x(ii + 1) = nubar(nutype, i)
        enddo
        N = 2 * NP(MF, MT)
        call xwrite(N, x, MAT, MF, MT, NS)
      endif
      write(2, fmt = SEND) blank2, MAT, MF, 0, NS
    enddo
!
! Components of energy release
!
    MT = 458
    if (mtexist(MF, MT)) then
      call hrwrite(ZA, AWR, 0, 0, 0, 0, MAT, MF, MT, NS)
      call hrwrite(0., 0., 0, 0, 18, 9, MAT, MF, MT, NS)
      call rwrite(EFR, DEFR, ENP, DENP, END1, DEND, MAT, MF, MT, NS)
      call rwrite(EGP, DEGP, EGD, DEGD, EBDEL, DEBDEL, MAT, MF, MT, NS)
      call rwrite(ENU, DENU, ERN, DERN, ET, DET, MAT, MF, MT, NS)
      write(2, fmt = SEND) blank2, MAT, MF, 0, NS
    endif
  endif
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write1
! Copyright A.J. Koning 2021
