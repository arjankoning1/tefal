subroutine read1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF1 (fission neutrons) from existing ENDF-6 data
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
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   adoptfile    ! name of library for MF information (per MT)
! Variables for TALYS info
!   k0           ! index of incident particle
! Variables for initialization of ENDF format
!   mtexist      ! flag for existence of MT - number
! Variables for ENDF format
!   INTER        ! interpolation scheme
!   NBT          ! separation value for interpolation scheme
!   NP           ! number of incident energies
!   NR           ! number of interpolation ranges
! Variables for MF1
!   Cnubar       ! coefficients of the polynomial expansion
!   DEBDEL       ! uncertainty of total energy of delayed beta's
!   DEFR         ! uncertainty of kinetic energy of fragments
!   DEGD         ! uncertainty of total energy of delayed gamma rays
!   DEGP         ! uncertainty of total energy of prompt gamma rays
!   DEND         ! uncertainty of kinetic energy of delayed fission neutrons
!   DENP         ! uncertainty of kinetic energy of prompt fission neutrons
!   DENU         ! uncertainty of energy of neutrinos
!   DERN         ! uncertainty of total energy minus energy of neutrinos
!   DET          ! uncertainty of sum of all partial energies
!   E1           ! incident energy for MF1 (in ENDF - 6 format)
!   EBDEL        ! total energy of delayed beta's
!   EFR          ! kinetic energy of fragments
!   EGD          ! total energy of delayed gamma rays
!   EGP          ! total energy of prompt gamma rays
!   END1         ! kinetic energy of delayed fission neutrons
!   ENP          ! kinetic energy of prompt fission neutrons
!   ENU          ! energy of neutrinos
!   ERN          ! total energy minus energy of neutrinos
!   ET           ! sum of all partial energies
!   lambda       ! decay constant for precursor
!   LNU          ! polynomial(LNU = 1) or tabulated (LNU = 2) representation
!   NC           ! number of lines for MF / MT number
!   NCnubar      ! number of terms in polynomial expansion
!   NNF          ! number of precursor families
!   nubar        ! average number of neutrons per fission
!   nuspont      ! nu for spontaneous fission
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)   :: MTstr     ! string for MT number
  character(len=80)  :: string    ! line with parameter value
  character(len=132) :: afile    ! name of library for MF information (per MT)
  integer            :: i         ! counter
  integer            :: ii        ! counter
  integer            :: istat     ! error code
  integer            :: j         ! counter
  integer            :: MF        ! MF-number
  integer            :: MT        ! MT-number
  integer            :: nlin      ! number of lines
  integer            :: nutype    ! help variable
!
! ***************************** Find ENDF file *************************
!
  MF = 1
  Loop1: do MT = 452, 458
    if (MT == 453 .or. MT == 454 .or. MT == 457) cycle
    if ( .not. adopt(MF, MT)) cycle
    afile = adoptfile(MF, MT)
    open (unit = 3, file = adoptfile(MF, MT), status = 'old', iostat = istat)
    if (istat /= 0) call read_error(afile, istat)
!
! Read until MF1/MT452 is found
!
    MTstr = '     '
    write(MTstr(1:2), '(i2)') MF
    write(MTstr(3:5), '(i3)') MT
    do
      read(3, '(a80)', iostat = istat) string
      if (istat ==  -1) cycle Loop1
      if (string(71:75) == MTstr) exit
    enddo
    read(string(34:44), '(i11)') LNU(MT)
!
! LNU=1: polynomial representation
!
    if (MT == 452) nutype = 1
    if (MT == 455) nutype = 2
    if (MT == 456) nutype = 3
    NC(MF, MT) = 1
    if (LNU(MT) == 1) then
      if (MT == 452) then
        read(3, '(44x, i11)', iostat = istat) NCnubar
        if (istat /= 0) call read_error(afile, istat)
        nlin = 1 + (NCnubar - 1) / 6
        do i = 1, nlin
          ii = 6 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (Cnubar(j), j = ii+1, ii+6)
          if (istat /= 0) call read_error(afile, istat)
        enddo
        NC(MF, MT) = NC(MF, MT) + nlin + 1
      endif
      if (MT == 455) then
        read(3, '(44x, i11)', iostat = istat) NNF
        if (istat /= 0) call read_error(afile, istat)
        nlin = 1 + (NNF - 1) / 6
        do i = 1, nlin
          ii = 6 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (lambda(j), j = ii+1, ii+6)
          if (istat /= 0) call read_error(afile, istat)
        enddo
        NC(MF, MT) = NC(MF, MT) + nlin + 1
      endif
      if (MT == 455 .or. MT == 456) then
        read(3, '(/e11.6)', iostat = istat) nuspont
        if (istat /= 0) call read_error(afile, istat)
        NC(MF, MT) = NC(MF, MT) + 2
      endif
    endif
!
! LNU=2: tabulated representation
!
    if (LNU(MT) == 2) then
      if (MT == 455) then
        read(3, '(44x, i11)', iostat = istat) NNF
        if (istat /= 0) call read_error(afile, istat)
        nlin = 1 + (NNF - 1) / 6
        do i = 1, nlin
          ii = 6 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (lambda(j), j = ii+1, ii+6)
          if (istat /= 0) call read_error(afile, istat)
        enddo
        NC(MF, MT) = NC(MF, MT) + nlin + 1
      endif
      read(3, '(44x, 2i11)', iostat = istat) NR(MF, MT), NP(MF, MT)
      if (istat /= 0) call read_error(afile, istat)
      nlin = 1 + (NR(MF, MT) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6i11)', iostat = istat) (NBT(MF, MT, j), INTER(MF, MT, j), j = ii+1, ii+3)
        if (istat /= 0) call read_error(afile, istat)
      enddo
      NC(MF, MT) = NC(MF, MT) + nlin + 1
      nlin = 1 + (NP(MF, MT) - 1) / 3
      do i = 1, nlin
        ii = 3 * (i - 1)
        read(3, '(6e11.6)', iostat = istat) (E1(nutype, j), nubar(nutype, j), j = ii + 1, ii + 3)
        if (istat /= 0) call read_error(afile, istat)
      enddo
      NC(MF, MT) = NC(MF, MT) + nlin
    endif
!
! Components of energy release
!
    if (MT == 458) then
      if (k0 == 0) return
      read(3, '(/6e11.6)', iostat = istat) EFR, DEFR, ENP, DENP, END1, DEND
      if (istat /= 0) call read_error(afile, istat)
      read(3, '(6e11.6)', iostat = istat) EGP, DEGP, EGD, DEGD, EBDEL, DEBDEL
      if (istat /= 0) call read_error(afile, istat)
      read(3, '(6e11.6)', iostat = istat) ENU, DENU, ERN, DERN, ET, DET
      if (istat /= 0) call read_error(afile, istat)
      NC(MF, MT) = 5
    endif
    close (unit = 3)
    mtexist(MF, MT) = .true.
  enddo Loop1
  if (mtexist(MF, 452) .and. mtexist(MF, 455)) mtexist(MF, 456) = .true.
  if (mtexist(MF, 452) .and. mtexist(MF, 456)) mtexist(MF, 455) = .true.
  return
end subroutine read1
! Copyright A.J. Koning 2021
