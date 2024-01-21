subroutine read4(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF4 from existing ENDF-6 data library
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
!   sgl          ! single precision kind
! All global variables
!   numenang     ! number of incident energies for angular distributions
!   numenin      ! number of incident energies
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   adoptfile    ! name of library for MF information (per MT)
!   Eahigh       ! upper energy of MT values to be adopted
!   Ealow        ! lower energy of MT values to be adopted
! Variables for input of ENDF library type
!   flaghigh     ! flag for high energies ( > 20 MeV)
! Variables for initialization of ENDF format
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
! Variables for MF1
!   EMAX         ! upper limit of energy range for evaluation
! Variables for ENDF format
!   INTER        ! interpolation scheme
!   NBT          ! separation value for interpolation scheme
!   NR           ! number of interpolation ranges
! Variables for MF4
!   E4hr         ! incident energy for MF4 (in ENDF - 6 format)
!   E4r          ! incident energy for MF4 (in ENDF - 6 format)
!   f4r          ! angular distribution
!   INTER4       ! interpolation scheme
!   INTERh       ! interpolation scheme
!   LCT          ! LAB / CM flag
!   legr         ! Legendre coefficients (in ENDF - 6 format)
!   LI4          ! isotropy flag
!   LTT          ! representation
!   LVT          ! specification of transformation matrix
!   NBT4         ! separation value for interpolation scheme
!   NBTh         ! separation value for interpolation scheme
!   NE4r         ! number of incident energies (MF4 only)
!   NEhr         ! number of incident energies (MF4 only)
!   NL4r         ! number of Legendre coefficients
!   NP4r         ! number of angles (MF4 only)
!   NR4r         ! number of interpolation ranges (MF4 only)
!   NRh          ! number of interpolation ranges
!   x4r          ! cosine of the angle
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: iE        ! energy counter
  integer           :: ii        ! counter
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: k         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: nen       ! energy counter
  integer           :: nlin      ! number of lines
  integer           :: NLx       ! number of Legendre coefficients
  integer           :: NPx       ! number of angles
  integer           :: NRx       ! number of interpolation ranges
  real(sgl)         :: Eah       ! incident energy
  real(sgl)         :: E4x       ! incident energy for MF4
!
! ********* Read data from particular MT-numbers from MF4 files ********
!
  MF = 4
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
!
! Read until MT is found
!
  MTstr = '     '
  write(MTstr(1:2), '(i2)') MF
  write(MTstr(3:5), '(i3)') MT
  do
    read(3, '(a80)', iostat = istat) string
    if (istat /= 0) then
      close (unit = 3)
      return
    endif
    if (string(71:75) == MTstr) exit
  enddo
  read(string(23:44), '(2i11)') LVT, LTT
  read(3, '(22x, 2i11)', iostat = istat) LI4, LCT
  if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  if (LTT == 2) goto 110
  read(3, '(44x, 2i11)', iostat = istat) NR(MF, MT), nen
  if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  if (LTT == 0) goto 150
!
! ********************** Legendre coefficients *************************
!
  nlin = 1 + (NR(MF, MT) - 1) / 3
  do i = 1, nlin
    ii = 3 * (i - 1)
    read(3, '(6i11)', iostat = istat) (NBT(MF, MT, j), INTER(MF, MT, j), j = ii+1, ii+3)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  enddo
  iE = 0
  do k = 1, nen
    read(3, '(11x, e11.6, 22x, i11)', iostat = istat) E4x, NLx
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    iE = iE + 1
    if (iE > numen4) then
      write(*, '(" TEFAL-error: there are more than ", i4, " incident energies in file ", a)') numen4, trim(adoptfile(MF, MT))
      write(*, '(" numen4 in tefal.cmb should be increased")')
      stop
    endif
    E4r(iE) = E4x
    NL4r(iE) = NLx
    nlin = 1 + (NLx - 1) / 6
    do i = 1, nlin
      ii = 6 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (legr(iE, j), j = ii+1, ii+6)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    if (E4x < Ealow(MF, MT) .or. E4x > Eahigh(MF, MT) .or. E4x > EMAX) iE = iE - 1
  enddo
  NE4r = iE
  if (IE > 0) then
    Eah = E4r(iE)
  else
    Eah = 0.
  endif
!
! *********** Tabulated angular distributions (high energies) **********
!
  110 nen = 0
  if (LTT == 2 .or. (flaghigh .and. LTT == 3)) then
    read(3, '(44x, 2i11)', iostat = istat) NRh, nen
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    if (nen > numenang) then
      write(*, '(" TEFAL-warning: there are more than", i3, " incident energies in file ", a)') numenang, trim(adoptfile(MF, MT))
      write(*, '(" numenang in A0_tefal.mod should be increased, now ", i4, " adopted")') numenang
    endif
    nlin = 1 + (NRh - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBTh(j), INTERh(j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
  endif
  iE = 0
  do k = 1, min(nen, numenang)
    read(3, '(11x, e11.6, 22x, 2i11)', iostat = istat) E4x, NRx, NPx
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    iE = iE + 1
    E4hr(iE) = E4x
    NR4r(iE) = NRx
    NP4r(iE) = NPx
    do i = 1, NR4r(iE)
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBT4(iE, j), INTER4(iE, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    nlin = 1 + (NP4r(iE) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (x4r(iE, j), f4r(iE, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    if (E4x < Ealow(MF, MT) .or. E4x > Eahigh(MF, MT) .or. E4x > EMAX) iE = iE - 1
  enddo
  NEhr = iE
  if (iE > 0) Eah = E4hr(iE)
  150 mtexist(MF, MT) = .true.
  mfexist(MF) = .true.
  Eahigh(MF, MT) = Eah
  if ( .not. flaghigh .and. LTT /= 2) LTT = 1
  close (unit = 3)
  return
end subroutine read4
! Copyright A.J. Koning 2021
