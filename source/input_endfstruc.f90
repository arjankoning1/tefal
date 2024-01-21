subroutine input_endfstruc
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for ENDF structure
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
! All global variables
!   nummf           ! number of MF numbers
!   nummt           ! number of MT numbers
! Variables for TEFAL input
!   flagaddlow      ! flag to add low - energy Q>0 reactions to total cross section
!   flagcapt6       ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagdisc6       ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
!   flagexclude     ! flag to specifically exclude an MT number (per MT)
!   flagfis10       ! flag to put (subactinide) fission cross sections in MF10
!   flaggam13       ! flag to use MF13 for gamma prod. instead of MF12 (if not in MF6)
!   flaggamdisc     ! flag to store gamma prod. per level (MT51..) instead of MT4
!   flaggamspec     ! flag to store gamma prod. only as a spectrum and not per level
!   flaginclude     ! flag to specifically include an MT number (per MT)
!   flagmtall       ! flag to include all defined MT numbers from ENDF manual
!   flagmtextra     ! flag to include extra MT numbers up to MT200
!   flagmulti       ! flag to include multi - chance fission
!   flagngn         ! flag to include (n, gamma n) data
!   flagpara        ! flag to include partial cross sections for alphs
!   flagpard        ! flag to include partial cross sections for deuterons
!   flagparh        ! flag to include partial cross sections for helions
!   flagparn        ! flag to include partial cross sections for neutrons
!   flagparp        ! flag to include partial cross sections for protons
!   flagpart        ! flag to include partial cross sections for tritons
!   flagpart6       ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
!   flagrecoil      ! flag to include recoil information
!   flagrenorm      ! flag for renormalization of spectra
!   flagres         ! flag to include resonance parameters
!   flagrp10        ! flag to put residual production cross sections in MF10
!   flagrp6         ! flag to put residual production cross sections in MF6
!   flagsubfis      ! flag to include subactinide fission
!   flagtabddx      ! flag to give explicit DDX in MF6
!   nomf            ! flag to exclude an entire MF
! Variables for TALYS info
!   flagtalysrec    ! flag for recoil information from TALYS
!   k0              ! index of incident particle
! Variables for reading TEFAL input lines
!   inline          ! input line
!   nlines          ! number of input lines
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch          ! character
  character(len=132) :: key         ! keyword
  character(len=132) :: val         ! value or string
  character(len=132) :: word(40)    ! words on input line
  character(len=132) :: line        ! input line
  integer            :: i           ! counter
  integer            :: imf         ! MF counter
  integer            :: imt         ! MT counter
  integer            :: istat       ! logical for file access
!
! ********************************* Defaults ***************************
!
  flagpart6 = .true.
  flagparp = .false.
  flagpard = .false.
  flagpart = .false.
  flagparh = .false.
  flagpara = .false.
  flagcapt6 = .true.
  if (k0 <= 1) then
    flagparn = .true.
  else
    flagparn = .false.
  endif
  if (k0 == 1) then
    flagdisc6 = .false.
  else
    flagdisc6 = .true.
  endif
  flagrp6 = .true.
  flagres = .true.
  flagrp10 = .false.
  flagfis10 = .false.
  flagngn = .false.
  flaggam13 = .false.
  flaggamdisc = .true.
  flaggamspec = .false.
  flagmtall = .true.
  flagmtextra = .false.
  flagmulti = .true.
  flagrecoil = flagtalysrec
  flagsubfis = .false.
  flagtabddx = .false.
  flagrenorm = .true.
  flagaddlow = .true.
  nomf = .false.
  flaginclude = .false.
  flagexclude = .false.
!
! ************************ Read input variables ************************
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified and the corresponding values are read.
! Erroneous input is immediately checked.
! The keywords and number of values on each line are retrieved from the input.
!
  istat = 0
  do i = 1, nlines
    line = inline(i)
    call getkeywords(line, word)
    key = word(1)
    val = word(2)
    ch = word(2)(1:1)
!
! Test for keywords
!
    if (key == 'disc6') then
      if (ch == 'n') flagdisc6 = .false.
      if (ch == 'y') flagdisc6 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'part6') then
      if (ch == 'n') flagpart6 = .false.
      if (ch == 'y') flagpart6 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rp6') then
      if (ch == 'n') flagrp6 = .false.
      if (ch == 'y') flagrp6 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'resonance') then
      if (ch == 'n') flagres = .false.
      if (ch == 'y') flagres = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rp10') then
      if (ch == 'n') flagrp10 = .false.
      if (ch == 'y') flagrp10 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fis10') then
      if (ch == 'n') flagfis10 = .false.
      if (ch == 'y') flagfis10 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ngn') then
      if (ch == 'n') flagngn = .false.
      if (ch == 'y') flagngn = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'capt6') then
      if (ch == 'n') flagcapt6 = .false.
      if (ch == 'y') flagcapt6 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gam13') then
      if (ch == 'n') flaggam13 = .false.
      if (ch == 'y') flaggam13 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gamdiscrete') then
      if (ch == 'n') flaggamdisc = .false.
      if (ch == 'y') flaggamdisc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partialn') then
      if (ch == 'n') flagparn = .false.
      if (ch == 'y') flagparn = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partialp') then
      if (ch == 'n') flagparp = .false.
      if (ch == 'y') flagparp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partiald') then
      if (ch == 'n') flagpard = .false.
      if (ch == 'y') flagpard = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partialt') then
      if (ch == 'n') flagpart = .false.
      if (ch == 'y') flagpart = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partialh') then
      if (ch == 'n') flagparh = .false.
      if (ch == 'y') flagparh = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partiala') then
      if (ch == 'n') flagpara = .false.
      if (ch == 'y') flagpara = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'mtall') then
      if (ch == 'n') flagmtall = .false.
      if (ch == 'y') flagmtall = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'mtextra') then
      if (ch == 'n') flagmtextra = .false.
      if (ch == 'y') flagmtextra = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'multichance') then
      if (ch == 'n') flagmulti = .false.
      if (ch == 'y') flagmulti = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'recoil') then
      if (ch == 'n') flagrecoil = .false.
      if (ch == 'y') flagrecoil = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'subfission') then
      if (ch == 'n') flagsubfis = .false.
      if (ch == 'y') flagsubfis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'tabddx') then
      if (ch == 'n') flagtabddx = .false.
      if (ch == 'y') flagtabddx = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'renorm') then
      if (ch == 'n') flagrenorm = .false.
      if (ch == 'y') flagrenorm = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'addlow') then
      if (ch == 'n') flagaddlow = .false.
      if (ch == 'y') flagaddlow = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'include') then
      read(val, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      flaginclude(imt) = .true.
      cycle
    endif
    if (key == 'exclude') then
      read(val, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      flagexclude(imt) = .true.
      cycle
    endif
    if (key == 'nomf') then
      read(val, * , iostat = istat) imf
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imf, 0, nummf)
      nomf(imf) = .true.
      cycle
    endif
  enddo
  return
end subroutine input_endfstruc
! Copyright A.J. Koning 2021
