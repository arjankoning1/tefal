subroutine input_endflimits
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for ENDF limits, switches and tolerances
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
! Variables for ENDF limits, switches and tolerances
!   maxrp       ! maximum number of residual products
!   NMTmax      ! maximum number of MT numbers
!   cuteps      ! energy shift at MT5 cutoff energy (in eV)
!   disclim     ! limit for specific MT numbers for discrete levels
!   Eswitch     ! energy where ENDF-6 representation is switched (in MeV)
!   Eswitch4    ! energy where MF4 representation is switched (in MeV)
! Variables for info from TALYS
!   flagtalysdet ! flag for detailed ENDF - 6 information from TALYS
! Variables for reading input lines
!   inline          ! input line
!   nlines          ! number of input lines
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
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
  integer            :: istat       ! logical for file access
!
! ********************************* Defaults ***************************
!
  if (flagtalysdet) then
    Eswitch = 30.
    Eswitch4 = 50.
  else
    Eswitch = 0.
    Eswitch4 = 0.
  endif
  disclim = 10.
  cuteps = 0.
  maxrp = 100
  NMTmax = 400
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
    if (key == 'disclim') then
      read(val, * , iostat = istat) disclim
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'eswitch') then
      read(val, * , iostat = istat) Eswitch
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'eswitch4') then
      read(val, * , iostat = istat) Eswitch4
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'cuteps') then
      read(val, * , iostat = istat) cuteps
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxrp') then
      read(val, * , iostat = istat) maxrp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'nmtmax') then
      read(val, * , iostat = istat) NMTmax
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_endflimits
! Copyright A.J. Koning 2021
