subroutine input_endftype
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for ENDF library type
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
! Variables for input of ENDF library type
!   flageaf        ! flag for EAF-formatted activation library
!   flagendfdet    ! flag for detailed ENDF-6 information per channel
!   flagclean      ! flag to clean up double points
!   flaggpf        ! flag for general purpose library
!   flagbreakup    ! breakup flag
!   flaghigh       ! flag for high energies ( > 20 MeV)
!   endffile       ! name of ENDF file
! Variables for TALYS info
!   flagtalysdet ! flag for detailed ENDF - 6 information from TALYS
! Variables for reading input lines
!   inline          ! input line
!   nlines          ! number of input lines
! Error handling
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
  flaggpf = .true.
  flageaf = .false.
  flaghigh = .true.
  flagbreakup = .false.
  flagclean = .true.
  flagendfdet = flagtalysdet
  endffile = ' '
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
    if (key == 'gpf') then
      if (ch == 'n') flaggpf = .false.
      if (ch == 'y') flaggpf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'eaf') then
      if (ch == 'n') flageaf = .false.
      if (ch == 'y') flageaf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'high') then
      if (ch == 'n') flaghigh = .false.
      if (ch == 'y') flaghigh = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'breakup') then
      if (ch == 'n') flagbreakup = .false.
      if (ch == 'y') flagbreakup = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'clean') then
      if (ch == 'n') flagclean = .false.
      if (ch == 'y') flagclean = .true.
      if (ch /= 'y' .and. ch /= 'n')  call read_error(line, istat)
      cycle
    endif
    if (key == 'endfdetail') then
      if (ch == 'n') flagendfdet = .false.
      if (ch == 'y') flagendfdet = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'endffile') then
      endffile = val
      cycle
    endif
  enddo
  return
end subroutine input_endftype
! Copyright A.J. Koning 2021
