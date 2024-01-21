subroutine input_endfcovar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for ENDF covariance keywords
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
!   nummt           ! number of MT numbers
! Variables for info from TALYS
!   k0              ! index of incident particle
! Variables for ENDF covariance input
!   covdiscrete     ! number of disc. inelastic levels with covariances
!   flagcovar       ! flag for covariances
!   flagcovrp       ! flag for covariance of residual production c.s.
!   flagcross       ! flag for covariance cross - channel correlation
!   flagintercor    ! flag for inter - MT covariance data
!   flagparcov      ! flag to include covariances for MT600 - 849
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
  integer            :: imt         ! MT counter
  integer            :: istat       ! logical for file access
!
! ********************************* Defaults ***************************
!
  flagparcov = .false.
  flagintercor = .true.
  flagcovrp = .false.
  flagcovar = .false.
  covdiscrete = 4
  do imt = 1, nummt
    flagcross(imt) = .false.
  enddo
  flagcross(1) = .true.
  flagcross(2) = .true.
  flagcross(3) = .true.
  flagcross(4) = .true.
  flagcross(16) = .true.
  flagcross(18) = .true.
  flagcross(51) = .true.
  flagcross(91) = .true.
  flagcross(102) = .true.
  flagcross(103) = .true.
  if (k0 == 3) flagcross(104) = .true.
  if (k0 == 4) flagcross(105) = .true.
  if (k0 == 5) flagcross(106) = .true.
  flagcross(107) = .true.
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
    if (key == 'partialcov') then
      if (ch == 'n') flagparcov = .false.
      if (ch == 'y') flagparcov = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'intercor') then
      if (ch == 'n') flagintercor = .false.
      if (ch == 'y') flagintercor = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'covrp') then
      if (ch == 'n') flagcovrp = .false.
      if (ch == 'y') flagcovrp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'covariance') then
      if (ch == 'n') flagcovar = .false.
      if (ch == 'y') flagcovar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'covdiscrete') then
      read(val, * , iostat = istat) covdiscrete
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'cross') then
      read(val, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      flagcross(imt) = .true.
      cycle
    endif
    if (key == 'nocross') then
      read(val, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      flagcross(imt) = .false.
      cycle
    endif
  enddo
  return
end subroutine input_endfcovar
! Copyright A.J. Koning 2021
