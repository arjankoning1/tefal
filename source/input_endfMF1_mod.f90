subroutine input_endfMF1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for ENDF MF1 keywords
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
! Variables for input of ENDF MF1
!   author          ! author
!   diffweight      ! differential weight to be put in MF1
!   endftext        ! file with MF1 information
!   identifier      ! library identifier
!   lab             ! laboratory
! Variables for TALYS info
!   k0              ! index of incident particle
! Variables for reading TEFAL input lines
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
  diffweight = -1.
  if (k0 == 1) then
    author = 'A.J. Koning and D. Rochman'
  else
    author = 'A.J. Koning               '
  endif
  lab = 'IAEA       '
  identifier = 'TENDL-2023'
  endftext = ''
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
    if (key == 'author') then
      author = val(1:33)
      cycle
    endif
    if (key == 'lab') then
      lab = val(1:11)
      cycle
    endif
    if (key == 'identifier') then
      identifier = val(1:10)
      cycle
    endif
    if (key == 'endftext') then
      endftext = val
      cycle
    endif
    if (key == 'diffweight') then
      read(val, * , iostat = istat) diffweight
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_endfMF1
! Copyright A.J. Koning 2021
