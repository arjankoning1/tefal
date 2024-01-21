subroutine input_endfspec
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for specific ENDF data
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
!   nummf          ! number of MF numbers
!   nummt          ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt         ! logical for existence of MF information (per MT)
!   adoptfile     ! name of library for MF information (per MT)
!   background    ! file with background cross sections
!   lssfinp       ! 0: URR cross section from MF2, 1: URR cross section
!   urrcomp       ! mode for competition in the URR, 0: none, 1:MT4, 2:all
!   urrmode       ! 0: no URR, 1: URR from TALYS, 2: URR from data library
!   Eahigh        ! upper energy of MT values to be adopted
!   Ealow         ! lower energy of MT values to be adopted
!   urrenergy     ! upper energy of the URR in MeV
! Variables for TALYS info
!   Atarget      ! mass number of nucleus
!   flagtalysurr ! flag for URR information from TALYS
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
  integer            :: imf         ! MF counter
  integer            :: imt         ! MT counter
  integer            :: istat       ! logical for file access
!
! ********************************* Defaults ***************************
!
  adopt = .false.
  adoptfile = ' '
  Ealow = 1.e-5
  Eahigh = 2.e+8
  background = ' '
  urrcomp = 1
  urrmode = 0
  urrenergy = -1.
  lssfinp = -1
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
    if (key == 'urrcomp') then
      read(val, * , iostat = istat) urrcomp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'urrmode') then
      read(val, * , iostat = istat) urrmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'lssf') then
      read(val, * , iostat = istat) lssfinp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'urrenergy') then
      read(val, * , iostat = istat) urrenergy
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'adopt') then
      read(val, * , iostat = istat) imf
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imf, 1, nummf)
      if (word(3)(1:1) >= '0' .and. word(3)(1:1) <= '9') then
        read(word(3), * , iostat = istat) imt
        if (istat /= 0) call read_error(line, istat)
        call range_integer_error(key, imt, 1, nummt)
        adopt(imf, imt) = .true.
        adoptfile(imf, imt) = word(4)
        if (word(5) /= '') then
          read(word(5), * , iostat = istat) Ealow(imf, imt)
          if (istat /= 0) call read_error(word(5), istat)
        endif
        if (word(6) /= '') then
          read(word(6), * , iostat = istat) Eahigh(imf, imt)
          if (istat /= 0) call read_error(word(6), istat)
        endif
      else
        do imt = 1, nummt
          adopt(imf, imt) = .true.
          adoptfile(imf, imt) = word(3)
          if (word(4) /= '') then
            read(word(4), * , iostat = istat) Ealow(imf, imt)
            if (istat /= 0) call read_error(word(4), istat)
          endif
          if (word(5) /= '') then
            read(word(5), * , iostat = istat) Eahigh(imf, imt)
            if (istat /= 0) call read_error(word(5), istat)
          endif
        enddo
      endif
      cycle
    endif
    if (key == 'background') then
      background = val
      cycle
    endif
  enddo
  return
end subroutine input_endfspec
! Copyright A.J. Koning 2021
