subroutine readinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read user input
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
!   numlines      ! number of input lines
! Variables to read TEFAL input
!   inline        ! input line
!   nlines        ! number of input lines
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: str   ! input line
  integer            :: i     ! counter
  integer            :: istat ! logical for file access
  integer            :: k     ! counter
!
! ************************** User Input ********************************
!
! We read the complete input file first as a set of character strings.
! The actual keywords will be read from these later on.
!
  i = 1
  do
    read(*, '(a132)', iostat = istat) inline(i)
    if (istat ==  -1) exit
    if (istat /= 0) call read_error(inline(i), istat)
    i = i + 1
    call range_integer_error('inline', i, 1, numlines)
  enddo
  nlines = i - 1
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input is converted to lowercase characters, with the exception of
! filenames or other character strings.
!
Loop1:  do i = 1, nlines
    str = inline(i)
    do k = 1, 132
      if (inline(i)(k:k) >= 'A' .and. inline(i)(k:k) <= 'Z') inline(i)(k:k) = achar(iachar(inline(i)(k:k)) + 32)
    enddo
    do k = 0, 60
      if (inline(i)(k+1:k+3) == 'lab') then
        inline(i)(k + 4:132) = str(k + 4:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+5) == 'adopt') then
        inline(i)(k + 6:132) = str(k + 6:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+6) == 'author') then
        inline(i)(k + 7:132) = str(k + 7:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+8) == 'endffile') then
        inline(i)(k + 9:132) = str(k + 9:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+8) == 'endftext') then
        inline(i)(k + 9:132) = str(k + 9:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+10) == 'identifier') then
        inline(i)(k + 11:132) = str(k + 11:132)
        cycle Loop1
      endif
      if (inline(i)(k+1:k+10) == 'background') then
        inline(i)(k + 11:132) = str(k + 11:132)
        cycle Loop1
      endif
    enddo
  enddo Loop1
  return
end subroutine readinput
! Copyright A.J. Koning 2021
