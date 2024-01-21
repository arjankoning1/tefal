subroutine finalwrite
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write final ENDF-6 file
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
! All global variables
!   nummf       ! number of MF numbers
! Variables for input of ENDF library type
!   endffile    ! name of ENDF file
!   flageaf     ! flag for EAF - formatted activation library
! Variables for initialization of ENDF format
!   blank2      ! blank string
!   MEND        ! ENDF - 6 format
!   mfexist     ! flag for existence of MF - number
!   rec0        ! first line of ENDF - 6 file
!   TEND        ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=4)  :: mffile    ! MF filename
  character(len=80) :: string    ! line with parameter value
  integer           :: imf       ! MF counter
  integer           :: imt       ! MT counter
  integer           :: istat     ! logical for file access
  integer           :: MF        ! MF-number
  integer           :: NS        ! line number
!
! ************************* Write total file ***************************
!
  open (unit = 10, file = trim(endffile), status = 'replace')
  if ( .not. flageaf) write(10, '(a80)') rec0
  do MF = 1, nummf
    if ( .not. mfexist(MF)) cycle
    if (flageaf .and. MF >= 30) cycle
    mffile = 'MF  '
    if (MF < 10) then
      write(mffile(3:3), '(i1)') MF
    else
      write(mffile(3:4), '(i2)') MF
    endif
    inquire (file = mffile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = mffile, status = 'old', iostat = istat)
      NS = 0
      do
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit
        read(string(71:72), '(i2)') imf
        read(string(73:75), '(i3)') imt
!
! The EAF-library does not contain line numbers.
!
        if (flageaf) then
          write(string(76:80), '(5x)')
        else
!
! Begin of new MF.
!
          if (imf == 0) then
            write(string(76:80), '(i5)') 0
            write(10, '(a80)') string
            exit
          endif
!
! End of MT-section.
!
          if (imt == 0) then
            write(string(76:80), '(i5)') 99999
            NS = 0
          else
!
! Increasing line numbers per MT-section.
!
            NS = NS + 1
            write(string(76:80), '(i5)') NS
          endif
        endif
        write(10, '(a80)') string
      enddo
      close (unit = 1)
    endif
  enddo
  write(10, fmt = MEND) blank2, 0, 0, 0, 0
  if ( .not. flageaf) write(10, fmt = TEND) blank2, -1, 0, 0, 0
  close (unit = 10)
  return
end subroutine finalwrite
! Copyright A.J. Koning 2021
