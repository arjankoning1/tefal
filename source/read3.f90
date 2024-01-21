subroutine read3
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF3 from existing ENDF-6 data library
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
!   numen6       ! number of incident energies
!   nummt        ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   adoptfile    ! name of library for MF information (per MT)
! Variables for MF3
!   E3adopt      ! incident energy adopted from other library
!   NE3adopt     ! number of adopted incident energies
!   xs3adopt     ! cross section adopted from other library
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none

  integer, parameter :: numdat=1000000   ! number of data points
  character(len=5)   :: MTstr            ! string for MT number
  character(len=80)  :: string           ! line with parameter value
  character(len=132) :: afile            ! name of library for MF information (per MT)
  integer            :: i                ! counter
  integer            :: ii               ! counter
  integer            :: istat            ! logical for file access
  integer            :: j                ! counter
  integer            :: k                ! counter
  integer            :: MF               ! MF-number
  integer            :: MT               ! MT-number
  integer            :: Nmo              ! number of value on line
  integer            :: nlin             ! number of lines
  integer            :: NP3              ! number of data points
  integer            :: NR3              ! number of interpolation ranges
  real(sgl)          :: x1(numdat)       ! coordinates of intersection points inside the bin
  real(sgl)          :: y1(numdat)       ! coordinates of the 1st summit of the triangle
!
! ********* Read data from particular MT-numbers from MF3 files ********
!
  MF = 3
  Loop1: do MT = 1, nummtres
    if ( .not. adopt(MF, MT)) cycle
    afile = adoptfile(MF, MT)
    open (unit = 3, file = afile, status = 'old', iostat= istat)
    if (istat /= 0) call read_error(afile, istat)
!
! Read until MT is found
!
    MTstr = '     '
    write(MTstr(1:2), '(i2)') MF
    write(MTstr(3:5), '(i3)') MT
    do
      read(3, '(a80)', iostat = istat) string
      if (istat ==  -1) cycle Loop1
      if (string(71:75) == MTstr) exit
    enddo
    read(3, '(a80)', iostat = istat) string
    if (istat /= 0) call read_error(afile, istat)
    read(string(45:55), '(i11)') NR3
    read(string(56:66), '(i11)') NP3
    nlin = 1 + (NR3 - 1) / 3
    do i = 1, nlin
      read(3, '(a80)', iostat = istat) string
      if (istat /= 0) call read_error(afile, istat)
    enddo
    NP3 = min(NP3, numdat - 100)
    nlin = 1 + (NP3 - 1) / 3
    do ii = 1, nlin
      j = 3 * (ii - 1)
      read(3, '(6e11.6)', iostat = istat) (x1(j+k), y1(j+k), k = 1, 3)
      if (istat /= 0) call read_error(afile, istat)
    enddo
    close (unit = 3)
    if (NP3 <= numen6) then
      do i = 1, NP3
        E3adopt(MT, i) = x1(i)
        xs3adopt(MT, i) = y1(i)
      enddo
      NE3adopt(MT) = NP3
    else
      Nmo = NP3 / numen6 + 1
      E3adopt(MT, 1) = x1(1)
      xs3adopt(MT, 1) = y1(1)
      i = 1
      do k = 2, NP3 - 1
        if (mod(k, Nmo) == 0) then
          i = i + 1
          E3adopt(MT, i) = x1(k)
          xs3adopt(MT, i) = y1(k)
        endif
      enddo
      i = i + 1
      E3adopt(MT, i) = x1(NP3)
      xs3adopt(MT, i) = y1(NP3)
      NE3adopt(MT) = i
    endif
  enddo Loop1
  return
end subroutine read3
! Copyright A.J. Koning 2021
