subroutine readbackground
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF3 background from existing ENDF-6 data library
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
!   sgl            ! single precision kind
! All global variables
!   numen6         ! number of incident energies
! Variables for input of specific ENDF data
!   background     ! file with background cross sections
! Variables for MF3
!   E3res          ! energy in resonance range
!   NE3res         ! number of energies in resonance range
!   xs3res         ! cross section in resonance range
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
! Variables for reaction initialization
!   EHres          ! upper energy in resonance range
!   includeres     ! flag to include resonance parameters
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist            ! logical to determine existence
  character(len=5)  :: MTstr             ! string for MT number
  character(len=80) :: string            ! line with parameter value
  integer           :: i                 ! counter
  integer           :: istat             ! error code
  integer           :: j                 ! counter
  integer           :: MF                ! MF-number
  integer           :: MT                ! MT-number
  integer           :: nen               ! energy counter
  integer           :: nlin              ! number of lines
  integer           :: NR3               ! number of interpolation ranges
  real(sgl)         :: deltaE            ! help variable
  real(sgl)         :: ea                ! help variable
  real(sgl)         :: eb                ! help variable
  real(sgl)         :: ee                ! energy
  real(sgl)         :: eloc(0:numen6)    ! help variable
  real(sgl)         :: x1                ! coordinates of intersection points inside the bin
  real(sgl)         :: x2                ! coordinates of the 2nd summit of the triangle
  real(sgl)         :: x3                ! coordinates of the 3rd summit of the triangle
  real(sgl)         :: y1                ! coordinates of the 1st summit of the triangle
  real(sgl)         :: y2                ! coordinates of the 2nd summit of the triangle
  real(sgl)         :: y3                ! coordinates of the 3rd summit of the triangle
!
! ****************** Read data in resonance range from MF2 *************
!
  inquire (file = background, exist = lexist)
  if ( .not. lexist) return
  open (unit = 3, file = background, status = 'old')
!
! Read until MF3 is found
!
  MF = 3
  do MT = 1, 107
    if ( .not. (MT == 1 .or. MT == 2 .or. (flagfission .and. MT == 18) .or. MT == 102)) cycle
    if ( .not. includeres) cycle
    MTstr = '     '
    write(MTstr(1:2), '(i2)') MF
    write(MTstr(3:5), '(i3)') MT
   20   read(3, '(a80)', end = 70, err = 70) string
    if (string(71:75) /= MTstr) goto 20
    read(3, '(a80)') string
    read(string(45:55), '(i11)') NR3
    nlin = 1 + (NR3 - 1) / 3
    do i = 1, nlin
      read(3, '(a80)') string
    enddo
    i = 0
    do
      read(3, '(6e11.6)', iostat = istat) x1, y1, x2, y2, x3, y3
      if (istat == -1) exit
      if (istat /= 0) call read_error(background, istat)
      if (x1 <= EHres) then
        i = i + 1
        E3res(MT, i) = x1
        xs3res(MT, i) = y1
      else
        exit
      endif
      if (x2 <= EHres) then
        i = i + 1
        E3res(MT, i) = x2
        xs3res(MT, i) = y2
      else
        exit
      endif
      if (x3 <= EHres) then
        i = i + 1
        E3res(MT, i) = x3
        xs3res(MT, i) = y3
      else
        exit
      endif
    enddo
    NE3res(MT) = i
    if (i > 2 .and. E3res(MT, i) == E3res(MT, i - 1)) NE3res(MT) = i - 1
   70   if (NE3res(MT) == 0) then
      E3res(MT, 1) = EmineV
      xs3res(MT, 1) = 0.
      E3res(MT, 2) = EHres
      xs3res(MT, 2) = 0.
      NE3res(MT) = 2
    endif
  enddo
!
! Force same energy range for non-elastic cross section
!
  if (NE3res(3) <= 0) then
    NE3res(3) = NE3res(1)
    do i = 1, NE3res(3)
      E3res(3, i) = E3res(1, i)
      xs3res(3, i) = 0.
      do j = 0, NE3res(2)
        eloc(j) = E3res(2, j)
      enddo
      ee = E3res(3, i)
      call locate(eloc, 0, NE3res(2), ee, nen)
      ea = eloc(nen)
      eb = eloc(nen + 1)
      if (ea /= eb) then
        deltaE = (ee - ea) / (eb - ea)
        xs3res(3, i) = xs3res(1, i) - (xs3res(2, nen) + deltaE * (xs3res(2, nen + 1) - xs3res(2, nen)))
      else
        xs3res(3, i) = xs3res(1, i) - xs3res(2, nen)
      endif
    enddo
    xs3res(3, NE3res(3)) = xs3res(1, NE3res(1)) - xs3res(2, NE3res(2))
  endif
  close (unit = 3)
  return
end subroutine readbackground
! Copyright A.J. Koning 2021
