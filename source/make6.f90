subroutine make6
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF6
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
!   nummt          ! number of MT numbers
! Variables for input of ENDF structure
!   flagcapt6      ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagdisc6      ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
! Variables for input of ENDF library type
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
! Variables for info from TALYS
!   k0             ! index of incident particle
! Variables for initialization of ENDF format
!   MAT            ! MAT number
!   mtexist        ! flag for existence of MT - number
!   blank2        ! blank string
!   FEND           ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  integer :: MT              ! MT-number
!
! **************************** Make MF6 ********************************
!
! make6mt2     : subroutine to make MF6 for MT2
! make6mt5     : subroutine to make MF6 for MT5
! make6partial : subroutine to make MF6 for partial cross sections
! make6discrete: subroutine to make MF6 for discrete levels
!
  open (unit = 2, file = 'MF6', status = 'replace')
  if (k0 > 1 .and. flaggpf) then
    call make6mt2
    if (flagendfdet .and. mtexist(3, 4)) call make6partial(4)
  endif
  call make6mt5(5)
  if (flagendfdet .and. flaggpf) then
    do MT = 6, nummt
      if ( .not. mtexist(3, MT)) cycle
      if (mtexist(4, MT)) cycle
      if (mtexist(5, MT)) cycle
      if (flagdisc6) call make6discrete(MT)
      if (MT /= 102) call make6partial(MT)
      if (flagcapt6 .and. MT == 102) call make6partial(MT)
      if (k0 /= 1 .and. MT == 18) call make6mt5(MT)
    enddo
  endif
  write(2, fmt = FEND) blank2, MAT, 0, 0, 0
  close (unit = 2)
  return
end subroutine make6
! Copyright A.J. Koning 2021
