subroutine make3
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3
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
! Variables for input of ENDF library type
!   flagclean      ! flag to clean up double points
!   flageaf        ! flag for EAF - formatted activation library
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
! Variables for info from TALYS
!   k0             ! index of incident particle
!
! *** Declaration of local data
!
  implicit none
!
! **************************** Make MF3 ********************************
!
! read3         : subroutine to read MF3 from existing ENDF-6 data library
! readbackground: subroutine to read MF3 background from existing ENDF-6 data library
! make3total    : subroutine to make MF3 for total, elastic and reaction cross sections
! make3discrete : subroutine to make MF3 for discrete levels and continuum
! make3mt5      : subroutine to make MF3 for MT5
! make3partial  : subroutine to make MF3 for partial cross sections
! make3eaf      : subroutine to make MF3 for isomeric cross sections in EAF-format
! make3clean    : subroutine to clean up MF3
!
  call read3
  if (flaggpf) then
    if (k0 == 1) call readbackground
    call make3total
    if (flagendfdet) call make3discrete
  endif
  if ( .not. flageaf) call make3mt5
  call make3partial
  if (flageaf) call make3eaf
  if (flagclean) call make3clean
  return
end subroutine make3
! Copyright A.J. Koning 2021
