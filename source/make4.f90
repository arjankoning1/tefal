subroutine make4
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF4
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
!   nummt        ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
! Variables for input of ENDF library type
!   flagclean    ! flag to clean up double points
! Variables for input of ENDF structure
!   flagdisc6    ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
! Variables for initialization of ENDF format
!   blank2       ! blank string
!   FEND         ! ENDF - 6 format
!   MAT          ! MAT number
!   mtexist      ! flag for existence of MT - number
! Variables for MF4
!   NE           ! number of incident energies
!   NE4r         ! number of incident energies (MF4 only)
!   NEh          ! number of incident energies (MF4 only)
!   NEhr         ! number of incident energies (MF4 only)
!
! *** Declaration of local data
!
  implicit none
  integer :: MF              ! MF-number
  integer :: MT              ! MT-number
!
! **************************** Make MF4 ********************************
!
! make4elastic : subroutine to make MF4 for elastic scattering
! read4        : subroutine to read MF4 from existing ENDF-6 data library
! make4discrete: subroutine to make MF4 for discrete levels
! make4fission : subroutine to make MF4 for fission
! make4clean   : subroutine to clean up MF4
! write4       : subroutine to write MF4
!
  MF = 4
  open (unit = 2, file = 'MF4', status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(3, MT)) cycle
    NE = 0
    NEh = 0
    NE4r = 0
    NEhr = 0
    if (adopt(MF, MT)) call read4(MT)
    if (MT == 2) then
      call make4elastic
    else
      if ( .not. flagdisc6) call make4discrete(MT)
      if (MT == 18) call make4fission
    endif
    if (mtexist(MF, MT)) then
      if (flagclean) call make4clean(MT)
      call write4(MT)
     endif
  enddo
  write(2, fmt = FEND) blank2, MAT, 0, 0, 0
  close (unit = 2)
  return
end subroutine make4
! Copyright A.J. Koning 2021
