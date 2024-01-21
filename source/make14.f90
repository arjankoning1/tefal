subroutine make14
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF14
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
!   nummt      ! number of MT numbers
! Variables for initialization of ENDF format
!   blank2     ! blank string
!   FEND       ! ENDF - 6 format
!   MAT        ! MAT number
!   mfexist    ! flag for existence of MF - number
!   mtexist    ! flag for existence of MT - number
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF12_15
!   LI14       ! isotropy flag
!
! *** Declaration of local data
!
  implicit none
  integer :: MF              ! MF-number
  integer :: MT              ! MT-number
  integer :: NS              ! line number
!
! ***************************** Make MF14 ******************************
!
! write14: subroutine to write MF14
!
! By default, all photon angular distributions are assumed to be isotropic.
!
  MF = 14
  NS = 0
  open (unit = 2, file = 'MF14', status = 'replace')
  do MT = 1, nummt
    if (mtexist(12, MT) .or. mtexist(13, MT)) then
      mtexist(MF, MT) = .true.
      LI14(MT) = 1
      if (mtexist(12, MT)) NK(MF, MT) = NK(12, MT)
      if (mtexist(13, MT)) NK(MF, MT) = NK(13, MT)
      mfexist(MF) = .true.
      call write14(MT)
    else
      mtexist(MF, MT) = .false.
    endif
  enddo
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine make14
! Copyright A.J. Koning 2021
