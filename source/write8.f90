subroutine write8
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF8
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
!   nummt       ! number of MT numbers
! Variables for input of ENDF structure
!   flagrp10    ! flag to put residual production cross sections in MF10
! Variables for initialization of ENDF format
!   AWR         ! standard mass parameter
!   blank2     ! blank string
!   FEND       ! ENDF - 6 format
!   LIS         ! state number of target nucleus
!   LISO        ! isomeric state number
!   MAT         ! MAT number
!   mtexist     ! flag for existence of MT - number
!   SEND       ! ENDF - 6 format
!   ZA          ! standard charge parameter
! Variables for MF8_10
!   ELFS        ! excitation energy of final state
!   ErpZAiso    ! energy of isomer
!   IZAP        ! second IZAP - number
!   LFS         ! final state number
!   LFSZA       ! final state number
!   LMF         ! file number for information
!   MATP        ! material number for reaction product
!   NDec        ! number of decay branches
!   NO          ! flag for decay information
!   NSt         ! number of final states
!   ZAPr        ! designation of final nucleus
!
! *** Declaration of local data
!
  implicit none
  integer :: ist              ! counter for final state
  integer :: iza              ! counter for Z,A combinations
  integer :: MF               ! MF-number
  integer :: MT               ! MT-number
  integer :: NS               ! line number
!
! ***************************** Write MF8 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 8
  NS = 0
  open (unit = 2, file = 'MF8', status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(8, MT)) cycle
    call hrwrite(ZA, AWR, LIS, LISO, NSt(MT), NO, MAT, MF, MT, NS)
    if (flagrp10 .and. MT == 5) then
      do iza = 1, NSt(MT)
        call hrwrite(real(IZAP(iza)), ErpZAiso(iza), 10, LFSZA(iza), 0, 0, MAT, MF, MT, NS)
      enddo
    else
      do ist = 1, NSt(MT)
        call hrwrite(ZAPr(MT, ist), ELFS(MT, ist), LMF(MT, ist), LFS(MT, ist), 6 * NDec, MATP(MT, ist), MAT, MF, MT, NS)
      enddo
    endif
    write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  enddo
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write8
! Copyright A.J. Koning 2021
