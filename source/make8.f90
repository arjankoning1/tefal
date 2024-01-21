subroutine make8
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF8
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
!   nummt         ! number of MT numbers
! Variables for input of ENDF structure
!   flagrp10      ! flag to put residual production cross sections in MF10
!   flagrp6       ! flag to put residual production cross sections in MF6
! Constants
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Variables for info from TALYS
!   Ainit         ! mass number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
! Variables for initialization of ENDF format
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
!   MTid          ! channel identifier for MT - number
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF8_10
!   IZAP          ! second IZAP - number
!   LMF           ! file number for information
!   MATP          ! material number for reaction product
!   NDec          ! number of decay branches
!   NO            ! flag for decay information
!   NSt           ! number of final states
!   NZA           ! number of nuclides
!   ZAPi          ! designation of final nucleus
!   ZAPr          ! designation of final nucleus
!
! *** Declaration of local data
!
  implicit none
  integer :: A                   ! mass number of target nucleus
  integer :: ist                 ! counter for final state
  integer :: iyield              ! particle yield
  integer :: MT                  ! MT-number
  integer :: Nix                 ! neutron number index for residual nucleus
  integer :: type                ! particle type
  integer :: Z                   ! charge number of target nucleus
  integer :: Zix                 ! charge number index for residual nucleus
!
! **************************** Make MF8 ********************************
!
  do MT = 1, nummt
    if ( .not. (mtexist(9, MT) .or. mtexist(10, MT))) cycle
    if (MT == 5) then
      if (flagrp10) then
        do ist = 1, NZA
          LMF(MT, ist) = 10
          MATP(MT, ist) = 0
          ZAPr(MT, ist) = real(IZAP(ist))
        enddo
      endif
    else
      if (MT == 18) then
        ZAPi(MT) = -1
        ZAPr(MT, 1) = -1
        MATP(MT, 1) = 0
        LMF(MT, 1) = 10
      else
        Zix = 0
        Nix = 0
        do type = 1, 6
          iyield = mod(MTid(MT), 10 **(7 - type)) / (10 **(6 - type))
          Zix = Zix + iyield * parZ(type)
          Nix = Nix + iyield * parN(type)
        enddo
        Z = Zinit - Zix
        A = Ainit - Zix - Nix
        ZAPi(MT) = 1000 * Z + A
        do ist = 1, NSt(MT)
          if (mtexist(9, MT)) then
            LMF(MT, ist) = 9
          else
            LMF(MT, ist) = 10
          endif
! MATP is considered obsolete, according to CHECKR, therefore we put it equal to zero.
          MATP(MT, ist) = 0
          ZAPr(MT, ist) = 1000. * Z + A
        enddo
      endif
      mtexist(8, MT) = .true.
      mfexist(8) = .true.
    endif
  enddo
  NDec = 0
  NO = 1
!
! Add data for MT5
!
  MT = 5
  if (mtexist(6, MT)) then
    mtexist(8, MT) = .true.
    mfexist(8) = .true.
    if ( .not. flagrp10 .and. flagrp6) NSt(MT) = NK(6, MT) - 1
  endif
  return
end subroutine make8
! Copyright A.J. Koning 2021
