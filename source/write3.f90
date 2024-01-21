subroutine write3
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF3
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
! Definition of single and double precision variables
!   sgl          ! single precision kind
! All global variables
!   numen6       ! number of incident energies
!   numint       ! number of interpolation sections
!   nummt        ! number of MT numbers
! Variables for input of ENDF library type
!   flageaf      ! flag for EAF - formatted activation library
! Variables for info from TALYS
!   Lisomer      ! isomeric number of target
! Variables for initialization of ENDF format
!   AWR          ! standard mass parameter
!   blank2       ! blank string
!   FEND         ! ENDF - 6 format
!   MAT          ! MAT number
!   mtexist      ! flag for existence of MT - number
!   SEND         ! ENDF - 6 format
!   TEXT         ! ENDF - 6 format
!   ZA           ! standard charge parameter
! Variables for ENDF format
!   INTER        ! interpolation scheme
!   NBT          ! separation value for interpolation scheme
!   NP           ! number of incident energies
!   NR           ! number of interpolation ranges
! Variables for MF3
!   eafstring    ! string with reaction information for EAF format
!   E3           ! incident energy for MF3 (in ENDF - 6 format)
!   LFS3         ! isomeric state number (EAF only)
!   LR3          ! breakup flag
!   mtstring     ! string with reaction information
!   QI           ! Q - value (in ENDF - 6 format)
!   QM           ! Q - value (in ENDF - 6 format)
!   xs           ! cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: ii                        ! counter
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(2*numen6)               ! help variable
!
! ***************************** Write MF3 ******************************
!
! hrwrite  : subroutine to write header with real values
!
  MF = 3
  NS = 0
  open (unit = 2, file = 'MF3', status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(3, MT)) cycle
!
! Special header of section for EAF-format
!
    if (flageaf) then
      write(2, fmt = TEXT) eafstring(MT), MAT, MF, MT, NS
      write(2, fmt = TEXT) mtstring(MT), MAT, MF, MT, NS
      call hrwrite(ZA, AWR, Lisomer, LFS3(MT), 0, 0, MAT, MF, MT, NS)
      call hrwrite(0., QI(MT), 0, LR3(MT), NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
    else
!
! Header for normal ENDF-6 format
!
      call hrwrite(ZA, AWR, 0, 0, 0, 0, MAT, MF, MT, NS)
      call hrwrite(QM(MT), QI(MT), 0, LR3(MT), NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
    endif
!
! Write interpolation ranges
!
! kwrite: subroutine to write integer value block
!
    do i = 1, NR(MF, MT)
      ii = 2 * i - 1
      Nval(ii) = NBT(MF, MT, i)
      Nval(ii + 1) = INTER(MF, MT, i)
    enddo
    N = 2 * NR(MF, MT)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write cross sections
!
! xwrite: subroutine to write real value block
!
    do i = 1, NP(MF, MT)
      ii = 2 * i - 1
      x(ii) = E3(MT, i)
      x(ii + 1) = xs(MT, i)
    enddo
    N = 2 * NP(MF, MT)
    call xwrite(N, x, MAT, MF, MT, NS)
    write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  enddo
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write3
! Copyright A.J. Koning 2021
