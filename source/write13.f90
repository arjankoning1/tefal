subroutine write13(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF13
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
!   sgl           ! single precision kind
! All global variables
!   numenin       ! number of incident energies
!   numint        ! number of interpolation sections
!  Variables for partial cross sections in ENDF format
!   idchannel     ! identifier for channel
! Variables for initialization of ENDF format
!   AWR           ! standard mass parameter
!   blank2        ! blank string
!   idnum         ! number of different exclusive cross sections
!   MAT           ! MAT number
!   MTid          ! channel identifier for MT - number
!   SEND          ! ENDF - 6 format
!   ZA            ! standard charge parameter
! Variables for ENDF format
!   INTER         ! interpolation scheme
!   NBT           ! separation value for interpolation scheme
!   NK            ! number of subsections
!   NP            ! number of incident energies
!   NR            ! number of interpolation ranges
! Variables for MF12_15
!   E13           ! incident energy (in ENDF - 6 format)
!   Eg            ! gamma energy
!   Egk           ! gamma energy
!   Esk           ! starting level (in ENDF - 6 format)
!   INTERg        ! interpolation scheme
!   LFg           ! photo energy distribution law
!   LPg           ! primary photon flag
!   NBTg          ! separation value for interpolation scheme
!   NPg           ! number of incident energies
!   xsg           ! gamma - ray cross section (in ENDF - 6 format)
!   xsgtot        ! total discrete photon production cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: id                        ! counter for deuterons
  integer   :: idc                       ! help variable
  integer   :: igam                      ! counter for gammas
  integer   :: ii                        ! counter
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(2*numenin)              ! help variable
!
! ***************************** Write MF13 *****************************
!
! hrwrite  : subroutine to write header with real values
!
  MF = 13
  NS = 0
  do id = 0, idnum
    if (idchannel(id) == MTid(MT)) then
      idc = id
      exit
    endif
  enddo
  call hrwrite(ZA, AWR, 0, 0, NK(MF, MT), 0, MAT, MF, MT, NS)
  if (NK(MF, MT) /= 1) then
    call hrwrite(0., 0., 0, 0, NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
!
! 1. Total discrete gamma-ray production cross sections
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
      x(ii) = E13(idc, i)
      x(ii + 1) = xsgtot(idc, i)
    enddo
    N = 2 * NP(MF, MT)
    call xwrite(N, x, MAT, MF, MT, NS)
!
! 2. Gamma-ray production cross sections per discrete transition
!
  endif
  do igam = 1, NK(MF, MT)
    call hrwrite(Egk(idc, igam), Esk(idc, igam), LPg(MT, igam), LFg(MT, igam), NR(MF, MT), NPg(MT, igam), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
    do i = 1, NR(MF, MT)
      ii = 2 * i - 1
      Nval(ii) = NBTg(MT, igam, i)
      Nval(ii + 1) = INTERg(MT, igam, i)
    enddo
    N = 2 * NR(MF, MT)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write cross sections
!
    do i = 1, NPg(MT, igam)
      ii = 2 * i - 1
      x(ii) = Eg(idc, igam, i)
      x(ii + 1) = xsg(idc, igam, i)
    enddo
    N = 2 * NPg(MT, igam)
    call xwrite(N, x, MAT, MF, MT, NS)
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write13
! Copyright A.J. Koning 2021
