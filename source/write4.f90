subroutine write4(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF4
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
!   sgl       ! single precision kind
! All global variables
!   numang    ! number of angles
!   numint    ! number of interpolation sections
! Variables for initialization of ENDF format
!   AWR       ! standard mass parameter
!   blank2    ! blank string
!   MAT       ! MAT number
!   SEND      ! ENDF - 6 format
!   ZA        ! standard charge parameter
! Variables for ENDF format
!   INTER     ! interpolation scheme
!   NBT       ! separation value for interpolation scheme
!   NR        ! number of interpolation ranges
! Variables for MF4
!   E4        ! incident energy for MF4 (in ENDF - 6 format)
!   E4h       ! incident energy for MF4 (in ENDF - 6 format)
!   f4        ! angular distribution
!   INTER4    ! interpolation scheme
!   INTERh    ! interpolation scheme
!   LCT       ! LAB / CM flag
!   leg       ! Legendre coefficients (in ENDF - 6 format)
!   LI4       ! isotropy flag
!   LTT       ! representation
!   LVT       ! specification of transformation matrix
!   NBT4      ! separation value for interpolation scheme
!   NBTh      ! separation value for interpolation scheme
!   NE        ! number of incident energies
!   NEh       ! number of incident energies (MF4 only)
!   NL        ! Legendre order or number of cosines
!   NP4       ! number of incident energies
!   NR4       ! number of interpolation ranges
!   NRh       ! number of interpolation ranges
!   x4        ! cosine of the angle
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                          ! counter
  integer   :: iE                         ! energy counter
  integer   :: ii                         ! counter
  integer   :: MF                         ! MF-number
  integer   :: MT                         ! MT-number
  integer   :: N                          ! neutron number of residual nucleus
  integer   :: NS                         ! line number
  integer   :: Nval(2*numint)             ! value
  real(sgl) :: x(2*(numang+1))            ! help variable
!
! ***************************** Write MF4 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 4
  NS = 0
  call hrwrite(ZA, AWR, LVT, LTT, 0, 0, MAT, MF, MT, NS)
  call hrwrite(0., AWR, LI4, LCT, 0, 0, MAT, MF, MT, NS)
!
! ********************** Legendre coefficients *************************
!
! kwrite: subroutine to write integer value block
!
! Write interpolation ranges
!
  if (LTT == 1 .or. LTT == 3) then
    call hrwrite(0., 0., 0, 0, NR(MF, MT), NE, MAT, MF, MT, NS)
    do i = 1, NR(MF, MT)
      ii = 2 * i - 1
      Nval(ii) = NBT(MF, MT, i)
      Nval(ii + 1) = INTER(MF, MT, i)
    enddo
    N = 2 * NR(MF, MT)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write Legendre coefficients
!
! xwrite: subroutine to write real value block
!
    do iE = 1, NE
      N = NL(MF, MT, iE)
      call hrwrite(0., E4(iE), 0, 0, N, 0, MAT, MF, MT, NS)
      do i = 1, N
        x(i) = leg(iE, i)
      enddo
      call xwrite(N, x, MAT, MF, MT, NS)
    enddo
  endif
!
! *********** Tabulated angular distributions (high energies) **********
!
  if (LTT >= 2) then
    call hrwrite(0., 0., 0, 0, NRh, NEh, MAT, MF, MT, NS)
!
! Write interpolation ranges
!
    do i = 1, NRh
      ii = 2 * i - 1
      Nval(ii) = NBTh(i)
      Nval(ii + 1) = INTERh(i)
    enddo
    N = 2 * NRh
    call kwrite(N, Nval, MAT, MF, MT, NS)
  endif
  if (LTT > 0) then
!
! Write angular distributions
!
    do iE = 1, NEh
      call hrwrite(0., E4h(iE), 0, 0, NR4(iE), NP4(iE), MAT, MF, MT, NS)
!
! 1. Interpolation ranges
!
      do i = 1, NR4(iE)
        ii = 2 * i - 1
        Nval(ii) = NBT4(iE, i)
        Nval(ii + 1) = INTER4(iE, i)
      enddo
      N = 2 * NR4(iE)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! 2. Values
!
      do i = 1, NP4(iE)
        ii = 2 * i - 1
        x(ii) = x4(iE, i)
        x(ii + 1) = f4(iE, i)
      enddo
      N = 2 * NP4(iE)
      call xwrite(N, x, MAT, MF, MT, NS)
    enddo
  endif
!
! end of section
!
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write4
! Copyright A.J. Koning 2021
