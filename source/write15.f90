subroutine write15(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF15
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
!   numen2       ! number of emission energies
!   numint       ! number of interpolation sections
! Variables for initialization of ENDF format
!   AWR          ! standard mass parameter
!   blank2       ! blank string
!   MAT          ! MAT number
!   SEND         ! ENDF - 6 format
!   ZA           ! standard charge parameter
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF12_15
!   E15          ! incident energy (in ENDF - 6 format)
!   E15ge        ! secondary energy
!   EPy          ! incident energy for probabilities (in ENDF - 6 format)
!   ge           ! gamma distribution
!   INTER15      ! interpolation scheme
!   INTER15g     ! interpolation scheme
!   INTER15ge    ! interpolation scheme
!   NBT15        ! separation value for interpolation scheme
!   NBT15g       ! separation value for interpolation scheme
!   NBT15ge      ! separation value for interpolation scheme
!   NE15g        ! number of incident energies for distribution
!   NP15         ! number of incident energies
!   NP15ge       ! number of secondary energy point
!   NR15         ! number of interpolation ranges
!   NR15g        ! number of interpolation ranges
!   NR15ge       ! number of interpolation ranges
!   Pg           ! probability (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: iE                        ! energy counter
  integer   :: ii                        ! counter
  integer   :: k                         ! counter
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(3*numen2)               ! help variable
!
! ***************************** Write MF6 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 15
  NS = 0
  call hrwrite(ZA, AWR, 0, 0, NK(MF, MT), 0, MAT, MF, MT, NS)
!
! ************************ Loop over subsections ***********************
!
  do k = 1, NK(MF, MT)
    call hrwrite(0., 0., 0, 1, NR15(k), NP15(k), MAT, MF, MT, NS)
!
! 1. Probabilities
!
! Write interpolation ranges
!
! kwrite : subroutine to write integer value block
!
    do i = 1, NR15(k)
      ii = 2 * i - 1
      Nval(ii) = NBT15(k, i)
      Nval(ii + 1) = INTER15(k, i)
    enddo
    N = 2 * NR15(k)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write probabilities
!
! xwrite: subroutine to write real value block
!
    do i = 1, NP15(k)
      ii = 2 * i - 1
      x(ii) = EPy(k, i)
      x(ii + 1) = Pg(k, i)
    enddo
    N = 2 * NP15(k)
    call xwrite(N, x, MAT, MF, MT, NS)
!
! 2. Photon-energy distribution
!
    call hrwrite(0., 0., 0, 0, NR15g(k), NE15g(k), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
    do i = 1, NR15g(k)
      ii = 2 * i - 1
      Nval(ii) = NBT15g(k, i)
      Nval(ii + 1) = INTER15g(k, i)
    enddo
    N = 2 * NR15g(k)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write photon energy distributions
!
    do iE = 1, NE15g(k)
      call hrwrite(0., E15(k, iE), 0, 0, NR15ge(k, iE), NP15ge(k, iE), MAT, MF, MT, NS)
      do i = 1, NR15ge(k, iE)
        ii = 2 * i - 1
        Nval(ii) = NBT15ge(k, iE, i)
        Nval(ii + 1) = INTER15ge(k, iE, i)
      enddo
      N = 2 * NR15ge(k, iE)
      call kwrite(N, Nval, MAT, MF, MT, NS)
      do i = 1, NP15ge(k, iE)
        ii = 2 * i - 1
        x(ii) = E15ge(k, iE, i)
        x(ii + 1) = ge(k, iE, i)
      enddo
      N = 2 * NP15ge(k, iE)
      call xwrite(N, x, MAT, MF, MT, NS)
    enddo
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write15
! Copyright A.J. Koning 2021
