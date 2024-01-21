subroutine write5(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF5
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
!   sgl         ! single precision kind
! All global variables
!   numen2      ! number of emission energies
!   numint      ! number of interpolation sections
! Variables for initialization of ENDF format
!   AWR         ! standard mass parameter
!   MAT         ! MAT number
!   blank2      ! blank string
!   SEND        ! ENDF - 6 format
!   ZA          ! standard charge parameter
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF1
!   EMAX        ! upper limit of energy range for evaluation
! Variables for MF5
!   E5          ! incident energy for MF5 (in ENDF - 6 format)
!   E5p         ! incident energy for which tabulated distribution is given
!   EFH         ! constant used in the energy - dependent fission neutron spectrum
!   EFL         ! constant used in the energy - dependent fission neutron spectrum
!   gE5         ! energy - spectrum values
!   INTER5      ! interpolation scheme
!   INTER5e     ! interpolation scheme
!   INTER5e2    ! interpolation scheme
!   LF          ! flag for energy distribution law
!   NBT5        ! separation value for interpolation scheme
!   NBT5e       ! separation value for interpolation scheme
!   NBT5e2      ! separation value for interpolation scheme
!   NE5e        ! number of incident energies for distribution
!   NF          ! number of secondary energy points
!   NP5         ! number of incident energies
!   NR5         ! number of interpolation ranges
!   NR5e        ! number of interpolation ranges
!   NR5e2       ! number of interpolation ranges
!   pE          ! fractional part of cross section
!   TM5         ! Effective (7) or maximum (12) temperature FNS parameter, T(E)
!   U           ! constant for upper energy limit
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
  real(sgl) :: x(10*numen2)              ! help variable
!
! ***************************** Write MF5 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 5
  NS = 0
  call hrwrite(ZA, AWR, 0, 0, NK(MF, MT), 0, MAT, MF, MT, NS)
!
! ************************ Loop over subsections ***********************
!
  do k = 1, NK(MF, MT)
    call hrwrite(U(k), 0., 0, LF(k), NR5(k), NP5(k), MAT, MF, MT, NS)
!
! 1. Fractional parts of cross section
!
! Write interpolation ranges
!
! kwrite : subroutine to write integer value block
!
    do i = 1, NR5(k)
      ii = 2 * i - 1
      Nval(ii) = NBT5(k, i)
      Nval(ii + 1) = INTER5(k, i)
    enddo
    N = 2 * NR5(k)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write fractional parts of cross section
!
! xwrite: subroutine to write real value block
!
    do i = 1, NP5(k)
      ii = 2 * i - 1
      x(ii) = E5p(k, i)
      x(ii + 1) = pE(k, i)
    enddo
    N = 2 * NP5(k)
    call xwrite(N, x, MAT, MF, MT, NS)
!
! 2. Specific distributions
!
    if (LF(k) == 1) then
!
! LF=1: arbitrary tabulated function
!
      call hrwrite(0., 0., 0, 0, NR5e(k), NE5e(k), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
      do i = 1, NR5e(k)
        ii = 2 * i - 1
        Nval(ii) = NBT5e(k, i)
        Nval(ii + 1) = INTER5e(k, i)
      enddo
      N = 2 * NR5e(k)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Energy distributions
!
      do iE = 1, NE5e(k)
        call hrwrite(0., E5(k, iE), 0, 0, NR5e2(k, iE), NF(k, iE), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
        do i = 1, NR5e2(k, iE)
          ii = 2 * i - 1
          Nval(ii) = NBT5e2(k, iE, i)
          Nval(ii + 1) = INTER5e2(k, iE, i)
        enddo
        N = 2 * NR5e2(k, iE)
        call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write energy distributions
!
        N = 2 * NF(k, iE)
        do i = 1, N
          x(i) = gE5(k, iE, i)
        enddo
        call xwrite(N, x, MAT, MF, MT, NS)
      enddo
    endif
    if (LF(k) == 5 .or. LF(k) == 7 .or. LF(k) == 9 .or. LF(k) == 12) then
      if (LF(5) == 5) then
        ii = 6
        call hrwrite(0., 0., ii, ii, 1, 2, MAT, MF, MT, NS)
        Nval(1) = 2
        Nval(2) = 2
        call kwrite(2, Nval, MAT, MF, MT, NS)
        x(1) = EmineV
        x(2) = 1.
        x(3) = EMAX
        x(4) = 1.
        call xwrite(4, x, MAT, MF, MT, NS)
      endif
      call hrwrite(EFL, EFH, 0, 0, NR5e(k), NE5e(k), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
      do i = 1, NR5e(k)
        ii = 2 * i - 1
        Nval(ii) = NBT5e(k, i)
        Nval(ii + 1) = INTER5e(k, i)
      enddo
      N = 2 * NR5e(k)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write temperatures
!
      N = 2 * NE5e(k)
      do i = 1, N
        x(i) = TM5(k, i)
      enddo
      call xwrite(N, x, MAT, MF, MT, NS)
    endif
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write5
! Copyright A.J. Koning 2021
