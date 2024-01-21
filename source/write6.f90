subroutine write6(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF6
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
!   numen2        ! number of emission energies
!   numint        ! number of interpolation sections
! Variables for initialization of ENDF format
!   AWR           ! standard mass parameter
!   blank2        ! blank string
!   MAT           ! MAT number
!   SEND          ! ENDF - 6 format
!   ZA            ! standard charge parameter
! Variables for ENDF format
!   NK            ! number of subsections
! Variables for MF4
!   LCT           ! LAB / CM flag
!   leg           ! Legendre coefficients (in ENDF - 6 format)
!   NL            ! Legendre order or number of cosines
! Variables for MF6
!   AWP           ! product mass
!   b6            ! energy - angle values
!   b6gam         ! energy - angle values for photons
!   b6rec         ! energy - angle values for recoils
!   E6            ! incident energy (in ENDF - 6 format) for distribution
!   Ey            ! incident energy for yields (in ENDF - 6 format)
!   flagrec       ! flag to state that b6 element concerns recoil grid
!   INTER6ea      ! interpolation scheme
!   INTER6y       ! interpolation scheme
!   kpart         ! section number for particles
!   LANG          ! flag for angular representation
!   LAW           ! flag for distribution function
!   LEP           ! interpolation scheme for secondary energy
!   LIDP          ! indicates that particles are identical
!   LIP           ! product modifier flag
!   LTP           ! representation
!   NA            ! number of angular parameters
!   NBT6ea        ! separation value for interpolation scheme
!   NBT6y         ! separation value for interpolation scheme
!   ND            ! number of discrete energies
!   NE6ea         ! number of incident energies for distribution
!   NEP           ! number of secondary energy points
!   NP6y          ! number of incident energies for yields
!   NR6ea         ! number of interpolation ranges for distribution
!   NR6y          ! number of interpolation ranges for yields
!   NW            ! number of words
!   SPIpar        ! particle spin
!   Y             ! product yield (in ENDF - 6 format)
!   ZAP           ! product identifier
!
! *** Declaration of local data
!
  implicit none
  integer   :: dtype                     ! data type
  integer   :: i                         ! counter
  integer   :: iE                        ! energy counter
  integer   :: ii                        ! counter
  integer   :: k                         ! counter
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(40*numen2)              ! help variable
!
! ***************************** Write MF6 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 6
  NS = 0
  call hrwrite(ZA, AWR, 0, LCT, NK(MF, MT), 0, MAT, MF, MT, NS)
!
! ************************ Loop over subsections ***********************
!
  do k = 1, NK(MF, MT)
    call hrwrite(ZAP(k), AWP(k), LIP(k), LAW(k), NR6y(k), NP6y(k), MAT, MF, MT, NS)
!
! 1. Product yields
!
! Write interpolation ranges
!
! kwrite : subroutine to write integer value block
!
    do i = 1, NR6y(k)
      ii = 2 * i - 1
      Nval(ii) = NBT6y(k, i)
      Nval(ii + 1) = INTER6y(k, i)
    enddo
    N = 2 * NR6y(k)
    call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write product yields
!
! xwrite: subroutine to write real value block
!
    do i = 1, NP6y(k)
      ii = 2 * i - 1
      x(ii) = Ey(k, i)
      x(ii + 1) = Y(k, i)
    enddo
    N = 2 * NP6y(k)
    call xwrite(N, x, MAT, MF, MT, NS)
!
! 2. Energy-angle distributions
!
    if (LAW(k) == 1 .or. LAW(k) == 2 .or. LAW(k) == 5) then
!
! LAW=1: continuum energy-angle distributions
!
      if (LAW(k) == 1) call hrwrite(0., 0., LANG(k), LEP(k), NR6ea(k), NE6ea(k), MAT, MF, MT, NS)
!
! LAW=2: discrete two-body scattering
!
      if (LAW(k) == 2) call hrwrite(0., 0., 0, 0, NR6ea(k), NE6ea(k), MAT, MF, MT, NS)
!
! LAW=5: charged particle-elastic scattering
!
      if (LAW(k) == 5) call hrwrite(SPIpar, 0., LIDP, 0, NR6ea(k), NE6ea(k), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
      do i = 1, NR6ea(k)
        ii = 2 * i - 1
        Nval(ii) = NBT6ea(k, i)
        Nval(ii + 1) = INTER6ea(k, i)
      enddo
      N = 2 * NR6ea(k)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write energy-angle distributions
!
      dtype = 1
      if (MT == 5 .or. MT == 18) then
        if (k > kpart .and. k < NK(MF, MT) .and. flagrec(k)) dtype = 2
        if (k == NK(MF, MT)) dtype = 3
      endif
      do iE = 1, NE6ea(k)
        if (LAW(k) == 1) call hrwrite(0., E6(k, iE), ND(k, iE), NA(k, iE), NW(k, iE), NEP(k, iE), MAT, MF, MT, NS)
        if (LAW(k) == 2) call hrwrite(0., E6(k, iE), LANG(k), 0, NW(k, iE), NL(6, MT, iE), MAT, MF, MT, NS)
        if (LAW(k) == 5) call hrwrite(0., E6(k, iE), LTP, 0, NW(k, iE), NL(6, MT, iE), MAT, MF, MT, NS)
        N = NW(k, iE)
        if (LAW(k) == 2) then
          do i = 1, N
            x(i) = leg(iE, i)
          enddo
          call xwrite(N, x, MAT, MF, MT, NS)
        else
          if (dtype == 1) then
            do i = 1, N
              x(i) = b6(k, iE, i)
            enddo
          endif
          if (dtype == 2) then
            do i = 1, N
              x(i) = b6rec(k, iE, i)
            enddo
          endif
          if (dtype == 3) then
            do i = 1, N
              x(i) = b6gam(iE, i)
            enddo
          endif
          call xwrite(N, x, MAT, MF, MT, NS)
        endif
      enddo
    endif
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write6
! Copyright A.J. Koning 2021
