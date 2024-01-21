subroutine make15
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF15
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
!   nummt        ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
! Variables for reaction initialization
!   Egamindex    ! enegy index for gamma cross sections
!   Nengam       ! number of incident energies for gamna c.s.
!   numcut       ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc        ! incident energy
! Constants
!   xsepslow     ! lower limit for cross sections in millibarns
! Variables for initialization of ENDF format
!   blank2       ! blank string
!   FEND         ! ENDF - 6 format
!   idnum        ! number of different exclusive cross sections
!   MAT          ! MAT number
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
!   MTid         ! channel identifier for MT - number
! Variables for partial cross sections in ENDF format
!   Ehist        ! histogram emission energy
!   f0ex         ! energy distribution for exclusive channel
!   idchannel    ! identifier for channel
!   nbeg         ! first outgoing energy
!   nend         ! last outgoing energy
!   Qexcl        ! Q - value
!   xsexcl       ! exclusive cross section
! Variables for photon production in ENDF format
!   yieldgam     ! gamma yield
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF1
!   EMAX         ! upper limit of energy range for evaluation
! Variables for MF3
!   E3           ! incident energy for MF3 (in ENDF - 6 format)
!   EthMT        ! threshold energy
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
  integer   :: ib                ! counter
  integer   :: id                ! counter for deuterons
  integer   :: idc               ! help variable
  integer   :: iE                ! energy counter
  integer   :: k                 ! counter
  integer   :: MF                ! MF-number
  integer   :: MT                ! MT-number
  integer   :: nen               ! energy counter
  integer   :: nenout            ! counter for outgoing energy
  integer   :: nin               ! counter for incident energy
  real(sgl) :: Ein               ! incident energy
!
! **************************** Make MF15 *******************************
!
! read15   : subroutine to read MF15 from existing ENDF-6 data library
! write15  : subroutine to write MF15
!
! Read existing MF15
!
  MF = 15
  open (unit = 2, file = 'MF15', status = 'replace')
  do MT = 1, nummt
    if (MTid(MT) <=  -1) cycle
    if ( .not. mtexist(12, MT) .and. .not. mtexist(13, MT)) cycle
    if (adopt(MF, MT)) then
      call read15(MT)
      goto 100
    endif
    do id = 0, idnum
      if (idchannel(id) == MTid(MT)) then
        idc = id
        exit
      endif
    enddo
!
! Make new MF15
!
! 1. Probabilities
!
    NK(15, MT) = 1
    do k = 1, NK(15, MT)
      NR15(k) = 1
      INTER15(k, 1) = 2
      EPy(k, 1) = E3(MT, 1)
      Pg(k, 1) = 1.
      EPy(k, 2) = EMAX
      Pg(k, 2) = 1.
      NP15(k) = 2
      NBT15(k, 1) = 2
!
! 2. Photon energy distributions
!
! For non-threshold reactions (e.g. MT102), the discrete gamma yields
! for the first energy are not zero.
!
      if (Qexcl(idc) <= 0..or.yieldgam(idc, 1) == 0.) then
        E15(k, 1) = E3(MT, 1)
        iE = 1
        NP15ge(k, iE) = 2
        NBT15ge(k, iE, 1) = 2
        INTER15ge(k, iE, 1) = 2
        NR15ge(k, iE) = 1
        E15ge(k, iE, 1) = 0.
        ge(k, iE, 1) = 1.
        E15ge(k, iE, 2) = 1.
        ge(k, iE, 2) = 0.
      else
        iE = 0
      endif
      do nen = 1, Nengam
        nin = Egamindex(nen)
        if (nin > numcut) cycle
        Ein = eninc(nin)
        if (Ein <= EthMT(MT) .and. Qexcl(idc) <= 0.) cycle
        if (xsexcl(idc, nin) <= xsepslow) cycle
        if (yieldgam(idc, nen) == 0.) cycle
        iE = iE + 1
        E15(k, iE) = Ein * 1.e6
        ib = 1
        E15ge(k, iE, ib) = 0.
        ge(k, iE, ib) = f0ex(idc, nen, 0, nbeg(idc, 0, nen)) * 1.e-6
        do nenout = nbeg(idc, 0, nen), nend(idc, 0, nen)
          ib = ib + 1
          E15ge(k, iE, ib) = Ehist(idc, nen, 0, nenout) * 1.e6
          ge(k, iE, ib) = f0ex(idc, nen, 0, nenout) * 1.e-6
        enddo
        NP15ge(k, iE) = ib
        NBT15ge(k, iE, 1) = ib
        INTER15ge(k, iE, 1) = 2
        NR15ge(k, iE) = 1
      enddo
      if (EMAX > E15(k, iE)) then
        iE = iE + 1
        E15(k, iE) = EMAX
        NP15ge(k, iE) = ib
        NBT15ge(k, iE, 1) = ib
        INTER15ge(k, iE, 1) = 2
        NR15ge(k, iE) = 1
        do ib = 1, NP15ge(k, iE)
          E15ge(k, iE, ib) = E15ge(k, iE - 1, ib)
          ge(k, iE, ib) = ge(k, iE - 1, ib)
        enddo
      endif
      NR15g(k) = 1
      INTER15g(k, 1) = 2
      NE15g(k) = iE
      NBT15g(k, 1) = iE
    enddo
  100   call write15(MT)
    mfexist(MF) = .true.
    mtexist(MF, MT) = .true.
  enddo
  write(2, fmt = FEND) blank2, MAT, 0, 0, 0
  close (unit = 2)
  return
end subroutine make15
! Copyright A.J. Koning 2021
