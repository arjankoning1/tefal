subroutine processdiscrete
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process discrete level and continuum data
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
!   numlevels      ! maximum number of discrete levels
! Variables for input of ENDF structure
!   flagngn        ! flag to include (n, gamma n) data
! Variables for ENDF limits, switches and tolerances
!   disclim        ! limit for specific MT numbers for discrete levels
! Variables for info from TALYS
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for reaction initialization
!   nlevmax        ! number of included discrete levels
!   numcut         ! number of energies before high - energy format
! Constants
!   xsepshigh      ! upper limit for cross sections in millibarns
!   xsepslow       ! lower limit for cross sections in millibarns
! Variables for discrete state cross sections in ENDF format
!   edisc          ! incident energy for discrete level cross section
!   Ethdisc        ! threshold energy
!   levexist       ! flag for existence of discrete level
!   ndisc          ! number of discrete levels
!   numendisc      ! number of incident energies including discrete energy grid
!   Qdisc          ! Q - value
!   xsbin          ! binary cross section
!   xscont         ! continuum cross section
!   xsdisc         ! discrete state cross section
!   xsintdisc      ! interpolated discrete cross section
!   xsngn          ! (projectile, gamma - ejectile) cross section
!   xsngnsum       ! total (projectile, gamma - ejectile) cross section
!   yieldngn       ! (projectile, gamma - ejectile) particle yield
! Variables for photon production in ENDF format
!   branchlevel    ! level to which branching takes place
!   Nbranch        ! number of branches for level
!   nlev           ! number of excited levels for nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                 ! counter
  integer   :: ie                ! counter
  integer   :: j                 ! counter
  logical   :: limit             ! help variable
  integer   :: nen               ! energy counter
  integer   :: nex               ! discrete level
  integer   :: nex2              ! counter
  integer   :: nexbeg            ! start of discrete levels numbering
  integer   :: nin               ! counter for incident energy
  integer   :: type              ! particle type
  real(sgl) :: deltaE            ! help variable
  real(sgl) :: ea                ! help variable
  real(sgl) :: ealog             ! help variable
  real(sgl) :: eb                ! help variable
  real(sgl) :: eblog             ! help variable
  real(sgl) :: ee                ! energy
  real(sgl) :: eelog             ! help variable
  real(sgl) :: Eth               ! help variable
  real(sgl) :: etmp              ! help variable
  real(sgl) :: xs1               ! help variable
  real(sgl) :: xs1log            ! help variable
  real(sgl) :: xs2               ! help variable
  real(sgl) :: xs2log            ! help variable
  real(sgl) :: xsc               ! interpolated cross section
!
! ** Determine limits for discrete state and continuum cross sections **
!
  do type = 1, 6
    do nex = 0, nlevmax
      limit = .false.
      do nin = 1, numcut
        if (xsdisc(type, nex, nin) < xsepslow) then
          if (Qdisc(type, nex) > 0.) then
            xsdisc(type, nex, nin) = xsepslow
          else
            xsdisc(type, nex, nin) = 0.
          endif
        endif
        if (xsdisc(type, nex, nin) >= xsepshigh) limit = .true.
      enddo
      if ( .not. limit) then
        levexist(type, nex) = .false.
        do nin = 1, numcut
          xsdisc(type, nex, nin) = 0.
        enddo
      endif
    enddo
    limit = .false.
    do nin = 1, numcut
      if (xscont(type, nin) < xsepslow) then
        if (Qdisc(type, ndisc(type)) > 0.) then
          xscont(type, nin) = xsepslow
        else
          xscont(type, nin) = 0.
        endif
      endif
      if (xscont(type, nin) >= xsepshigh) limit = .true.
    enddo
    if ( .not. limit) then
      do nin = 1, numcut
        xscont(type, nin) = 0.
      enddo
    endif
!
! Lump discrete levels into continuum in case of too small cross section (we put the limit at 10 mb).
! We always include 4 discrete levels.
!
    if (type == k0) then
      nexbeg = 1
    else
      nexbeg = -1
    endif
Loop1: do nex = nlevmax, 0, -1
      do nin = 1, numcut
        if (xsdisc(type, nex, nin) > disclim) then
          nexbeg = nex
          exit Loop1
        endif
      enddo
    enddo Loop1
    if (type == k0) then
      nexbeg = max(nexbeg, 2)
    else
      nexbeg = max(nexbeg, 1)
    endif
    do nex = nexbeg + 1, nlevmax
      levexist(type, nex) = .false.
      do nin = 1, numcut
        xscont(type, nin) = xscont(type, nin) + xsdisc(type, nex, nin)
        xsdisc(type, nex, nin) = 0.
      enddo
    enddo
  enddo
!
! ************* (n,gn), (n,gp)...,(n,galpha) cross sections ************
!
  if (flagngn) then
    do nin = 1, numcut
      xsngnsum(nin) = 0.
    enddo
    do type = 1, 6
      limit = .false.
      do nin = 1, numcut
        yieldngn(type, nin) = 0.
        if (xsngn(type, nin) < xsepslow) xsngn(type, nin) = 0.
        if (xsngn(type, nin) >= xsepshigh) limit = .true.
        xsngnsum(nin) = xsngnsum(nin) + xsngn(type, nin)
      enddo
      if ( .not. limit) then
        do nin = 1, numcut
          xsngn(type, nin) = 0.
        enddo
      endif
    enddo
    do type = 1, 6
      do nin = 1, numcut
        if (xsngnsum(nin) < xsepslow) xsngnsum(nin) = 0.
        if (xsngnsum(nin) /= 0.) yieldngn(type, nin) = xsngn(type, nin) / xsngnsum(nin)
      enddo
    enddo
  endif
!
! ************* Unification of discrete and total energy points ********
!
  do type = 1, 6
    ie = 0
    do nin = 1, numcut
      ie = ie + 1
      edisc(type, ie) = eninc(nin)
    enddo
    do nex = 0, nlevmax
      if (type == k0 .and. nex == 0) cycle
      if (Ethdisc(type, nex) == 0.) cycle
      ie = ie + 1
      edisc(type, ie) = Ethdisc(type, nex)
    enddo
    numendisc(type) = ie
    do i = 1, numendisc(type)
      do j = 1, i
        if (edisc(type, i) < edisc(type, j)) then
          etmp = edisc(type, j)
          edisc(type, j) = edisc(type, i)
          edisc(type, i) = etmp
        endif
      enddo
    enddo
!
! Interpolation
!
    do i = 1, numendisc(type)
      xsbin(type, i) = 0.
      ee = edisc(type, i)
      do nex = 0, nlevmax
        xsintdisc(type, nex, i) = 0.
        if (ee >= eninc(numcut)) xsintdisc(type, nex, i) = xsdisc(type, nex, numcut)
      enddo
      if (ee >= eninc(numcut)) cycle
      call locate(eninc, 1, numcut, ee, nen)
      ea = eninc(nen)
      eb = eninc(nen + 1)
      do nex = 0, nlevmax
        if (type == k0 .and. nex == Ltarget) cycle
        xs1 = xsdisc(type, nex, nen)
        xs2 = xsdisc(type, nex, nen + 1)
        if (Qdisc(type, 0) > 0..and.xs1 /= 0..and.xs2 /= 0.) then
          ealog = log(ea)
          eblog = log(eb)
          eelog = log(ee)
          xs1log = log(xs1)
          xs2log = log(xs2)
          deltaE = (eelog - ealog) / (eblog - ealog)
          xsintdisc(type, nex, i) = exp(xs1log + deltaE * (xs2log - xs1log))
        else
          Eth = Ethdisc(type, nex)
          if (Eth > ea .and. Eth <= eb) then
            if (ee > Eth) then
              deltaE = (ee - Eth) / (eb - Eth)
            else
              deltaE = 0.
            endif
          else
            deltaE = (ee - ea) / (eb - ea)
          endif
          xsintdisc(type, nex, i) = xs1 + deltaE * (xs2 - xs1)
        endif
      enddo
    enddo
    do nex = 0, nlevmax
      do i = 1, numendisc(type)
        xsbin(type, i) = xsbin(type, i) + xsintdisc(type, nex, i)
      enddo
    enddo
    do i = 1, numendisc(type)
      ee = edisc(type, i)
      if (ee >= eninc(numcut)) then
        xsbin(type, i) = xsbin(type, i) + xscont(type, numcut)
        cycle
      endif
      call locate(eninc, 1, numcut, ee, nen)
      ea = eninc(nen)
      eb = eninc(nen + 1)
      xs1 = xscont(type, nen)
      xs2 = xscont(type, nen + 1)
      if (Qdisc(type, 0) > 0..and.xs1 /= 0..and.xs2 /= 0.) then
        ealog = log(ea)
        eblog = log(eb)
        eelog = log(ee)
        xs1log = log(xs1)
        xs2log = log(xs2)
        deltaE = (eelog - ealog) / (eblog - ealog)
        xsc = exp(xs1log + deltaE * (xs2log - xs1log))
      else
        Eth = Ethdisc(type, ndisc(type))
        if (Eth > ea .and. Eth <= eb) then
          if (ee > Eth) then
            deltaE = (ee - Eth) / (eb - Eth)
          else
            deltaE = 0.
          endif
        else
          deltaE = (ee - ea) / (eb - ea)
        endif
        xsc = xs1 + deltaE * (xs2 - xs1)
      endif
      xsbin(type, i) = xsbin(type, i) + xsc
    enddo
  enddo
!
! Correct branching ratios for levels which have been left out
!
  do type = 1, 6
    do nex = 0, nlevmax
      if ( .not. levexist(type, nex)) then
        do i = nex + 1, min(nlev(type), numlevels)
Loop2:    do j = 1, Nbranch(type, i)
            if (branchlevel(type, i, j) == nex) then
              do nex2 = nex - 1, 0, -1
                if (levexist(type, nex2)) then
                  branchlevel(type, i, j) = nex2
                  cycle Loop2
                endif
              enddo
            endif
          enddo Loop2
        enddo
      endif
    enddo
  enddo
  return
end subroutine processdiscrete
! Copyright A.J. Koning 2021
