subroutine processtotal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process total cross sections
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
!   sgl            ! single precision kind
!   dbl            ! double precision kind
! All global variables
!   numchan        ! maximum number of exclusive channels
!   numenin        ! number of incident energies
!   nummt          ! number of MT numbers
! Variables for ENDF limits, switches and tolerances
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
! Variables for input of ENDF structure
!   flagngn        ! flag to include (n, gamma n) data
! Variables for input of ENDF library type
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
! Variables for input of specific ENDF data
!   urrenergy      ! upper energy of the URR in MeV
! Variables for info from TALYS
!   eninc          ! incident energy
!   flagelastic    ! flag to give priority to elastic cross section for normalization
!   k0             ! index of incident particle
!   numinc         ! number of incident energies
! Constants
!   xsepslow       ! lower limit for cross sections in millibarns
! Variables for initialization of ENDF format
!   MTid           ! channel identifier for MT - number
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   idchannel      ! identifier for channel
! Variables for total cross sections in ENDF format
!   eninc6         ! incident energy for total cross section
!   numinc6        ! number of incident energies for total cross sections
!   xsany          ! (x, anything) cross section (MT5)
!   xsany6         ! (x, anything) cross section (MT5)
!   xselas6        ! total elastic cross section
!   xsnonel        ! nonelastic cross section
!   xsnonel6       ! nonelastic cross section
!   xstot6         ! total cross section
! Variables for partial cross sections in ENDF format
!   Ethexcl        ! threshold energy
!   xsexcl         ! exclusive cross section
! Variables for discrete state cross sections in ENDF format
!   edisc          ! incident energy for discrete level cross section
!   Ethdisc        ! threshold energy
!   numendisc      ! number of incident energies including discrete energy grid
!   Qdisc          ! Q - value
!   xsbin          ! binary cross section
!   xsngn          ! (projectile, gamma - ejectile) cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: nen                            ! energy counter
  integer   :: nin                            ! counter for incident energy
  integer   :: type                           ! particle type
  integer   :: i                              ! counter
  integer   :: id                             ! counter for deuterons
  integer   :: idc                            ! help variable
  integer   :: MT                             ! MT-number
  real(sgl) :: deltaE                         ! help variable
  real(sgl) :: ea                             ! help variable
  real(sgl) :: ealog                          ! help variable
  real(sgl) :: eb                             ! help variable
  real(sgl) :: eblog                          ! help variable
  real(sgl) :: ee                             ! energy
  real(sgl) :: eelog                          ! help variable
  real(sgl) :: eloc(0:numenin+150)            ! help variable
  real(sgl) :: Eth                            ! help variable
  real(dbl) :: xs1                            ! help variable
  real(dbl) :: xs1log                         ! help variable
  real(dbl) :: xs2                            ! help variable
  real(dbl) :: xs2log                         ! help variable
  real(dbl) :: xshigh                         ! high energy cross section
  real(dbl) :: xsMT                           ! MT cross section
  real(dbl) :: xsngx                          ! (projectile,gamma-x) cross section
  real(dbl) :: xsnonMT                        ! non-MT cross section
  real(dbl) :: xsr                            ! help variable
  real(dbl) :: xss                            ! help variable
!
! ************************* Set URR energy *****************************
!
  if (urrenergy ==  -1.) then
    if (flagfission) then
      if (Ethdisc(k0, 2) > 0.) then
        urrenergy = Ethdisc(k0, 2)
      else
        urrenergy = 0.1
      endif
    else
      if (Ethdisc(k0, 1) > 0.) then
        urrenergy = Ethdisc(k0, 1)
      else
        urrenergy = 0.1
      endif
    endif
  endif
!
! ************* Interpolation of total cross sections ******************
!
! locate        : subroutine to find value in ordered table
!
  if (k0 == 1 .or. (k0 /= 1 .and. flagendfdet)) then
    do nin = 1, numinc6
      xsany6(nin) = 0.
      ee = eninc6(nin)
      if (ee > Eswitch .or. ee > eninc(numinc)) then
        xsany6(nin) = xsnonel6(nin)
      else
        call locate(eninc, 1, numinc, ee, nen)
        ea = eninc(nen)
        eb = eninc(nen + 1)
        deltaE = (ee - ea) / (eb - ea)
!
! ngn etc. If no (n,gamma n) is included, the difference will end up in the elastic cross section.
!
        xsngx = 0.
        if (flagngn) then
          do type = 1, 6
            xsngx = xsngx + xsngn(type, nen) + deltaE * (xsngn(type, nen + 1) - xsngn(type, nen))
          enddo
        endif
!
! Exclusive channels except binary particle channels
!
        xsr = xsngx
        xsMT = 0.
        xsnonMT = 0.
Loop1:  do idc = -1, numchan
          id = idchannel(idc)
          if (id == 100000) cycle
          if (id == 10000) cycle
          if (id == 1000) cycle
          if (id == 100) cycle
          if (id == 10) cycle
          if (id == 1) cycle
          xs1 = xsexcl(idc, nen)
          xs2 = xsexcl(idc, nen + 1)
          if (xs1 == 0..and.xs2 == 0.) cycle
          if (id == 0 .and. xs1 /= 0..and.xs2 /= 0..and.ea /= eb .and. ee < 0.001) then
            ealog = log(ea)
            eblog = log(eb)
            eelog = log(ee)
            xs1log = log(xs1)
            xs2log = log(xs2)
            deltaE = (eelog - ealog) / (eblog - ealog)
            xss = exp(xs1log + deltaE * (xs2log - xs1log))
          else
            Eth = Ethexcl(idc)
            if (Eth > ea .and. Eth <= eb) then
              if (ee > Eth) then
                deltaE = (ee - Eth) / (eb - Eth)
              else
                deltaE = 0.
              endif
            else
              deltaE = (ee - ea) / (eb - ea)
            endif
            xss = xs1 + deltaE * (xs2 - xs1)
          endif
          xsr = xsr + xss
          do MT = 1, nummt
            if (MTid(MT) == id) then
              xsMT = xsMT + xss
              cycle Loop1
            endif
          enddo
          xsnonMT = xsnonMT + xss
        enddo Loop1
!
! Binary particle channels
!
        do type = 1, 6
          if (ee < edisc(type, 1) .or. ee > edisc(type, numendisc(type))) cycle
          eloc(0) = 0.
          do i = 1, numendisc(type)
            eloc(i) = edisc(type, i)
          enddo
          call locate(eloc, 1, numendisc(type), ee, nen)
          ea = eloc(nen)
          eb = eloc(nen + 1)
          xs1 = xsbin(type, nen)
          xs2 = xsbin(type, nen + 1)
          if (Qdisc(type, 0) > 0..and.xs1 /= 0..and.xs2 /= 0..and. ea /= eb .and. ee < 0.001) then
            ealog = log(ea)
            eblog = log(eb)
            eelog = log(ee)
            xs1log = log(xs1)
            xs2log = log(xs2)
            deltaE = (eelog - ealog) / (eblog - ealog)
            xss = exp(xs1log + deltaE * (xs2log - xs1log))
          else
            if (type == k0) then
              Eth = Ethdisc(type, 1)
            else
              Eth = Ethdisc(type, 0)
            endif
            if (Eth > ea .and. Eth <= eb) then
              if (ee > Eth) then
                deltaE = (ee - Eth) / (eb - Eth)
              else
                deltaE = 0.
              endif
            else
              deltaE = (ee - ea) / (eb - ea)
            endif
            xss = xs1 + deltaE * (xs2 - xs1)
          endif
          xsr = xsr + xss
          xsMT = xsMT + xss
        enddo
!
! Total nonelastic cross section
!
        xshigh = xsnonel6(nin) - xsr
        if (ee <= Eswitch) then
          xsnonel6(nin) = xsr
          xsany6(nin) = xsnonMT + xsngx
        else
          xsany6(nin) = xsnonMT + xshigh + xsngx
        endif
        if (k0 == 1) then
          if (flagelastic) then
            xstot6(nin) = xselas6(nin) + xsnonel6(nin)
          else
            xselas6(nin) = xstot6(nin) - xsnonel6(nin)
          endif
        endif
      endif
    enddo
  else
    do nin = 1, numinc6
      xsany6(nin) = xsnonel6(nin)
    enddo
  endif
  do nin = 1, numinc
    xsany(nin) = 0.
    ee = eninc(nin)
    call locate(eninc6, 1, numinc6, ee, nen)
    if (ee < eninc6(numinc6)) then
      ea = eninc6(nen)
      eb = eninc6(nen + 1)
      xs1 = xsany6(nen)
      xs2 = xsany6(nen + 1)
      deltaE = (ee - ea) / (eb - ea)
      xsany(nin) = xs1 + deltaE * (xs2 - xs1)
    else
      xsany(nin) = xsany6(numinc6)
    endif
    if (flagfission .and. ee > Eswitch .and. .not. flagfis10) then
      xsany(nin) = xsany(nin) - dble(xsexcl(-1, nin))
    endif
  enddo
!
! *********** Determine limits for non-elastic cross section ***********
!
  do nin = 1, numinc
    if (xsnonel(nin) < xsepslow) xsnonel(nin) = 0.
    if (xsany(nin) < xsepslow) xsany(nin) = 0.
  enddo
  do nin = 1, numinc6
    if (xsnonel6(nin) < xsepslow) xsnonel6(nin) = 0.
    if (xsany6(nin) < xsepslow) xsany6(nin) = 0.
  enddo
  return
end subroutine processtotal
! Copyright A.J. Koning 2021
