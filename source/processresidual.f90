subroutine processresidual
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process residual production data
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
! All global variables
!   numN           ! maximal number of neutrons from initial compound nucleus
!   numrp          ! nunber of residual products
!   numZ           ! maximal number of protons from initial compound nucleus
! Variables for input of ENDF structure
!   flagrecoil     ! flag to include recoil information
! Variables for ENDF limits, switches and tolerances
!   maxrp          ! maximum number of residual products
! Constants
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
!   xsepshigh      ! upper limit for cross sections in millibarns
!   xsepslow       ! lower limit for cross sections in millibarns
! Variables for info from TALYS
!   numinc         ! number of incident energies
! Variables for reaction initialization
!   Especindex     ! enegy index for spectra
!   Nenspec        ! number of incident energies for spectra
!   nlevmax        ! number of included discrete levels
!   numcut         ! number of energies before high - energy format
! Variables for initialization of ENDF format
!   idnum          ! number of different exclusive cross sections
!   MTid           ! channel identifier for MT - number
!   MTmax          ! highest MT number for exclusive channels
! Variables for total cross sections in ENDF format
!   xsnonel        ! nonelastic cross section
! Variables for partial cross sections in ENDF format
!   idchannel      ! identifier for channel
!   xsexcl         ! exclusive cross section
! Variables for residual production cross sections in ENDF format
!   Ehistcumrec    ! histogram emission energy for recoil spectra
!   Erecrp         ! recoil energy of residual nucleus
!   f0cumrec       ! energy distribution for recoil spectra
!   isorpexist     ! flag for existence of isomer of residual nuclide
!   nbegcumrec     ! first outgoing energy
!   nendcumrec     ! last outgoing energy
!   noutrecrp      ! number of recoil energies for residual nucleus
!   recrp          ! recoils or residual nucleus
!   rpexist        ! flag for existence of residual nuclide
!   xsrp           ! residual production cross section
!   xsrpiso        ! residual production cross section for isomer
!   Yrp            ! residual production yield
!   Yrpiso         ! residual production yield for isomer
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagex                              ! flag for existence
  logical   :: limit                               ! help variable
  integer   :: i                                   ! counter
  integer   :: id                                  ! counter for deuterons
  integer   :: idc                                 ! help variable
  integer   :: imax                                ! help variable
  integer   :: iyield                              ! particle yield
  integer   :: j                                   ! counter
  integer   :: MT                                  ! MT-number
  integer   :: nbeg0                               ! first outgoing energy
  integer   :: nend0                               ! last outgoing energy
  integer   :: nen                                 ! energy counter
  integer   :: nenout                              ! counter for outgoing energy
  integer   :: nex                                 ! discrete level
  integer   :: nin                                 ! counter for incident energy
  integer   :: Nix                                 ! neutron number index for residual nucleus
  integer   :: Nix1                                ! neutron index of residual nucleus
  integer   :: Nmax(numrp)                         ! maximum N of residual nucleus
  integer   :: Nrp                                 ! number of residual products
  integer   :: nrpchan(0:numZ, 0:numN)             ! number of MT numbers for residual product
  integer   :: Ntmp                                ! help variable
  integer   :: rpex(0:numZ, 0:numN, 20)            ! index for residual product
  integer   :: type                                ! particle type
  integer   :: Zix                                 ! charge number index for residual nucleus
  integer   :: Zix1                                ! proton index of residual nucleus
  integer   :: Zmax(numrp)                         ! maximum Z of residual nucleus
  integer   :: Ztmp                                ! help variable
  real(sgl) :: xsmax(numrp)                        ! maximum cross sections
  real(sgl) :: xsnonin                             ! sum of partial cross sections with MT number
  real(sgl) :: xsnonout                            ! sum of partial cross sections without MT number
  real(sgl) :: xsrpin(0:numZ, 0:numN)              ! sum of residual production cross sections with MT number
  real(sgl) :: xsrpout                             ! sum of residual production cross sections without MT  number
  real(sgl) :: xss                                 ! help variable
  real(sgl) :: xstmp                               ! help variable
  real(sgl) :: xsy                                 ! help variable
!
! ******* Determine limits for residual production cross sections ******
!
  Nrp = 0
  do Zix = 0, numZ
    do Nix = 0, numN
!
! Total
!
      rpexist(Zix, Nix) = .false.
      limit = .false.
      do nin = 1, numinc
        if (xsrp(Zix, Nix, nin) < xsepslow) xsrp(Zix, Nix, nin) = 0.
        if (nin >= numcut .and. xsrp(Zix, Nix, nin) >= xsepshigh) limit = .true.
      enddo
      if (limit) then
        rpexist(Zix, Nix) = .true.
        Nrp = Nrp + 1
      else
        do nin = 1, numinc
          xsrp(Zix, Nix, nin) = 0.
        enddo
      endif
      do nin = 1, numinc
        if (xsnonel(nin) > 0.) then
          Yrp(Zix, Nix, nin) = xsrp(Zix, Nix, nin) / real(xsnonel(nin))
        else
          Yrp(Zix, Nix, nin) = 0.
        endif
      enddo
!
! Ground state and isomer
!
      do nex = 0, nlevmax
        isorpexist(Zix, Nix, nex) = .false.
        limit = .false.
        do nin = 1, numinc
          if (xsrpiso(Zix, Nix, nex, nin) < xsepslow) xsrpiso(Zix, Nix, nex, nin) = 0.
          if (nin >= numcut .and. xsrpiso(Zix, Nix, nex, nin) >= xsepshigh) limit = .true.
        enddo
        if (limit) then
          isorpexist(Zix, Nix, nex) = .true.
        else
          do nin = 1, numinc
            xsrpiso(Zix, Nix, nex, nin) = 0.
          enddo
        endif
      enddo
      do nex = 0, nlevmax
        do nin = 1, numinc
          if (xsnonel(nin) > 0.) then
            Yrpiso(Zix, Nix, nex, nin) = xsrpiso(Zix, Nix, nex, nin) / real(xsnonel(nin))
          else
            Yrpiso(Zix, Nix, nex, nin) = 0.
          endif
        enddo
      enddo
    enddo
  enddo
!
! ************* Only include maxrp largest residual products ***********
!
  if (Nrp > maxrp) then
    imax = 0
    do Zix = 0, numZ
      do Nix = 0, numN
        if ( .not. rpexist(Zix, Nix)) cycle
        imax = imax + 1
        xsmax(imax) = 0.
        Zmax(imax) = Zix
        Nmax(imax) = Nix
        do nin = 1, numinc
          xsmax(imax) = max(xsrp(Zix, Nix, nin), xsmax(imax))
        enddo
      enddo
    enddo
    do i = 1, imax
      do j = 1, i
        if (xsmax(i) >= xsmax(j)) then
          xstmp = xsmax(i)
          Ztmp = Zmax(i)
          Ntmp = Nmax(i)
          xsmax(i) = xsmax(j)
          Zmax(i) = Zmax(j)
          Nmax(i) = Nmax(j)
          xsmax(j) = xstmp
          Zmax(j) = Ztmp
          Nmax(j) = Ntmp
        endif
      enddo
    enddo
    do i = maxrp + 1, imax
      rpexist(Zmax(i), Nmax(i)) = .false.
    enddo
  endif
!
! ************* Determine yields for (n,anything) channel **************
!
! This is a correction for the residual production yields for the case where specific MT numbers run to high energies
!
  do Zix = 0, numZ
    do Nix = 0, numN
      nrpchan(Zix, Nix) = 0
      do i = 1, 20
        rpex(Zix, Nix, i) = 0
      enddo
      do idc = -1, idnum
        id = idchannel(idc)
        do MT = 4, MTmax
          if (MT == 18 .or. (MT >= 51 .and. MT <= 91)) cycle
          if (MTid(MT) == id) then
            Zix1 = 0
            Nix1 = 0
            do type = 1, 6
              iyield = mod(id, 10 **(7 - type)) / (10 **(6 - type))
              Zix1 = Zix1 + iyield * parZ(type)
              Nix1 = Nix1 + iyield * parN(type)
            enddo
            if (Zix1 == Zix .and. Nix1 == Nix) then
              nrpchan(Zix, Nix) = nrpchan(Zix, Nix) + 1
              rpex(Zix, Nix, nrpchan(Zix, Nix)) = idc
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  do nin = 1, numcut - 1
    xsnonin = 0.
    do idc = -1, idnum
       id = idchannel(idc)
       do MT = 4, MTmax
         if (MT >= 51 .and. MT <= 91) cycle
         if (MTid(MT) == id) then
           xss = xsexcl(idc, nin)
           xsnonin = xsnonin + xss
           exit
         endif
       enddo
    enddo
    do Zix = 0, numZ
      do Nix = 0, numN
        xsrpin(Zix, Nix) = 0.
        do i = 1, nrpchan(Zix, Nix)
          idc = rpex(Zix, Nix, i)
          xsrpin(Zix, Nix) = xsrpin(Zix, Nix) + xsexcl(idc, nin)
        enddo
        xsrpout = max(xsrp(Zix, Nix, nin) - xsrpin(Zix, Nix), 0.)
        xsnonout = real(xsnonel(nin)) - xsnonin
        if (xsnonout > 0.) Yrp(Zix, Nix, nin) = xsrpout / xsnonout
      enddo
    enddo
  enddo
!
! Set residual production yields to zero for cases where all exclusive channels exist.
!
  do nin = 1, numcut - 1
    if (MTid(4) /=  -1) Yrp(0, 1, nin) = 0.
    if (MTid(16) /=  -1) Yrp(0, 2, nin) = 0.
    if (MTid(17) /=  -1) Yrp(0, 3, nin) = 0.
    if (MTid(37) /=  -1) Yrp(0, 4, nin) = 0.
    if (MTid(102) /=  -1) Yrp(0, 0, nin) = 0.
    if (MTid(103) /=  -1) Yrp(1, 0, nin) = 0.
    if (MTid(104) /=  -1 .and. MTid(28) /=  -1) Yrp(1, 1, nin) = 0.
    if (MTid(105) /=  -1 .and. MTid(32) /=  -1 .and. MTid(41) /=  -1) Yrp(1, 2, nin) = 0.
    if (MTid(106) /=  -1 .and. MTid(44) /=  -1 .and. MTid(115) /=  -1) Yrp(2, 1, nin) = 0.
  enddo
  if (numcut == numinc) then
    do Zix = 0, numZ
      do Nix = 0, numN
        flagex = .false.
        do nin = 1, numcut
          if (Yrp(Zix, Nix, nin) > 0.) flagex = .true.
        enddo
        if ( .not. flagex) rpexist(Zix, Nix) = .false.
      enddo
    enddo
  endif
!
! **************** Cumulative recoil spectra ***************************
!
  if (flagrecoil) then
    do Zix = 0, numZ
      do Nix = 0, numN
        if ( .not. rpexist(Zix, Nix)) cycle
        do nen = 1, Nenspec
          nin = Especindex(nen)
!
! Determine first and last outgoing energy
!
          nbeg0 = 1
          nend0 = noutrecrp(Zix, Nix, nen)
          do nenout = noutrecrp(Zix, Nix, nen), 2, -1
            if (recrp(Zix, Nix, nen, nenout) > 0..and. recrp(Zix, Nix, nen, nenout - 1) > 0.) then
              nend0 = nenout
              exit
            endif
          enddo
          nbegcumrec(Zix, Nix, nen) = nbeg0
          nendcumrec(Zix, Nix, nen) = nend0
!
! Histogram recoil energies
!
          Ehistcumrec(Zix, Nix, nen, 1) = 0.
          do nenout = 2, noutrecrp(Zix, Nix, nen)
            Ehistcumrec(Zix, Nix, nen, nenout) = 0.5 * (Erecrp(Zix, Nix, nen, nenout - 1) + Erecrp(Zix, Nix, nen, nenout))
          enddo
!
! Normalization of spectra
!
          xsy = 0.
          do nenout = nbeg0 - 1, nend0 - 1
            xsy = xsy + (Ehistcumrec(Zix, Nix, nen, nenout + 1) - Ehistcumrec(Zix, Nix, nen, nenout)) * recrp(Zix, Nix, nen, nenout)
          enddo
          if (xsy < xsepslow) then
            nbegcumrec(Zix, Nix, nen) = 0
            nendcumrec(Zix, Nix, nen) = 0
            cycle
          endif
          do nenout = nbeg0, nend0
            f0cumrec(Zix, Nix, nen, nenout) = max(recrp(Zix, Nix, nen, nenout) / xsy, 1.e-13)
          enddo
        enddo
      enddo
    enddo
  endif
  return
end subroutine processresidual
! Copyright A.J. Koning 2021
