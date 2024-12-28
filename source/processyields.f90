subroutine processyields
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process yields and spectra
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
! Variables for input of ENDF structure
!   flagngn       ! flag to include (n, gamma n) data
!   flagrenorm    ! flag for renormalization of spectra
!   flagtabddx    ! flag to give explicit DDX in MF6
! Variables for reacion initialization
!   Especindex    ! energy index for spectra
!   massN         ! mass of nucleus in neutron units
!   Nenspec       ! number of incident energies for spectra
!   numcut        ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc         ! incident energy
!   k0            ! index of incident particle
!   numinc        ! number of incident energies
! Constants
!   parA          ! mass number of particle
!   xsepshigh     ! upper limit for cross sections in millibarns
!   xsepslow      ! lower limit for cross sections in millibarns
! Variables for initialization of ENDF format
!   idnum         ! number of different exclusive cross sections
!   MTid          ! channel identifier for MT - number
!   MTmax         ! highest MT number for exclusive channels
! Variables for total cross sections in ENDF format
!   xsnonel       ! nonelastic cross section
! Variables for partial cross sections in ENDF format
!   Eout          ! emission energy
!   idchannel     ! identifier for channel
!   nout          ! number of emission energies
!   specexcl      ! exclusive spectra
!   xsexcl        ! exclusive cross section
! Variables for discrete state cross sections in ENDF format
!   Qdisc         ! Q - value
!   xsngn         ! (projectile, gamma - ejectile) cross section
! Variables for spectra in ENDF format
!   buratio       ! break - up ratio
!   ddxemis       ! double - differential emission spectra
!   Eocum         ! emission energies for total production spectra
!   ncumout       ! number of emission energies for total production spectra
!   Nddx          ! number of angles
!   preeqratio    ! pre - equilibrium ratio
!   xsemis        ! total production emission spectra
! Variables for yields in ENDF format
!   Ehistcum      ! histogram emission energy for total production spectra
!   f0cum         ! energy distribution for cumulative spectra
!   f0ddx         ! energy distribution for DDX
!   nbegcum       ! first outgoing energy
!   nendcum       ! last outgoing energy
!   xsprod        ! particle production cross section
!   yieldany      ! yield for (n, anything) channel
!   yieldp        ! particle production yield
!
! *** Declaration of local data
!
  implicit none
  logical   :: limit                ! help variable
  integer   :: i                    ! counter
  integer   :: id                   ! counter for deuterons
  integer   :: idc                  ! help variable
  integer   :: iyield               ! particle yield
  integer   :: j                    ! counter
  integer   :: k                    ! counter
  integer   :: l                    ! counter
  integer   :: MT                   ! MT-number
  integer   :: MT0                  ! MT number
  integer   :: N                    ! neutron number of residual nucleus
  integer   :: nbeg0                ! first outgoing energy
  integer   :: nen                  ! energy counter
  integer   :: nen2                 ! energy counter
  integer   :: nend0                ! last outgoing energy
  integer   :: nenout               ! counter for outgoing energy
  integer   :: nin                  ! counter for incident energy
  integer   :: type                 ! particle type
  real(sgl) :: denom                ! help variable
  real(sgl) :: eb                   ! help variable
  real(sgl) :: ee                   ! energy
  real(sgl) :: Ein                  ! incident energy
  real(sgl) :: Elast                ! help variable
  real(sgl) :: enum                 ! enumerator of Lorentzian
  real(dbl) :: QQ                   ! Q-value
  real(sgl) :: specmass             ! specific mass
  real(sgl) :: xb                   ! begin x value
  real(sgl) :: xe                   ! end x value
  real(sgl) :: xsnonin              ! sum of partial cross sections with MT number
  real(sgl) :: xsnonout             ! sum of partial cross sections without MT number
  real(sgl) :: xsprodin             ! sum of production cross sections with MT number
  real(sgl) :: xsprodout            ! sum of production cross sections without MT number
  real(sgl) :: xss                  ! help variable
  real(sgl) :: xsy                  ! help variable
!
! ***************** Determine limits for particle yields ***************
!
  do type = 0, 6
    limit = .false.
    do nin = 1, numinc
      if (xsprod(type, nin) < xsepslow) then
        xsprod(type, nin) = 0.
        yieldp(type, nin) = 0.
      endif
      if (xsprod(type, nin) >= xsepshigh) limit = .true.
    enddo
    if ( .not. limit) then
      do nin = 1, numinc
        xsprod(type, nin) = 0.
        yieldp(type, nin) = 0.
      enddo
    endif
  enddo
!
! ******* Determine limits for particle spectra and normalization ******
!
! Cumulative particle and photon spectra
!
  do type = 0, 6
    do nen = 1, Nenspec
      nin = Especindex(nen)
      nbegcum(type, nen) = 0
      nendcum(type, nen) = 0
      if (xsprod(type, nin) == 0.) cycle
!
! Determine first and last outgoing energy
!
      nbeg0 = 1
      nend0 = ncumout(type, nen)
      do nenout = 2, ncumout(type, nen)
        if (xsemis(type, nen, nenout) > 0.) then
          nbeg0 = nenout - 1
          exit
        endif
      enddo
      do nenout = ncumout(type, nen), 1, -1
        if (xsemis(type, nen, nenout) > 0.) then
          nend0 = nenout
          exit
        endif
      enddo
      if (type == k0) then
        Ein = eninc(nin)
        k = nend0
        do nenout = k, 1, -1
          if (Ein > Eocum(type, nen, nenout)) then
            nend0 = nenout
            exit
          endif
        enddo
      endif
!
! Histogram emission energies
! Corrections applied thanks to observations of Alexander konobeyev and JC Sublet.
!
      Ehistcum(type, nen, 1) = 0.
      specmass = massN / (massN + parA(type))
      QQ = Qdisc(type, 0)
      Elast = specmass * (eninc(nin) * specmass + QQ) * (1. - 1.e-5)
      N = ncumout(type, nen)
      do nenout = 2, N
        Ehistcum(type, nen, nenout) = 0.5 * (Eocum(type, nen, nenout - 1) + Eocum(type, nen, nenout)) * specmass
        if (Elast > 0. .and. Ehistcum(type, nen, nenout) > Elast) then
          Ehistcum(type, nen, nenout) = Elast
          ncumout(type, nen) = nenout - 1
          nend0 = nenout - 1
          exit
        endif
      enddo
      nbegcum(type, nen) = nbeg0
      nendcum(type, nen) = nend0
!
! Pre-equilibrium ratios
!
      do nenout = 2, ncumout(type, nen) + 1
        if (preeqratio(type, nen, nenout) == 0..and. preeqratio(type, nen, nenout - 1) > 0.) &
 &        preeqratio(type, nen, nenout) = preeqratio(type, nen, nenout - 1)
        if (buratio(type, nen, nenout) == 0..and. buratio(type, nen, nenout - 1) > 0.) &
 &        buratio(type, nen, nenout) = buratio(type, nen, nenout - 1)
      enddo
!
! Normalization of spectra
!
      nend0 = nend0 + 1
      nendcum(type, nen) = nend0
      if (nend0 > 2 .and. Elast > 0.) then
        Ehistcum(type, nen, nend0) = Elast
      else
        Ehistcum(type, nen, nend0) = Ehistcum(type, nen, nend0 - 1)
      endif
      f0cum(type, nen, nend0) = 0.
      xsy = 0.
      do nenout = nbeg0 - 1, nend0 - 1
        xsy = xsy + (Ehistcum(type, nen, nenout + 1) - Ehistcum(type, nen, nenout)) * xsemis(type, nen, nenout)
      enddo
      if (xsy <= 1.e-30) then
        xsprod(type, nen) = 0.
        nbegcum(type, nen) = 0
        nendcum(type, nen) = 0
        cycle
      endif
      if ( .not. flagrenorm) xsy = xsprod(type, nin)
      do nenout = nbeg0, nend0
        f0cum(type, nen, nenout) = max(xsemis(type, nen, nenout) / xsy, 1.e-13)
      enddo
      f0cum(type, nen, nend0) = 0.
    enddo
  enddo
!
! ******* Determine yields and spectra for (n,anything) channel ********
!
! This is a correction for the total yields and spectra for the case where specific MT numbers run to high energies
!
  do type = 0, 6
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (k0 == 1 .and. type == 1) then
        yieldany(type, nin) = 1.
      else
        yieldany(type, nin) = 0.
      endif
      if (nin >= numcut) cycle
      xsnonin = 0.
      xsprodin = 0.
      nbeg0 = nbegcum(type, nen)
      nend0 = nendcum(type, nen)
      do idc = -1, idnum
        id = idchannel(idc)
        iyield = 0
        do MT = 4, MTmax
          if (MT >= 51 .and. MT <= 91) cycle
          if (MTid(MT) == id) then
            xss = xsexcl(idc, nin)
            MT0 = MT
            iyield = mod(id, 10 **(7 - type)) / (10 **(6 - type))
            xsnonin = xsnonin + xss
            if (iyield > 0 .and. MT0 /= 18) then
              xsprodin = xsprodin + iyield * xss
    Loop1:    do nenout = nbeg0, nend0
                do nen2 = 1, nout(idc, nen)
                  if (abs(Eocum(type, nen, nenout) - Eout(idc, nen, nen2)) <= 0.001) then
                    xsemis(type, nen, nenout) = xsemis(type, nen, nenout) - specexcl(idc, nen, type, nen2)
                    cycle Loop1
                  endif
                enddo
              enddo Loop1
            endif
            exit
          endif
        enddo
      enddo
      xsprodout = xsprod(type, nin) - xsprodin
      xsnonout = real(xsnonel(nin)) - xsnonin
      if (flagngn) then
        enum = xsprodout + xsngn(type, nin)
        denom = xsnonout + xsngn(type, nin)
      else
        enum = xsprodout
        denom = xsnonout
      endif
      if (enum > 0.1) then
        if (denom > 0.1) yieldany(type, nin) = enum / denom
        nend0 = min(numen2, nend0 + 1)
        f0cum(type, nen, nend0) = 0.
        xsy = 0.
        if (nbeg0 >= 1) then
          do nenout = nbeg0 - 1, nend0 - 1
            xsy = xsy + (Ehistcum(type, nen, nenout + 1) - Ehistcum(type, nen, nenout)) * xsemis(type, nen, nenout)
          enddo
        endif
        if (xsy == 0.) xsy = xsprodout
        if (xsy > 0.1) then
          do nenout = nbeg0, nend0
            f0cum(type, nen, nenout) = max(xsemis(type, nen, nenout) / xsy, 1.e-13)
          enddo
          f0cum(type, nen, nend0) = 0.
        endif
      endif
    enddo
  enddo
  do type = 0, 6
    do nin = 2, numcut - 1
      if (yieldany(type, nin) == 0.) then
        eb = 0.
        ee = 0.
        xb = 0.
        xe = 0.
        do j = 1, 3
          k = max(nin - j, 1)
          l = min(nin + j, numcut)
          if (xb == 0..and.yieldany(type, k) /= 0.) then
            eb = eninc(k)
            xb = yieldany(type, k)
          endif
          if (xe == 0..and.yieldany(type, l) /= 0.) then
            ee = eninc(l)
            xe = yieldany(type, l)
          endif
          if (xb /= 0. .and. xe /= 0. .and. ee > eb) then
            yieldany(type, nin) = xb + (eninc(nin) - eb) / (ee - eb) * (xe - xb)
            exit
          endif
        enddo
      endif
    enddo
  enddo
!
! **** For incident deuterons: Double-differential nucleon spectra *****
!
  if (flagtabddx) then
    do type = 1, 2
      do nen = 1, Nenspec
        do nenout = nbegcum(type, nen), nendcum(type, nen)
          if (f0cum(type, nen, nenout) > 0.) then
            do i = 1, Nddx
              f0ddx(type, nen, i, nenout) = ddxemis(type, nen, i, nenout) / f0cum(type, nen, nenout) * 0.5
            enddo
          endif
        enddo
      enddo
    enddo
  endif
  return
end subroutine processyields
! Copyright A.J. Koning 2021
