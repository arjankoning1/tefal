subroutine processchannels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process exclusive channel data
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
!   nummt         ! number of MT numbers
! Variables for input of ENDF library type
!   flaggpf       ! flag for general purpose library
! Variables for input of ENDF structure
!   flagaddlow    ! flag to add low - energy Q>0 reactions to total cross section
!   flagngn       ! flag to include (n, gamma n) data
!   flagpara      ! flag to include partial cross sections for alphs
!   flagpard      ! flag to include partial cross sections for deuterons
!   flagparp      ! flag to include partial cross sections for protons
!   flagpart      ! flag to include partial cross sections for tritons
!   flagrecoil    ! flag to include recoil information
!   flagrenorm    ! flag for renormalization of spectra
! Constants
!   parA          ! mass number of particle
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
!   xsepshigh     ! upper limit for cross sections in millibarns
!   xsepslow      ! lower limit for cross sections in millibarns
! Variables for reaction initialization
!   Especindex    ! enegy index for spectra
!   massN         ! mass of nucleus in neutron units
!   Nenspec       ! number of incident energies for spectra
!   nlevmax       ! number of included discrete levels
!   numcut        ! number of energies before high - energy format
! Variables for info from TALYS
!   eninc         ! incident energy
! Variables for initialization of ENDF format
!   idnum         ! number of different exclusive cross sections
!   MTid          ! channel identifier for MT - number
! Variables for channel cross sections in ENDF format
!   branchiso     ! branching ratio for isomer
!   Ehist         ! histogram emission energy
!   Ehistrec      ! histogram recoil energy
!   Eout          ! emission energy
!   Eparticles    ! total energy carried away by particles
!   Erec          ! recoil energy
!   Erecav        ! average recoil energy
!   f0ex          ! energy distribution for exclusive channel
!   f0exrec       ! energy distribution for recoil
!   idchannel     ! identifier for channel
!   isoexist      ! flag for existence of isomer
!   nbeg          ! first outgoing energy
!   nbegrec       ! first outgoing energy
!   nend          ! last outgoing energy
!   nendrec       ! last outgoing energy
!   nout          ! number of emission energies
!   noutrec       ! number of recoil energies
!   Qexcl         ! Q - value
!   Qexcliso      ! Q - value for isomer
!   recexcl       ! exclusive recoils
!   specexcl      ! exclusive spectra
!   xsexcl        ! exclusive cross section
!   xsexcliso     ! exclusive cross section for isomer
!   xsnonth       ! sum of all non - thr. reactions except (n, g) and (n, f)
! Variables for discrete state cross sections in ENDF format
!   Qdisc         ! Q - value
!   xscont        ! continuum cross section
!   xsngn         ! (projectile, gamma - ejectile) cross section
! Variables for spectra in ENDF format
!   Eocum         ! emission energies for total production spectra
!   ncumout       ! number of emission energies for total production spectra
!   xsemis        ! total production emission spectra
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagadd                  ! flag for addition of discrete states to spectra
  logical   :: limit                    ! help variable
  integer   :: id                       ! counter for deuterons
  integer   :: idc                      ! help variable
  integer   :: iyield                   ! particle yield
  integer   :: MT                       ! MT-number
  integer   :: N                        ! neutron number of residual nucleus
  integer   :: nbeg0                    ! first outgoing energy
  integer   :: nen                      ! energy counter
  integer   :: nen2                     ! energy counter
  integer   :: nend0                    ! last outgoing energy
  integer   :: nenout                   ! counter for outgoing energy
  integer   :: nex                      ! discrete level
  integer   :: nin                      ! counter for incident energy
  integer   :: nin2                     ! energy counter
  integer   :: Nix                      ! neutron number index for residual nucleus
  integer   :: type                     ! particle type
  integer   :: Zix                      ! charge number index for residual nucleus
  real(sgl) :: E0                       ! constant of temperature formula
  real(sgl) :: ea                       ! help variable
  real(sgl) :: Eaverage                 ! average emission energy
  real(sgl) :: eb                       ! help variable
  real(sgl) :: ee                       ! energy
  real(sgl) :: Elast                    ! help variable
  real(sgl) :: Eo(0:numen2)             ! help variable
  real(dbl) :: QQ                       ! Q-value
  real(sgl) :: specmass                 ! specific mass
  real(sgl) :: tot                      ! total
  real(sgl) :: totE                     ! total
  real(sgl) :: xsex                     ! help variable
  real(sgl) :: xsint                    ! interpolated cross section
  real(sgl) :: xso(0:numen2)            ! help variable
  real(sgl) :: xsy                      ! help variable
  real(sgl) :: y1                       ! coordinates of the 1st summit of the triangle
  real(sgl) :: y2                       ! coordinates of the 2nd summit of the triangle
!
! ******* Determine limits for exclusive channel cross sections ********
!
  do idc = 0, idnum
!
! Ensure that capture cross sections are always greater than zero.
!
    if (idchannel(idc) == 0) then
      do nin = numcut - 1, 1, -1
        if (xsexcl(idc, nin) == 0.) xsexcl(idc, nin) = xsexcl(idc, nin + 1)
      enddo
    endif
    limit = .false.
    do nin = 1, numcut
      if (xsexcl(idc, nin) < xsepslow) then
        if (Qexcl(idc) > 0.) then
          xsexcl(idc, nin) = xsepslow
        else
          xsexcl(idc, nin) = 0.
        endif
      endif
      if (xsexcl(idc, nin) >= xsepshigh) limit = .true.
    enddo
    if ( .not. limit) then
      do nin = 1, numcut
        xsexcl(idc, nin) = 0.
      enddo
    endif
    do nex = 0, nlevmax
      if ( .not. isoexist(idc, nex)) cycle
      limit = .false.
      do nin = 1, numcut
        if (xsexcliso(idc, nex, nin) < xsepslow) then
          if (Qexcliso(idc, nex) > 0.) then
            xsexcliso(idc, nex, nin) = xsepslow
          else
            xsexcliso(idc, nex, nin) = 0.
          endif
        endif
        if (xsexcliso(idc, nex, nin) /= 0..and. xsexcl(idc, nin) >= xsepshigh) limit = .true.
      enddo
      if ( .not. limit) then
        do nin = 1, numcut
          xsexcliso(idc, nex, nin) = 0.
        enddo
      endif
    enddo
!
! Sum all non-threshold reaction with the exception of capture and fission to be added to the nonelastic and total cross
! section as background in the resonance range.
!
    id = idchannel(idc)
    if (flagaddlow .and. Qexcl(idc) > 0 .and. id /= 0 .and. id /=  - 99) then
      do MT = 1, nummt
        if (id == MTid(MT) .and. MTid(MT) >= 0) then
          do nin = 1, numcut
            xsnonth(nin) = xsnonth(nin) + xsexcl(idc, nin)
          enddo
          exit
        endif
      enddo
    endif
  enddo
!
! ************************* Energy spectra *****************************
!
  if (flaggpf) then
!
! For lowest incident energies, relate the exclusive capture and (n,gx) spectra to those of the first energy calculated by TALYS.
!
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (nout(0, nen) /= 0) then
        do nen2 = 1, nen - 1
          nin2 = Especindex(nen2)
          nout(0, nen2) = nout(0, nen)
          do nenout = 1, nout(0, nen)
            do type = 0, 6
              Ehist(0, nen, type, nenout) = 0.5 * (Eout(0, nen, nenout - 1) + Eout(0, nen, nenout))
              Ehist(0, nen2, type, nenout) = Ehist(0, nen, type, nenout)
              specexcl(0, nen2, type, nenout) = specexcl(0, nen, type, nenout) * xsexcl(0, nin2) / xsexcl(0, nin)
            enddo
          enddo
        enddo
        exit
      endif
    enddo
!
! Exclusive particle spectra
!
    do idc = 0, idnum
      id = idchannel(idc)
      Zix = 0
      Nix = 0
      do type = 1, 6
        iyield = mod(id, 10 **(7 - type)) / (10 **(6 - type))
        if (iyield == 0) cycle
        Zix = Zix + iyield * parZ(type)
        Nix = Nix + iyield * parN(type)
      enddo
      do nen = 1, Nenspec
        if (nout(idc, nen) == 0) cycle
        nin = Especindex(nen)
        xsex = xsexcl(idc, nin)
!
! Add high-energy discrete peaks to (n,p)...(n,alpha) spectra,
! if required.
!
        if (id == 100000 .or. id == 10000 .or. id == 1000 .or. id == 100 .or. id == 10 .or. id == 1) then
          flagadd = .false.
          if (id == 100000) then
            type = 1
            if (flagparn) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (id == 10000) then
            type = 2
            if (flagparp) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (id == 1000) then
            type = 3
            if (flagpard) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (id == 100) then
            type = 4
            if (flagpart) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (id == 10) then
            type = 5
            if (flagpart) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (id == 1) then
            type = 6
            if (flagpara) then
              xsex = xscont(type, nin)
            else
              flagadd = .true.
            endif
          endif
          if (flagadd) then
            E0 = eninc(nin) + Qexcl(idc)
            Elast = max(E0 - 7., Eocum(type, nen, 1))
            Eo(0) = 0.
            do nen2 = 1, ncumout(type, nen)
              Eo(nen2) = Eocum(type, nen, nen2)
              xso(nen2) = xsemis(type, nen, nen2)
            enddo
            do nenout = 1, nout(idc, nen)
              ee = Eout(idc, nen, nenout)
              if (ee > Elast .and. ee <= E0 .and. ee < Eo(ncumout(type, nen))) then
                call locate(Eo, 1, ncumout(type, nen), ee, nen2)
                ea = Eo(nen2)
                eb = Eo(nen2 + 1)
                y1 = xso(nen2)
                y2 = xso(nen2 + 1)
                call pol1(ea, eb, y1, y2, ee, xsint)
                specexcl(idc, nen, type, nenout) = xsint
              endif
            enddo
          endif
        endif
!
! Determine first and last outgoing energy
!
        if (xsex == 0.) cycle
        do type = 1, 6
          nbeg(idc, type, nen) = 0
          nend(idc, type, nen) = 0
          iyield = mod(id, 10 **(7 - type)) / (10 **(6 - type))
          if (id /= 0 .and. iyield == 0) cycle
          nbeg0 = 1
          nend0 = nout(idc, nen)
          do nenout = 2, nout(idc, nen)
            if (specexcl(idc, nen, type, nenout) > 0.) then
              nbeg0 = nenout - 1
              exit
            endif
          enddo
          do nenout = nout(idc, nen), 1, -1
            if (specexcl(idc, nen, type, nenout) > 0.) then
              nend0 = nenout
              exit
            endif
          enddo
!
! Histogram emission energies
!
          specmass = (massN - Zix - Nix + parA(type)) / (massN - Zix - Nix + 2 * parA(type))
          QQ = Qdisc(type, 0)
          Elast = specmass * (eninc(nin) * specmass + QQ) * (1. - 1.e-5)
          Ehist(idc, nen, type, 1) = 0.
          N = nout(idc, nen)
          do nenout = 2, N
            Ehist(idc, nen, type, nenout) = 0.5 * (Eout(idc, nen, nenout - 1) + Eout(idc, nen, nenout)) * specmass
          enddo
          nbeg(idc, type, nen) = nbeg0
          nend(idc, type, nen) = nend0
!
! Normalization of spectra
!
          nend0 = nend0 + 1
          nend(idc, type, nen) = nend0
          if (nend0 > 2) then
            Ehist(idc, nen, type, nend0) = 2. * Ehist(idc, nen, type, nend0 - 1) - Ehist(idc, nen, type, nend0 - 2)
          else
            Ehist(idc, nen, type, nend0) = Ehist(idc, nen, type, nend0 - 1)
          endif
          f0ex(idc, nen, type, nend0) = 0.
          xsy = 0.
          do nenout = nbeg0 - 1, nend0 - 1
            xsy = xsy + (Ehist(idc, nen, type, nenout + 1) - Ehist(idc, nen, type, nenout)) * specexcl(idc, nen, type, nenout)
          enddo
          if (xsy == 0.) then
            if (id == 0) then
              xsngn(type, nin) = 0.
            else
              if (id == 100000) xscont(1, nin) = 0.
              if (id == 10000) xscont(2, nin) = 0.
              if (id == 1000) xscont(3, nin) = 0.
              if (id == 100) xscont(4, nin) = 0.
              if (id == 10) xscont(5, nin) = 0.
              if (id == 1) xscont(6, nin) = 0.
            endif
            nbeg(idc, type, nen) = 0
            nend(idc, type, nen) = 0
            cycle
          endif
          if ( .not. flagrenorm) then
            if (id == 0) then
              if (flagngn) xsy = xsngn(type, nin)
            else
              xsy = xsex * iyield
            endif
            if (xsy == 0.) then
              nbeg(idc, type, nen) = 0
              nend(idc, type, nen) = 0
              cycle
            endif
          endif
          do nenout = nbeg0, nend0
            f0ex(idc, nen, type, nenout) = max(specexcl(idc, nen, type, nenout) / xsy, 1.e-13)
          enddo
          f0ex(idc, nen, type, nend0) = 0.
!
! Determine average energy of particles for energy balance
! normalization
!
          if (flagrenorm) then
            tot = 0.
            totE = 0.
            do nenout = nbeg0 - 1, nend0 - 1
              tot = tot + (Ehist(idc, nen, type, nenout + 1) - Ehist(idc, nen, type, nenout)) * f0ex(idc, nen, type, nenout)
              totE = totE + 0.5 * (Ehist(idc, nen, type, nenout) + Ehist(idc, nen, type, nenout + 1)) * &
                (Ehist(idc, nen, type, nenout + 1) - Ehist(idc, nen, type, nenout)) * f0ex(idc, nen, type, nenout)
            enddo
            if (tot > 0.) then
              Eaverage = totE / tot
              Eparticles(idc, nen) = Eparticles(idc, nen) + iyield * Eaverage
            endif
          endif
        enddo
        if (idc /= 0 .or. Ehist(idc, nen, 0, 1) == 0.) then
          do nenout = 1, nout(idc, nen)
            Ehist(idc, nen, 0, nenout) = 0.5 * (Eout(idc, nen, nenout - 1) + Eout(idc, nen, nenout))
          enddo
        endif
!
! Recoils
!
        if (flagrecoil) then
          if (noutrec(idc, nen) == 0) cycle
          nbeg0 = 1
          nend0 = noutrec(idc, nen)
          do nenout = noutrec(idc, nen), 1, -1
            if (recexcl(idc, nen, nenout) > 0.) then
              nend0 = nenout
              exit
            endif
          enddo
          nbegrec(idc, nen) = nbeg0
          nendrec(idc, nen) = nend0
          Ehistrec(idc, nen, 1) = 0.
          do nenout = 2, noutrec(idc, nen)
            Ehistrec(idc, nen, nenout) = 0.5 * (Erec(idc, nen, nenout - 1) + Erec(idc, nen, nenout))
          enddo
          tot = 0.
          totE = 0.
          do nenout = nbeg0 - 1, nend0 - 1
            tot = tot + (Ehistrec(idc, nen, nenout + 1) - Ehistrec(idc, nen, nenout)) * recexcl(idc, nen, nenout)
            totE = totE + 0.5 * (Ehistrec(idc, nen, nenout) + Ehistrec(idc, nen, nenout + 1)) * &
              (Ehistrec(idc, nen, nenout + 1) - Ehistrec(idc, nen, nenout)) * recexcl(idc, nen, nenout)
          enddo
          if (tot > 1.e-30) then
            Erecav(idc, nen) = totE / tot
          else
            nbegrec(idc, nen) = 0
            nendrec(idc, nen) = 0
            cycle
          endif
          do nenout = nbeg0, nend0
            f0exrec(idc, nen, nenout) = max(recexcl(idc, nen, nenout) / tot, 1.e-13)
          enddo
        endif
      enddo
    enddo
  endif
!
! ************ Branching ratios for ground state and isomers ***********
!
  do idc = 0, idnum
    do nex = 0, nlevmax
      if ( .not. isoexist(idc, nex)) cycle
      do nin = numcut, 1, -1
        if (xsexcl(idc, nin) > 0.) then
          branchiso(idc, nex, nin) = xsexcliso(idc, nex, nin) / xsexcl(idc, nin)
          if (nin < numcut .and. xsexcliso(idc, nex, nin) < 1.01*xsepslow) branchiso(idc, nex, nin) = branchiso(idc, nex, nin + 1)
        else
          branchiso(idc, nex, nin) = 0.
        endif
      enddo
    enddo
  enddo
  return
end subroutine processchannels
! Copyright A.J. Koning 2021
