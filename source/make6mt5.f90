subroutine make6mt5(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF6 for MT5
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
!   numenin        ! number of incident energies
!   numiso         ! number of isomers
!   numN           ! maximal number of neutrons from initial compound nucleus
!   numZ           ! maximal number of protons from initial compound nucleus
! Variables for ENDF limits, switches and tolerances
!   cuteps         ! energy shift at MT5 cutoff energy (in eV)
! Variables for input of ENDF library type
!   flagbreakup    ! breakup flag
!   flagclean      ! flag to clean up double points
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for input of ENDF structure
!   flagrecoil     ! flag to include recoil information
!   flagrp10       ! flag to put residual production cross sections in MF10
!   flagrp6        ! flag to put residual production cross sections in MF6
!   flagtabddx     ! flag to give explicit DDX in MF6
! Constants
!   parA           ! mass number of particle
!   parmass        ! mass of particle in a.m.u.
!   parZ           ! charge number of particle
! Variables for reaction initialization
!   eninccut       ! last incident energy before high - energy format
!   Especindex     ! enegy index for spectra
!   Nenspec        ! number of incident energies for spectra
!   nlevmax        ! number of included discrete levels
!   numcut         ! number of energies before high - energy format
!   relmass        ! mass relative to neutron mass
! Variables for info from TALYS
!   Ainit          ! mass number of initial compound nucleus
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   numinc         ! number of incident energies
!   Zinit          ! charge number of initial compound nucleus
! Variables for initialization of ENDF format
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
! Variables for total cross sections in ENDF format
!   xsany          ! (x, anything) cross section (MT5)
!   xsnonel        ! nonelastic cross section
! Variables for residual production cross sections in ENDF format
!   nucmass        ! mass of nucleus
! Variables for partial cross sections in ENDF format
!   nout           ! number of emission energies
! Variables for spectra in ENDF format
!   buratio        ! break - up ratio
!   Nddx           ! number of angles
!   preeqratio     ! pre - equilibrium ratio
!   rmuddx         ! cosine of angle
! Variables for yields in ENDF format
!   Ehistcum       ! histogram emission energy for total production spectra
!   f0cum          ! energy distribution for cumulative spectra
!   f0ddx          ! energy distribution for DDX
!   nbegcum        ! first outgoing energy
!   nendcum        ! last outgoing energy
!   yieldany       ! yield for (n, anything) channel
!   yieldp         ! particle production yield
! Variables for residual production cross sections in ENDF format
!   Ehistcumrec    ! histogram emission energy for recoil spectra
!   Erpiso         ! energy of isomer
!   f0cumrec       ! energy distribution for recoil spectra
!   isorpexist     ! flag for existence of isomer of residual nuclide
!   nbegcumrec     ! first outgoing energy
!   nendcumrec     ! last outgoing energy
!   Nisorp         ! number of isomers
!   rpexist        ! flag for existence of residual nuclide
!   Yrp            ! residual production yield
!   Yrpiso         ! residual production yield for isomer
! Variables for ENDF format
!   NK             ! number of subsections
! Variables for MF3
!   E3             ! incident energy for MF3 (in ENDF - 6 format)
!   xs             ! cross section
! Variables for MF4
!   LCT            ! LAB / CM flag
! Variables for MF8
!   ELFS           ! excitation energy of final state
!   LFS            ! final state number
!   LMF            ! file number for information
!   MATP           ! material number for reaction product
!   ZAPr           ! designation of final nucleus
! Variables for MF6
!   AWP            ! product mass
!   b6             ! energy - angle values
!   b6gam          ! energy - angle values for photons
!   b6rec          ! energy - angle values for recoils
!   E6             ! incident energy (in ENDF - 6 format) for distribution
!   Ey             ! incident energy for yields (in ENDF - 6 format)
!   flagrec        ! flag to state that b6 element concerns recoil grid
!   INTER6ea       ! interpolation scheme
!   INTER6y        ! interpolation scheme
!   kpart          ! section number for particles
!   LANG           ! flag for angular representation
!   LAW            ! flag for distribution function
!   LEP            ! interpolation scheme for secondary energy
!   LIP            ! product modifier flag
!   NA             ! number of angular parameters
!   NBT6ea         ! separation value for interpolation scheme
!   NBT6y          ! separation value for interpolation scheme
!   ND             ! number of discrete energies
!   NE6ea          ! number of incident energies for distribution
!   NEP            ! number of secondary energy points
!   NP6y           ! number of incident energies for yields
!   NR6ea          ! number of interpolation ranges for distribution
!   NR6y           ! number of interpolation ranges for yields
!   NW             ! number of words
!   Y              ! product yield (in ENDF - 6 format)
!   ZAP            ! product identifier
!
! *** Declaration of local data
!
  implicit none
  integer   :: A                  ! mass number of target nucleus
  integer   :: i                  ! counter
  integer   :: ib                 ! counter
  integer   :: idc                ! help variable
  integer   :: iE                 ! energy counter
  integer   :: iNEP               ! counter for energy
  integer   :: ipos               ! help variable
  integer   :: iso                ! counter for isomer
  integer   :: j                  ! counter
  integer   :: k                  ! counter
  integer   :: kbeg               ! index to mark begin of word
  integer   :: kk                 ! counter
  integer   :: MT                 ! MT-number
  integer   :: nen                ! energy counter
  integer   :: nenout             ! counter for outgoing energy
  integer   :: nex                ! discrete level
  integer   :: nin                ! counter for incident energy
  integer   :: Nix                ! neutron number index for residual nucleus
  integer   :: type               ! particle type
  integer   :: Z                  ! charge number of target nucleus
  integer   :: Zix                ! charge number index for residual nucleus
  real(sgl) :: Eev                ! energy in eV
  real(sgl) :: Ein                ! incident energy
  real(sgl) :: enencut            ! incident energy at cutoff
  real(sgl) :: Erpcut             ! incident energy at cutoff
  real(sgl) :: xsy                ! help variable
  real(sgl) :: yy                 ! help variable
!
! ************ Inclusive yields and spectra for high energies **********
!
! 1. Product yields
!
  k = 0
  if ( .not. mtexist(3, MT)) return
  if (flaggpf) then
!
! A. Particles
!
    idc = 0
    do type = 1, 6
      k = k + 1
!
! 1. Product yields
!
      ZAP(k) = 1000. * parZ(type) + parA(type)
      AWP(k) = real(relmass(type))
      Ey(k, 1) = E3(MT, 1)
      if (flagendfdet .and. k0 == 1 .and. type == 1 .and. xs(MT, 1) > 0.) then
        Y(k, 1) = 1.
      else
        Y(k, 1) = 0.
      endif
      iE = 1
!
! Low energy (n,gamma z) yields
!
      if (flagendfdet) then
        do nin = 1, numcut
          if (xsany(nin) == 0..and. (nin < numcut .and. xsany(min(nin + 1, numenin)) == 0.)) &
            cycle
          Eev = eninc(nin) * 1.e6
          if (Eev <= Ey(k, 1)) cycle
          iE = iE + 1
          Ey(k, iE) = Eev
          Y(k, iE) = yieldany(type, nin)
          if (xs(MT, iE) == 0.) Y(k, iE) = 0.
        enddo
        mtexist(6, MT) = .true.
        mfexist(6) = .true.
      endif
!
! High energy yields
!
      if (flaghigh) then
        if (flagendfdet .and. numcut > 1) then
          iE = iE + 1
          Ey(k, iE) = eninccut * 1.e6 + cuteps
          Y(k, iE) = yieldp(type, numcut)
        endif
        do nin = numcut + 1, numinc
          if (nin < numinc) then
            if (yieldp(type, nin) == 0..and. yieldp(type, nin + 1) == 0.) cycle
          endif
          yy = yieldp(type, nin)
          iE = iE + 1
          Ey(k, iE) = eninc(nin) * 1.e6
          Y(k, iE) = yy
          mtexist(6, MT) = .true.
          mfexist(6) = .true.
        enddo
      endif
      if (iE <= 2) then
        k = k - 1
        cycle
      endif
      ipos = 0
      do i = 1, iE
        if (Y(k, i) > 0.) ipos = ipos + 1
      enddo
      if (ipos == 0) then
        k = k - 1
        cycle
      endif
      NP6y(k) = iE
      NR6y(k) = 1
      NBT6y(k, 1) = NP6y(k)
      INTER6y(k, 1) = 2
      LAW(k) = 1
      LIP(k) = 0
      if (flagbreakup .and. (type == 1 .or. type == 2)) then
        LANG(k) = 3
      else
        LANG(k) = 2
      endif
      LEP(k) = 1
!
! 2. Energy-angle distributions
!
! First two distributions are zero
!
      E6(k, 1) = Ey(k, 1)
      ND(k, 1) = 0
      if (flagbreakup .and. (type == 1 .or. type == 2)) then
        NA(k, 1) = 2
      else
        NA(k, 1) = 1
      endif
      NEP(k, 1) = 2
      b6(k, 1, 1) = 0.
      b6(k, 1, 2) = 1.
      if (flagtabddx .and. (type == 1 .or. type == 2) .and. Nddx > 0) then
        b6(k, 1, 3) = 1.
        b6(k, 1, 4) = 0.
        NA(k, 1) = 0
      else
        b6(k, 1, 3) = 0.
        b6(k, 1, 4) = 1.
      endif
      b6(k, 1, 5) = 0.
      b6(k, 1, 6) = 0.
      NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
!
! Loop over incident energies
!
! To avoid very big files, we do not give secondary information at each incident energy, but instead skip various energy points.
!
      iE = 1
      if (flagendfdet) then
        enencut = eninc(Especindex(Nenspec))
        do nen = 1, Nenspec
          nin = Especindex(nen)
          if (nin > numcut) then
            enencut = Ein
            exit
          endif
          Ein = eninc(nin)
          if (yieldany(type, nin) == 0..and. (nin < numcut .and. yieldany(type, min(nin + 1, numenin)) == 0.)) cycle
          if (xsany(nin) == 0..and. (nin < numcut .and. xsany(min(nin + 1, numenin)) == 0.)) cycle
          if (nout(idc, nen) == 0) cycle
          Eev = Ein * 1.e6
          if (Eev <= E6(k, 1)) cycle
          iE = iE + 1
          E6(k, iE) = Eev
          ND(k, iE) = 0
          if (flagbreakup .and. (type == 1 .or. type == 2)) then
            NA(k, iE) = 2
          else
            NA(k, iE) = 1
          endif
!
! Store energy-angle distributions
!
! First we do an extra spectrum normalization
!
          ib = 0
          xsy = 0.
          do nenout = nbegcum(type, nen) - 1, nendcum(type, nen) - 1
            xsy = xsy + (Ehistcum(type, nen, nenout + 1) - Ehistcum(type, nen, max(nenout, 1))) * f0cum(type, nen, max(nenout, 1))
          enddo
          if (xsy > 0.) then
            do nenout = nbegcum(type, nen), nendcum(type, nen)
              f0cum(type, nen, nenout) = f0cum(type, nen, nenout) / xsy
            enddo
          endif
          do nenout = nbegcum(type, nen), nendcum(type, nen)
            ib = ib + 1
            b6(k, iE, ib) = Ehistcum(type, nen, nenout) * 1.e6
            ib = ib + 1
              b6(k, iE, ib) = f0cum(type, nen, nenout) * 1.e-6
            ib = ib + 1
            if (f0cum(type, nen, nenout) == 0.) then
              b6(k, iE, ib) = 0.
            else
              b6(k, iE, ib) = preeqratio(type, nen, nenout)
            endif
            if (flagbreakup .and. (type == 1 .or. type == 2)) then
              ib = ib + 1
              if (f0cum(type, nen, nenout) == 0.) then
                b6(k, iE, ib) = 0.
              else
                b6(k, iE, ib) = buratio(type, nen, nenout)
              endif
            endif
          enddo
          NEP(k, iE) = nendcum(type, nen) - nbegcum(type, nen) + 1
          if (NEP(k, iE) <= 1) then
            NEP(k, iE) = 2
            b6(k, iE, 1) = 0.
            b6(k, iE, 2) = 1.
            if (LANG(k) > 10) then
              b6(k, iE, 3) = 1.
              b6(k, iE, 4) = 0.
              NA(k, iE) = 0
            else
              b6(k, iE, 3) = 0.
              b6(k, iE, 4) = 1.
            endif
            b6(k, iE, 5) = 0.
            b6(k, iE, 6) = 0.
          endif
          NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
        enddo
      endif
!
! Loop over high incident energies
!
      if (flaghigh) then
        if (flagtabddx .and. (type == 1 .or. type == 2) .and. Nddx > 0) LANG(k) = 12
        do nen = 1, Nenspec
          nin = Especindex(nen)
          if (nin < numcut) cycle
          Ein = eninc(nin)
          if (nin == numcut) then
            if (flagendfdet .and. numcut > 1) then
              iE = iE + 1
              E6(k, iE) = Ein * 1.e6 + cuteps
            else
              cycle
            endif
          else
            if (nin < numinc) then
              if (yieldp(type, nin) == 0..and. yieldp(type, nin + 1) == 0.) cycle
            endif
            iE = iE + 1
            E6(k, iE) = Ein * 1.e6
          endif
          ND(k, iE) = 0
!
! Explicit double-differential spectra for deuteron-induced nucleon spectra (break-up effects not covered by Kalbach systematics).
!
          if (flagtabddx .and. (type == 1 .or. type == 2) .and. Nddx > 0) then
            NA(k, iE) = 2 * Nddx
          else
            if (flagbreakup .and. (type == 1 .or. type == 2)) then
              NA(k, iE) = 2
            else
              NA(k, iE) = 1
            endif
          endif
          ib = 0
!
! Store energy-angle distributions
!
          iNEP = 0
          do nenout = nbegcum(type, nen), nendcum(type, nen)
            iNEP = iNEP + 1
            ib = ib + 1
            b6(k, iE, ib) = Ehistcum(type, nen, nenout) * 1.e6
            ib = ib + 1
            b6(k, iE, ib) = f0cum(type, nen, nenout) * 1.e-6
            if (flagtabddx .and. (type == 1 .or. type == 2) .and. Nddx > 0) then
              do i = Nddx, 1, -1
                ib = ib + 1
                b6(k, iE, ib) = rmuddx(i)
                if (f0cum(type, nen, nenout) == 0.) then
                  ib = ib + 1
                  b6(k, iE, ib) = 0.
                else
                  ib = ib + 1
                  b6(k, iE, ib) = f0ddx(type, nen, i, nenout)
                endif
              enddo
            else
              ib = ib + 1
              if (f0cum(type, nen, nenout) == 0.) then
                b6(k, iE, ib) = 0.
              else
                b6(k, iE, ib) = preeqratio(type, nen, nenout)
              endif
              if (flagbreakup .and. (type == 1 .or. type == 2)) then
                ib = ib + 1
                if (f0cum(type, nen, nenout) == 0.) then
                  b6(k, iE, ib) = 0.
                else
                  b6(k, iE, ib) = buratio(type, nen, nenout)
                endif
              endif
            endif
          enddo
          NEP(k, iE) = iNEP
          if (NEP(k, iE) <= 1) then
            NEP(k, iE) = 2
            b6(k, iE, 1) = 0.
            b6(k, iE, 2) = 1.
            if (LANG(k) > 10) then
              b6(k, iE, 3) = 1.
              b6(k, iE, 4) = 0.
              NA(k, iE) = 0
            else
              b6(k, iE, 3) = 0.
              b6(k, iE, 4) = 1.
            endif
            b6(k, iE, 5) = 0.
            b6(k, iE, 6) = 0.
          endif
          NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
        enddo
      endif
!
! ENDF parameters
!
      NE6ea(k) = iE
      NR6ea(k) = 1
      NBT6ea(k, 1) = NE6ea(k)
      INTER6ea(k, 1) = 2
!
! ENDF parameters for MF8
!
      ZAPr(MT, k) = ZAP(k)
      if ( .not. flagrp10) LMF(MT, k) = 6
      LFS(MT, k) = 0
      ELFS(MT, k) = 0.
      MATP(MT, k) = 0
    enddo
    kpart = k
  endif
!
! B. Residual production cross sections
!
  if (flagrp6 .and. MT /= 18) then
    do Zix = numZ, 0, -1
      do Nix = numN, 0, -1
        if ( .not. rpexist(Zix, Nix)) cycle
!
! Total only
!
        iso = 0
        Erpcut = eninccut
        if (Nisorp(Zix, Nix) == 0) then
          k = k + 1
          Z = Zinit - Zix
          A = Ainit - Zix - Nix
          ZAP(k) = 1000. * Z + A
          AWP(k) = real(nucmass(Zix, Nix) / parmass(1))
          Ey(k, 1) = E3(MT, 1)
          Y(k, 1) = 0.
          iE = 1
          if (flaghigh) then
            if (flagendfdet .and. numcut > 1) then
              Ey(k, 2) = eninccut * 1.e6
              Y(k, 2) = 0.
              iE = 2
            endif
            do nin = numcut, numinc
              if (xsnonel(nin) == 0.) cycle
              if (nin < numinc) then
                if (Yrp(Zix, Nix, nin) == 0..and. Yrp(Zix, Nix, nin + 1) == 0.) cycle
              endif
              iE = iE + 1
              if (nin == numcut .and. flagendfdet) then
                Ey(k, iE) = eninc(nin) * 1.e6 + cuteps
              else
                Ey(k, iE) = eninc(nin) * 1.e6
              endif
              Y(k, iE) = Yrp(Zix, Nix, nin)
            enddo
          else
            do nin = 1, numcut
              if (xsnonel(nin) == 0.) cycle
              if (Yrp(Zix, Nix, nin) == 0..and.nin < numcut .and. Yrp(Zix, Nix, min(nin + 1, numenin)) == 0.) cycle
              if (eninc(nin) * 1.e6 < Ey(k, 1)) cycle
              iE = iE + 1
              Ey(k, iE) = eninc(nin) * 1.e6
              if (Ey(k, iE) == Ey(k, iE - 1)) iE = iE - 1
              Y(k, iE) = Yrp(Zix, Nix, nin)
              if (iE == 2 .and. Yrp(Zix, Nix, nin) == 0.) Erpcut = eninc(nin)
            enddo
          endif
          NP6y(k) = iE
          NR6y(k) = 1
          NBT6y(k, 1) = NP6y(k)
          INTER6y(k, 1) = 2
          LEP(k) = 1
          LAW(k) = 0
          LIP(k) = 0
          ZAPr(MT, k) = ZAP(k)
          if ( .not. flagrp10) LMF(MT, k) = 6
          LFS(MT, k) = 0
          ELFS(MT, k) = 0.
          MATP(MT, k) = 0
        else
!
! Ground state and isomers
!
          iso = -1
          do nex = 0, nlevmax
            if ( .not. isorpexist(Zix, Nix, nex)) cycle
            if (iso == (numiso - 1)) cycle
            iso = iso + 1
            k = k + 1
            Z = Zinit - Zix
            A = Ainit - Zix - Nix
            ZAP(k) = 1000. * Z + A
            AWP(k) = real(nucmass(Zix, Nix) / parmass(1))
            Ey(k, 1) = E3(MT, 1)
            Y(k, 1) = 0.
            iE = 1
            if (flaghigh) then
              if (flagendfdet) then
                Ey(k, 2) = eninccut * 1.e6
                Y(k, 2) = 0.
                iE = 2
              endif
              do nin = numcut, numinc
                if (xsnonel(nin) == 0.) cycle
                if (nin < numinc) then
                  if (Yrpiso(Zix, Nix, nex, nin) == 0..and. Yrpiso(Zix, Nix, nex, nin + 1) == 0.) cycle
                endif
                iE = iE + 1
                if (nin == numcut .and. flagendfdet) then
                  Ey(k, iE) = eninc(nin) * 1.e6 + cuteps
                else
                  Ey(k, iE) = eninc(nin) * 1.e6
                endif
                Y(k, iE) = Yrpiso(Zix, Nix, nex, nin)
              enddo
            else
              do nin = 1, numcut
                if (xsnonel(nin) == 0.) cycle
                if (Yrpiso(Zix, Nix, nex, nin) == 0..and.nin < numcut &
 &                .and. Yrpiso(Zix, Nix, nex, min(nin + 1, numenin)) == 0.) cycle
                if (eninc(nin) * 1.e6 < Ey(k, 1)) cycle
                if (iE < numenin) iE = iE + 1
                Ey(k, iE) = eninc(nin) * 1.e6
                Y(k, iE) = Yrpiso(Zix, Nix, nex, nin)
              enddo
            endif
            NP6y(k) = iE
            NR6y(k) = 1
            NBT6y(k, 1) = NP6y(k)
            INTER6y(k, 1) = 2
            LEP(k) = 1
            LAW(k) = 0
            LIP(k) = iso
            ZAPr(MT, k) = ZAP(k)
            if ( .not. flagrp10) LMF(MT, k) = 6
            LFS(MT, k) = nex
            ELFS(MT, k) = Erpiso(Zix, Nix, nex) * 1.e6
            MATP(MT, k) = 0
          enddo
        endif
!
! 1. Recoil energy distributions
!
        if (flagrecoil) then
          if (iso > 0) then
            kbeg = k - iso
          else
            kbeg = k
          endif
          if (kbeg == kpart) kbeg = kbeg + 1
          do kk = kbeg, k
            flagrec(kk) = .true.
            LAW(kk) = 1
            LANG(kk) = 1
            E6(kk, 1) = Ey(k, 1)
            iE = 1
            do j = 1, 2
              ND(kk, j) = 0
              NA(kk, j) = 0
              NEP(kk, j) = 2
              NW(kk, j) = NEP(kk, j) * (NA(kk, j) + 2)
              b6rec(kk, j, 1) = 0.
              b6rec(kk, j, 2) = 1.
              b6rec(kk, j, 3) = 1.
              b6rec(kk, j, 4) = 0.
            enddo
            do nen = 1, Nenspec
              nin = Especindex(nen)
              Ein = eninc(nin)
              if (Ein < Erpcut) cycle
              iE = iE + 1
              E6(kk, iE) = Ein * 1.e6
              ib = 0
              do nenout = nbegcumrec(Zix, Nix, nen), nendcumrec(Zix, Nix, nen)
                ib = ib + 1
                b6rec(kk, iE, ib) = Ehistcumrec(Zix, Nix, nen, nenout) * 1.e6
                ib = ib + 1
                b6rec(kk, iE, ib) = f0cumrec(Zix, Nix, nen, nenout) * 1.e-6
              enddo
              NEP(kk, iE) = nendcumrec(Zix, Nix, nen) - nbegcumrec(Zix, Nix, nen) + 1
              if (NEP(kk, iE) <= 1) then
                NEP(kk, iE) = 2
                b6rec(kk, iE, 1) = 0.
                b6rec(kk, iE, 2) = 1.
                b6rec(kk, iE, 3) = 1.
                b6rec(kk, iE, 4) = 0.
              endif
              NA(kk, iE) = 0
              ND(kk, iE) = 0
              NW(kk, iE) = NEP(kk, iE) * (NA(kk, iE) + 2)
            enddo
            NE6ea(kk) = iE
            NR6ea(kk) = 1
            NBT6ea(kk, 1) = NE6ea(kk)
            INTER6ea(kk, 1) = 2
          enddo
        endif
      enddo
    enddo
  endif
!
! B. Photons
!
! 1. Product yields
!
  k = k + 1
  ZAP(k) = 0.
  AWP(k) = 0.
  if (flaghigh) then
    if (flagendfdet) then
      Ey(k, 1) = E3(MT, 1)
      Y(k, 1) = yieldp(0, 1)
      iE = 1
      if (numcut > 1) then
        Ey(k, 2) = eninccut * 1.e6
        Y(k, 2) = yieldp(0, 1)
        Ey(k, 3) = eninccut * 1.e6 + cuteps
        Y(k, 3) = yieldp(0, numcut)
        iE = 3
      endif
    endif
    do nin = numcut + 1, numinc
      if (nin < numinc) then
        if (yieldp(0, nin) == 0..and. yieldp(0, nin + 1) == 0.) cycle
      endif
      iE = iE + 1
      Ey(k, iE) = eninc(nin) * 1.e6
      Y(k, iE) = yieldp(0, nin)
    enddo
  else
    iE = 0
    do nin = 1, numcut
      if (yieldp(0, nin) == 0..and. nin < numcut .and. yieldp(0, min(nin + 1, numenin)) == 0.) cycle
      if (xsany(nin) == 0..and.xsany(min(nin + 1, numenin)) == 0.) cycle
      iE = iE + 1
      Ey(k, iE) = eninc(nin) * 1.e6
      Y(k, iE) = yieldp(0, nin)
    enddo
  endif
  NP6y(k) = iE
  NR6y(k) = 1
  NBT6y(k, 1) = NP6y(k)
  INTER6y(k, 1) = 2
  LAW(k) = 1
  LIP(k) = 0
  LANG(k) = 1
  LEP(k) = 1
!
! 2. Energy-angle distributions
!
  if (flaghigh) then
    E6(k, 1) = Ey(k, 1)
    E6(k, 2) = eninccut * 1.e6
    do j = 1, 2
      ND(k, j) = 0
      NA(k, j) = 0
      NEP(k, j) = 2
      NW(k, j) = NEP(k, j) * (NA(k, j) + 2)
      b6gam(j, 1) = 0.
      b6gam(j, 2) = 1.
      b6gam(j, 3) = 1.
      b6gam(j, 4) = 0.
    enddo
    if (flagendfdet .and. numcut > 1) then
      iE = 2
    else
      iE = 1
    endif
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (nin < numcut) cycle
      Ein = eninc(nin)
      if (nin == numcut) then
        if (flagendfdet .and. numcut > 1) then
          iE = iE + 1
          E6(k, iE) = enencut * 1.e6 + cuteps
        else
          cycle
        endif
      else
        if (nin < numinc) then
          if (yieldp(0, nin) == 0..and. yieldp(0, nin + 1) == 0.) cycle
        endif
        iE = iE + 1
        E6(k, iE) = Ein * 1.e6
      endif
      ND(k, iE) = 0
      NA(k, iE) = 0
      ib = 0
!
! Store energy-angle distributions
!
      do nenout = nbegcum(0, nen), nendcum(0, nen)
        ib = ib + 1
        b6gam(iE, ib) = Ehistcum(0, nen, nenout) * 1.e6
        ib = ib + 1
        b6gam(iE, ib) = f0cum(0, nen, nenout) * 1.e-6
      enddo
      NEP(k, iE) = nendcum(0, nen) - nbegcum(0, nen) + 1
      if (NEP(k, iE) <= 0 .or. nendcum(0, nen) == 0) then
        NEP(k, iE) = 2
        b6gam(iE, 1) = 0.
        b6gam(iE, 2) = 1.
        b6gam(iE, 3) = 1.
        b6gam(iE, 4) = 0.
      endif
      NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
    enddo
  else
    E6(k, 1) = Ey(k, 1)
    ND(k, 1) = 0
    NA(k, 1) = 0
    NEP(k, 1) = 2
    NW(k, 1) = NEP(k, 1) * (NA(k, 1) + 2)
    b6gam(1, 1) = 0.
    b6gam(1, 2) = 1.
    b6gam(1, 3) = 1.
    b6gam(1, 4) = 0.
    iE = 1
    do nen = 1, Nenspec
      nin = Especindex(nen)
      if (nin > numcut) cycle
      Ein = eninc(nin)
      if (yieldp(0, nin) == 0..and. yieldp(0, min(nin + 1, numenin)) == 0.) cycle
      if (xsany(nin) == 0..and.xsany(min(nin + 1, numenin)) == 0.) cycle
      iE = iE + 1
      E6(k, iE) = Ein * 1.e6
      ND(k, iE) = 0
      NA(k, iE) = 0
      ib = 0
      do nenout = nbegcum(0, nen), nendcum(0, nen)
        ib = ib + 1
        b6gam(iE, ib) = Ehistcum(0, nen, nenout) * 1.e6
        ib = ib + 1
        b6gam(iE, ib) = f0cum(0, nen, nenout) * 1.e-6
      enddo
      NEP(k, iE) = nendcum(0, nen) - nbegcum(0, nen) + 1
      if (NEP(k, iE) <= 0 .or. nendcum(0, nen) == 0) then
        NEP(k, iE) = 2
        b6gam(iE, 1) = 0.
        b6gam(iE, 2) = 1.
        b6gam(iE, 3) = 1.
        b6gam(iE, 4) = 0.
      endif
      NW(k, iE) = NEP(k, iE) * (NA(k, iE) + 2)
    enddo
  endif
  NE6ea(k) = iE
  NR6ea(k) = 1
  NBT6ea(k, 1) = NE6ea(k)
  INTER6ea(k, 1) = 2
  NK(6, MT) = k
!
! ENDF-6 parameters
!
! make3clean: subroutine to clean up MF6
! write6    : subroutine to write MF6
!
  if (k /= 0) mtexist(6, MT) = .true.
  NK(6, MT) = k
  if (NK(6, MT) /= 0) mfexist(6) = .true.
  if (flagrecoil .and. flaghigh) then
    LCT = 3
  else
    LCT = 2
  endif
  if (mtexist(6, MT)) then
    if (flagclean) call make6clean(MT)
    call write6(MT)
  endif
  return
end subroutine make6mt5
! Copyright A.J. Koning 2021
