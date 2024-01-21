subroutine make9_10
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF9 and MF10
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
!   numiso        ! number of isomers
!   nummt         ! number of MT numbers
!   numN          ! maximal number of neutrons from initial compound nucleus
!   numsec        ! number of sections
!   numZ          ! maximal number of protons from initial compound nucleus
! Variables for ENDF limits, switches and tolerances
!   cuteps        ! energy shift at MT5 cutoff energy (in eV)
! Variables for input of ENDF structure
!   flagfis10     ! flag to put (subactinide) fission cross sections in MF10
!   flagrp10      ! flag to put residual production cross sections in MF10
! Variables for input of ENDF library type
!   flagclean     ! flag to clean up double points
!   flaghigh      ! flag for high energies ( > 20 MeV)
! Constants
!   xsepslow      ! lower limit for cross sections in millibarns
! Variables for reaction initialization
!   eninccut      ! last incident energy before high - energy format
!   nlevmax       ! number of included discrete levels
!   numcut        ! number of energies before high - energy format
! Variables for info from TALYS
!   Ainit         ! mass number of initial compound nucleus
!   eninc         ! incident energy
!   k0            ! index of incident particle
!   Lisomer       ! isomeric number of target
!   numinc        ! number of incident energies
!   Zinit         ! charge number of initial compound nucleus
! Variables for initialization of ENDF format
!   idnum         ! number of different exclusive cross sections
!   mfexist       ! flag for existence of MF - number
!   mtexist       ! flag for existence of MT - number
!   MTid          ! channel identifier for MT - number
!   MTinel        ! MT - number for inelastic scattering
! Variables for channel cross sections in ENDF format
!   branchiso     ! branching ratio for isomer
!   Ethexcliso    ! threshold energy for isomer
!   idchannel     ! identifier for channel
!   isoexist      ! flag for existence of isomer
!   Nisomer       ! number of isomers
!   Qexcl         ! Q - value
!   Qexcliso      ! Q - value for isomer
!   xsexcl        ! exclusive cross section
!   xsexcliso     ! exclusive cross section for isomer
! Variables for discrete state cross sections in ENDF format
!   Qdisc         ! Q - value
! Variables for residual production cross sections in ENDF format
!   Erpiso        ! energy of isomer
!   Ethrp         ! threshold energy for residual product
!   Ethrpiso      ! threshold energy for isomer of residual product
!   isorpexist    ! flag for existence of isomer of residual nuclide
!   Qrp           ! Q - value for residual product
!   Qrpiso        ! Q - value for isomer of residual product
!   rpexist       ! flag for existence of residual nuclide
!   xsrp          ! residual production cross section
!   xsrpiso       ! residual production cross section for isomer
! Variables for MF1
!   EMAX          ! upper limit of energy range for evaluation
! Variables for MF3
!   QM            ! Q - value (in ENDF - 6 format)
! Variables for MF8_10
!   E10           ! incident energy (in ENDF - 6 format)
!   E10ZA         ! incident energy (in ENDF - 6 format)
!   ELFS          ! excitation energy of final state
!   ErpZAiso      ! energy of isomer
!   EthZA         ! threshold energy
!   INTER10       ! interpolation scheme
!   INTERZA       ! interpolation scheme
!   IZAP          ! second IZAP - number
!   LFS           ! final state number
!   LFSZA         ! final state number
!   NBT10         ! separation value for interpolation scheme
!   NBTZA         ! separation value for interpolation scheme
!   NP10          ! number of incident energies
!   NPZA          ! number of incident energies
!   NR10          ! number of interpolation ranges
!   NRZA          ! number of interpolation ranges
!   NSt           ! number of final states
!   NZA           ! number of nuclides
!   QIiso         ! Q - value for isomer (in ENDF - 6 format)
!   QIZA          ! Q - value (in ENDF - 6 format)
!   QMZA          ! Q - value (in ENDF - 6 format)
!   XMFZA         ! second MF - number
!   xsiso         ! cross section for isomer (in ENDF - 6 format)
!   xsrpZA        ! cross section for residual production (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  logical   :: yield            ! logical to use MF9 or MF10
  integer   :: A                ! mass number of target nucleus
  integer   :: idc              ! help variable
  integer   :: iE               ! energy counter
  integer   :: iso              ! counter for isomer
  integer   :: iza              ! counter for Z,A combinations
  integer   :: L                ! counter for Legendre coefficients
  integer   :: MT               ! MT-number
  integer   :: nex              ! discrete level
  integer   :: nin              ! counter for incident energy
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: Z                ! charge number of target nucleus
  integer   :: Zix              ! charge number index for residual nucleus
  real(sgl) :: eps              ! help variable
!
! **************************** Make MF10 *******************************
!
!
  do MT = 1, nummt
    if (flagfis10 .and. MT == 18) goto 300
    if ( .not. mtexist(3, MT)) cycle
!
! Residual production cross sections
!
    if (flagrp10 .and. MT == 5) then
      iza = 0
      do Zix = numZ, 0, -1
        do Nix = numN, 0, -1
          if ( .not. rpexist(Zix, Nix)) cycle
          if (Zix == 0 .and. Nix == 0) cycle
          Z = Zinit - Zix
          A = Ainit - Zix - Nix
          do nex = 0, nlevmax
            if (isorpexist(Zix, Nix, nex)) goto 150
          enddo
          if (iza > numsec) goto 200
          iza = iza + 1
          IZAP(iza) = 1000 * Z + A
          XMFZA(iza) = 0
          LFSZA(iza) = 0
          QMZA(iza) = Qrp(Zix, Nix) * 1.e6
          QIZA(iza) = Qrp(Zix, Nix) * 1.e6
          EthZA(iza) = Ethrp(Zix, Nix)
          E10ZA(iza, 1) = Ethrp(Zix, Nix) * 1.e6
          xsrpZA(iza, 1) = 0.
          iE = 1
          do nin = 1, numinc
            if (xsrp(Zix, Nix, nin) == 0.) cycle
            if (nin > 1) iE = iE + 1
            E10ZA(iza, iE) = eninc(nin) * 1.e6
            xsrpZA(iza, iE) = xsrp(Zix, Nix, nin) * 1.e-3
          enddo
          NPZA(iza) = iE
          NRZA(iza) = 1
          NBTZA(iza, 1) = iE
          INTERZA(iza, 1) = 2
          if (iE <= 1) iza = iza - 1
          cycle
  150         iso = -1
          do nex = 0, nlevmax
            if ( .not. isorpexist(Zix, Nix, nex)) cycle
            if (iza > numsec) goto 200
            iza = iza + 1
            iso = iso + 1
            IZAP(iza) = 1000 * Z + A
            LFSZA(iza) = nex
            QMZA(iza) = Qrp(Zix, Nix) * 1.e6
            QIZA(iza) = Qrpiso(Zix, Nix, nex) * 1.e6
            EthZA(iza) = Erpiso(Zix, Nix, nex)
            ErpZAiso(iza) = Erpiso(Zix, Nix, nex) * 1.e6
            E10ZA(iza, 1) = Ethrpiso(Zix, Nix, nex) * 1.e6
            xsrpZA(iza, 1) = 0.
            iE = 1
            do nin = 1, numinc
              if (xsrpiso(Zix, Nix, nex, nin) == 0.) cycle
              if (nin > 1) iE = iE + 1
              E10ZA(iza, iE) = eninc(nin) * 1.e6
              xsrpZA(iza, iE) = xsrpiso(Zix, Nix, nex, nin) * 1.e-3
            enddo
            NPZA(iza) = iE
            NRZA(iza) = 1
            NBTZA(iza, 1) = iE
            INTERZA(iza, 1) = 2
            if (iE <= 1) iza = iza - 1
          enddo
        enddo
      enddo
  200     NZA = iza
      NSt(MT) = iza
      mtexist(10, MT) = .true.
      mfexist(10) = .true.
      cycle
    endif
    if (MTid(MT) ==  -1) cycle
    if (MT == 91) cycle
    if (MT == 649) cycle
    if (MT == 699) cycle
    if (MT == 749) cycle
    if (MT == 799) cycle
    if (MT == 849) cycle
!
! Isomeric exclusive cross sections
!
    do idc = 0, idnum
      if (idchannel(idc) == MTid(MT)) then
        if (k0 == 1 .and. MT == 102) then
          yield = .true.
        else
          yield = .false.
        endif
        iso = 0
        do nex = 0, nlevmax
          if ( .not. isoexist(idc, nex)) cycle
          if (Ethexcliso(idc, nex) >= eninccut) cycle
          if (iso == numiso) cycle
          iso = iso + 1
          QM(MT) = Qexcl(idc) * 1.e6
          if ( k0 > 0 .and. MT == MTinel) then
            if (nex == 0 .and. Lisomer > 0) then
              QIiso(MT, iso) = Qdisc(k0, 0) * 1.e6
            else
              QIiso(MT, iso) = 0.
            endif
            if (nex > 0 .and. Lisomer == 0) QIiso(MT, iso) = Qdisc(k0, nex) * 1.e6
          else
            QIiso(MT, iso) = Qexcliso(idc, nex) * 1.e6
          endif
!
! Make sure the energy level is precisely the difference of Q-values
!
          ELFS(MT, iso) = QM(MT) - QIiso(MT, iso)
          if (nex > 0.) then
            L = int(log10(ELFS(MT, iso)))
            eps = 9. * 10 **real(L - 7)
            ELFS(MT, iso) = real(int(ELFS(MT, iso) + eps))
          endif
          if (k0 == 0 .and. Qexcl(idc) >= 0.) then
            iE = 0
          else
            if (MT == MTinel .and. nex == 0) then
              E10(MT, iso, 1) = EmineV
            else
              E10(MT, iso, 1) = max(Ethexcliso(idc, nex) * 1.e6, EmineV)
            endif
            if (MT == MTinel .and. nex > 0 .and. Lisomer > 0) E10(MT, iso, 1) = EmineV
            if (eninc(1) == EminMeV) then
              if (yield) then
                if (xsexcl(idc, 1) <= xsepslow) then
                  if (iso == 1) then
                    xsiso(MT, iso, 1) = 1.
                  else
                    xsiso(MT, iso, 1) = 0.
                  endif
                else
                  xsiso(MT, iso, 1) = branchiso(idc, nex, 1)
                endif
              else
                xsiso(MT, iso, 1) = xsexcliso(idc, nex, 1) * 1.e-3
              endif
            else
              if (yield) then
                xsiso(MT, iso, 1) = 1. / Nisomer(idc)
              else
                xsiso(MT, iso, 1) = 0.
              endif
            endif
            iE = 1
          endif
          do nin = 1, numcut
            if (eninc(nin) == EminMeV) cycle
            if (eninc(nin) <= Ethexcliso(idc, nex)) cycle
            if (xsexcl(idc, nin) == 0.) cycle
            iE = iE + 1
            E10(MT, iso, iE) = eninc(nin) * 1.e6
            if (yield) then
              if (xsexcl(idc, nin) <= xsepslow) then
                if (iso == 1) then
                  xsiso(MT, iso, iE) = 1.
                else
                  xsiso(MT, iso, iE) = 0.
                endif
              else
                xsiso(MT, iso, iE) = branchiso(idc, nex, nin)
              endif
            else
              xsiso(MT, iso, iE) = xsexcliso(idc, nex, nin) * 1.e-3
            endif
          enddo
          if (iE <= 1) then
            iso = iso - 1
            cycle
          endif
          if (xsexcliso(idc, nex, numcut) == 0.) then
            iE = iE + 1
            E10(MT, iso, iE) = eninccut * 1.e6
            if (yield) then
              xsiso(MT, iso, iE) = 1. / Nisomer(idc)
            else
              xsiso(MT, iso, iE) = 0.
            endif
          endif
!
! High energies
!
          if (flaghigh) then
            iE = iE + 1
            E10(MT, iso, iE) = eninccut * 1.e6 + cuteps
            if (yield) then
              xsiso(MT, iso, iE) = 1. / Nisomer(idc)
            else
              xsiso(MT, iso, iE) = 0.
            endif
            iE = iE + 1
            E10(MT, iso, iE) = EMAX
            if (yield) then
              xsiso(MT, iso, iE) = 1. / Nisomer(idc)
            else
              xsiso(MT, iso, iE) = 0.
            endif
          endif
!
! ENDF-6 parameters
!
          LFS(MT, iso) = nex
          NP10(MT, iso) = iE
          NR10(MT, iso) = 1
          NBT10(MT, iso, 1) = iE
          INTER10(MT, iso, 1) = 2
          if (yield) then
            mtexist(9, MT) = .true.
            mfexist(9) = .true.
          else
            mtexist(10, MT) = .true.
            mfexist(10) = .true.
          endif
        enddo
        if ((k0 == 1 .and. MT == 4) .or. (k0 == 2 .and. MT == 103) .or. (k0 == 3 .and. MT == 104) .or. &
 &        (k0 == 4 .and. MT == 105) .or. (k0 == 5 .and. MT == 106) .or. (k0 == 6 .and. MT == 107)) then
          ELFS(MT, 1) = 0.
        endif
        NSt(MT) = iso
        cycle
      endif
    enddo
    cycle
!
! Fission (usually only for subactinides)
!
  300   iso = 1
    iE = 0
    do nin = 1, numinc
      iE = iE + 1
      E10(MT, iso, iE) = eninc(nin) * 1.e6
      xsiso(MT, iso, iE) = xsexcl(-1, nin) * 1.e-3
    enddo
    Nst(MT) = 1
    QM(MT) = 2.e8
    QIiso(MT, iso) = 0.
    ELFS(MT, iso) = QM(MT)
    LFS(MT, iso) = 0
    NP10(MT, iso) = iE
    NR10(MT, iso) = 1
    NBT10(MT, iso, 1) = iE
    INTER10(MT, iso, 1) = 2
    mtexist(10, MT) = .true.
    mfexist(10) = .true.
  enddo
!
! Clean up
!
! make9_10clean: subroutine to clean up MF9 and MF10
!
  if (flagclean) call make9_10clean
  return
end subroutine make9_10
! Copyright A.J. Koning 2021
