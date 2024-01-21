subroutine make3total
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3 for total, elastic and reaction cross sections
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
!   numen6         ! number of incident energies
!   numenin        ! number of incident energies
! Variables for input of specific ENDF data
!   adopt          ! logical for existence of MF information (per MT)
!   Eahigh         ! upper energy of MT values to be adopted
!   Ealow          ! lower energy of MT values to be adopted
! Variables for input of ENDF library type
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps         ! energy shift at MT5 cutoff energy (in eV)
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
! Variables for reaction initialization
!   EHres          ! upper energy in resonance range
! Variables for info from TALYS
!   eninc          ! incident energy
!   k0             ! index of incident particle
!   numinc         ! number of incident energies
! Variables for initialization of ENDF format
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
! Variables for total cross sections in ENDF format
!   eninc6         ! incident energy for total cross section
!   numcut6        ! number of energies before high - energy format for total cross sections
!   numinc6        ! number of incident energies for total cross sections
!   xselas         ! total elastic cross section
!   xselas6        ! total elastic cross section
!   xsnonel        ! nonelastic cross section
!   xsnonel6       ! nonelastic cross section
!   xstot6         ! total cross section
! Variables for partial cross sections in ENDF format
!   xsnonth        ! sum of all non - thr. reactions except (n, g) and (n, f)
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF3
!   E3             ! incident energy for MF3 (in ENDF - 6 format)
!   E3adopt        ! incident energy adopted from other library
!   E3res          ! energy in resonance range
!   EthMT          ! threshold energy
!   LR3            ! breakup flag
!   NE3adopt       ! number of adopted incident energies
!   NE3res         ! number of energies in resonance range
!   NMT            ! total number of MT sections
!   QI             ! Q - value (in ENDF - 6 format)
!   QM             ! Q - value (in ENDF - 6 format)
!   xs             ! cross section
!   xs3adopt       ! cross section adopted from other library
!   xs3res         ! cross section in resonance range
!
! *** Declaration of local data
!
  implicit none
  logical   :: click                          ! logical to confirm existence
  integer   :: i                              ! counter
  integer   :: iE                             ! energy counter
  integer   :: MF                             ! MF-number
  integer   :: MT                             ! MT-number
  integer   :: N                              ! neutron number of residual nucleus
  integer   :: nen                            ! energy counter
  integer   :: nen1                           ! energy counter
  integer   :: nin                            ! counter for incident energy
  integer   :: nin2                           ! energy counter
  integer   :: Nlin                           ! number of linearized points
  real(sgl) :: dxs                            ! uncertainty of experimental cross section
  real(sgl) :: E                              ! incident energy
  real(sgl) :: ea                             ! help variable
  real(sgl) :: Eadd(10*numenin)               ! energy (to add)
  real(sgl) :: eb                             ! help variable
  real(sgl) :: Eback(0:numen6)                ! energy for background cross section
  real(sgl) :: ee                             ! energy
  real(sgl) :: eemev                          ! energy in MeV
  real(sgl) :: Eev                            ! energy in eV
  real(sgl) :: Efac                           ! help variable
  real(sgl) :: Eprev                          ! previous energy
  real(sgl) :: fac                            ! factor
  real(sgl) :: xsadd(10 * numenin)            ! cross section (to add)
  real(sgl) :: xsb                            ! help variable
  real(sgl) :: xsback(0:numen6)               ! background cross section
  real(sgl) :: xse                            ! experimental cross section
  real(sgl) :: xsint                          ! interpolated cross section
  real(sgl) :: xst                            ! help variable
  real(sgl) :: xstot                          ! total cross section (neutrons only) for
  real(sgl) :: y1                             ! coordinates of the 1st summit of the triangle
  real(sgl) :: y2                             ! coordinates of the 2nd summit of the triangle
!
! *********************** Make MF3 for MT1,2,3 *************************
!
  MF = 3
  NMT = 0
  do MT = 1, 3
    EthMT(MT) = EminMeV
    QM(MT) = 0.
    QI(MT) = 0.
    iE = 0
!
! Photonuclear reactions: we only use MT3
!
    if (k0 == 0) then
      if (MT <= 2) cycle
      do nin = 1, numcut6
        if (xsnonel(nin) == 0.) cycle
        iE = iE + 1
        E3(MT, iE) = eninc(nin) * 1.e6
        xs(MT, iE) = real(xsnonel(nin)) * 1.e-3
      enddo
      if (flaghigh) then
        iE = iE + 1
        E3(MT, iE) = eninc6(numcut6) * 1.e6 + cuteps
        xs(MT, iE) = real(xsnonel6(numcut6)) * 1.e-3
        do nin = numcut6 + 1, numinc6
          if (xsnonel6(nin) == 0.) cycle
          iE = iE + 1
          E3(MT, iE) = eninc6(nin) * 1.e6
          xs(MT, iE) = real(xsnonel6(nin)) * 1.e-3
        enddo
      endif
    endif
!
! Charged-particle induced reactions: we do not use MT1
!
    if (k0 > 1) then
      if (MT == 1) cycle
      if (MT == 2) then
        do nin = 1, numinc
          if (xselas(nin) == 0..and.xselas(nin - 1) /= 0.) then
            iE = iE + 1
            E3(MT, iE) = eninc(numinc) * 1.e6
            xs(MT, iE) = real(xselas(nin - 1)) * 1.e-3
            exit
          else
            if (xselas(nin) == 0.) cycle
            iE = iE + 1
            E3(MT, iE) = eninc(nin) * 1.e6
            xs(MT, iE) = real(xselas(nin)) * 1.e-3
          endif
        enddo
      endif
      if (MT == 3) then
        if (flagendfdet .and. Eswitch > 0.) then
          iE = 0
          do nin = 1, numinc
            if (xsnonel(nin) == 0.) cycle
            iE = iE + 1
            E3(MT, iE) = eninc(nin) * 1.e6
            xs(MT, iE) = real(xsnonel(nin)) * 1.e-3
          enddo
        else
          cycle
        endif
      endif
    endif
!
! Neutron induced reactions
!
! Adopt background cross sections for resonance range or set them to 0
!
    if (k0 == 1) then
      if (NE3res(MT) /= 0) then
        do i = 1, NE3res(MT)
          E3(MT, i) = E3res(MT, i)
          xs(MT, i) = xs3res(MT, i)
          Eback(i) = E3(MT, i)
          xsback(i) = xs(MT, i)
        enddo
        iE = NE3res(MT)
      else
        EthMT(MT) = EminMeV
        if (EHres > 0.) then
          E3(MT, 1) = EthMT(MT) * 1.e6
          xs(MT, 1) = 0.
          E3(MT, 2) = EHres
          xs(MT, 2) = 0.
          iE = 2
        endif
      endif
!
! Add exothermic partial cross sections other than capture and fission to nonelastic and total cross section.
!
      if ((MT == 1 .or. MT == 3) .and. xsnonth(1) > 0.) then
        Nlin = 100
        call locate(eninc, 1, numinc, EHres * 1.e-6, nen)
        if (nen >= 1) then
          nen1 = min(nen + 1, numinc)
          call linear(eninc, xsnonth, nen1, Eadd, xsadd, Nlin)
          iE = 0
          do i = 1, Nlin
            E = Eadd(i) * 1.e6
            N = NE3res(MT)
            xsb = 0.
            if (N > 0) then
              if (E >= Eback(1) .and. E <= Eback(N)) then
                call locate(Eback, 1, N, E, nen)
                call pol1(Eback(nen), Eback(nen + 1), xsback(nen), xsback(nen + 1), E, xsb)
              endif
            endif
            if (E <= EHres) then
              iE = iE + 1
              E3(MT, iE) = E
              xs(MT, iE) = (xsb + xsadd(i)) * 1.e-3
            endif
          enddo
          fac = (EHres - E3(MT, iE)) / (E3(MT, iE) - E3(MT, iE - 1))
          dxs = xs(MT, iE) - xs(MT, iE - 1)
          xse = xs(MT, iE) + fac * dxs
          iE = iE + 1
          E3(MT, iE) = Ehres
          if (xse >= 0.) then
            xs(MT, iE) = xse
          else
            xs(MT, iE) = xs(MT, iE - 1)
          endif
        endif
      endif
      eninc6(0) = 0.
      Eprev = 0.
      click = .true.
      do nin = 1, numcut6
        Eev = eninc6(nin) * 1.e6
        if (nin > 1 .and. Eprev < EHres .and. Eev - EHres >= 1.e-6 .and. E3(MT, iE) + 1.e-3 >= E3(MT, max(iE - 1, 0))) then
          Efac = (EHres - Eprev) / (Eev - Eprev)
          if (MT == 1) xst = real(xstot6(nin - 1)) + Efac * real(xstot6(nin) - xstot6(nin - 1))
          if (MT == 2) xst = real(xselas6(nin - 1)) + Efac * real(xselas6(nin) - xselas6(nin - 1))
          if (MT == 3) xst = real(xsnonel6(nin - 1)) + Efac * real(xsnonel6(nin) - xsnonel6(nin - 1))
          if (xst > 0.) then
            iE = iE + 1
            E3(MT, iE) = EHres
            xs(MT, iE) = xst * 1.e-3
          endif
        endif
        Eprev = Eev
        if (Eev < Ehres) cycle
        if (xselas6(nin) == 0.) cycle
!
! Include cross sections from other datafile
!
        if (click .and. adopt(3, 1) .and. Eev >= Ealow(3, 1)) then
          do nin2 = 1, NE3adopt(1)
            iE = iE + 1
            ee = E3adopt(1, nin2)
            xstot = xs3adopt(1, nin2)
            if (MT == 1) then
              E3(1, iE) = ee
              xs(1, iE) = xstot
            endif
            eemev = ee * 1.e-6
            call locate(eninc6, 1, numinc6 - 1, eemev, nen)
            ea = eninc6(nen)
            eb = eninc6(nen + 1)
            y1 = real(xsnonel6(nen)) * 1.e-3
            y2 = real(xsnonel6(nen + 1)) * 1.e-3
            call pol1(ea, eb, y1, y2, eemev, xsint)
            if (MT == 2) then
              E3(2, iE) = ee
              xs(2, iE) = xstot - xsint
            endif
            if (MT == 3) then
              E3(3, iE) = ee
              xs(3, iE) = xsint
            endif
          enddo
          click = .false.
        else
          if (adopt(3, 1) .and. Eev >= Ealow(3, 1) .and. Eev < Eahigh(3, 1)) cycle
          iE = iE + 1
          E3(MT, iE) = Eev
          if (MT == 1) xs(1, iE) = real(xstot6(nin)) * 1.e-3
          if (MT == 2) xs(2, iE) = real(xselas6(nin)) * 1.e-3
          if (MT == 3) xs(3, iE) = real(xsnonel6(nin)) * 1.e-3
        endif
      enddo
!
! High energies
!
      if (flaghigh) then
        if (eninc6(numcut6) > E3adopt(1, NE3adopt(1)) * 1.e-6) then
          iE = iE + 1
          E3(MT, iE) = eninc6(numcut6) * 1.e6 + cuteps
          if (MT == 1) xs(1, iE) = real(xstot6(numcut6)) * 1.e-3
          if (MT == 2) xs(2, iE) = real(xselas6(numcut6)) * 1.e-3
          if (MT == 3) xs(3, iE) = real(xsnonel6(numcut6)) * 1.e-3
        endif
        do nin = numcut6 + 1, numinc6
          if (xselas6(nin) == 0.) cycle
          if (eninc6(nin) * 1.e6 < EHres) cycle
          if (eninc6(nin) < E3adopt(1, NE3adopt(1)) * 1.e-6) cycle
          iE = iE + 1
          E3(MT, iE) = eninc6(nin) * 1.e6
          if (MT == 1) xs(1, iE) = real(xstot6(nin)) * 1.e-3
          if (MT == 2) xs(2, iE) = real(xselas6(nin)) * 1.e-3
          if (MT == 3) xs(3, iE) = real(xsnonel6(nin)) * 1.e-3
        enddo
      endif
    endif
!
! ENDF-6 parameters
!
    NP(MF, MT) = iE
    NR(MF, MT) = 1
    NBT(MF, MT, 1) = iE
    INTER(MF, MT, 1) = 2
    LR3(MT) = 0
    mtexist(MF, MT) = .true.
    mfexist(MF) = .true.
    NMT = NMT + 1
  enddo
  return
end subroutine make3total
! Copyright A.J. Koning 2021
