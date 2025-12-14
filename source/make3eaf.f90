subroutine make3eaf
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF3 for isomeric cross sections in EAF-format
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
!   numenin       ! number of incident energies
!   nummt         ! number of MT numbers
! Constants
!   parsym        ! symbol of particle
! Variables for reaction initialization
!   EHres         ! upper energy in resonance range
!   eninccut      ! last incident energy before high - energy format
!   nlevmax       ! number of included discrete levels
!   nuclid        ! nuclide
!   numcut        ! number of energies before high - energy format
! Variables for info from TALYS
!   Atarget       ! mass number of nucleus
!   eninc         ! incident energy
!   k0            ! index of incident particle
!   Lisomer       ! isomeric number of target
!   numinc        ! number of incident energies
! Variables for initialization of ENDF format
!   idnum         ! number of different exclusive cross sections
!   mtexist       ! flag for existence of MT - number
!   MTid          ! channel identifier for MT - number
!   MTinel        ! MT - number for inelastic scattering
! Variables for partial cross sections in ENDF format
!   branchiso     ! branching ratio for isomer
!   Ethexcliso    ! threshold energy for isomer
!   idchannel     ! identifier for channel
!   isoexist      ! flag for existence of isomer
!   Qexcl         ! Q - value
!   Qexcliso      ! Q - value for isomer
!   reacstring    ! string with reaction information
!   xsexcliso     ! exclusive cross section for isomer
! Variables for ENDF format
!   INTER         ! interpolation scheme
!   NBT           ! separation value for interpolation scheme
!   NP            ! number of incident energies
!   NR            ! number of interpolation ranges
! Variables for MF3
!   E3            ! incident energy for MF3 (in ENDF - 6 format)
!   E3adopt       ! incident energy adopted from other library
!   eafstring     ! string with reaction information for EAF format
!   LFS3          ! isomeric state number (EAF only)
!   LR3           ! breakup flag
!   mtstring      ! string with reaction information
!   NE3adopt      ! number of adopted incident energies
!   QI            ! Q - value (in ENDF - 6 format)
!   QM            ! Q - value (in ENDF - 6 format)
!   xs            ! cross section
!   xs3adopt      ! cross section adopted from other library
!
! *** Declaration of local data
!
  implicit none
  logical          :: flagiso               ! logical for isomer
  character(len=3) :: str                   ! input line
  integer          :: i                     ! counter
  integer          :: ibeg                  ! index to mark begin of word
  integer          :: idc                   ! help variable
  integer          :: iE                    ! energy counter
  integer          :: iend                  ! index to mark end of word
  integer          :: ilength               ! length of string
  integer          :: iMT                   ! counter for MT number
  integer          :: iso                   ! counter for isomer
  integer          :: MF                    ! MF-number
  integer          :: MT                    ! MT-number
  integer          :: nen                   ! energy counter
  integer          :: nen1                  ! energy counter
  integer          :: nex                   ! discrete level
  integer          :: nin                   ! counter for incident energy
  integer          :: Nlin                  ! number of linearized points
  real(sgl)        :: Eadd(10*numenin)      ! energy (to add)
  real(sgl)        :: ee                    ! energy
  real(sgl)        :: Eev                   ! energy in eV
  real(sgl)        :: Ehigh                 ! help variable
  real(sgl)        :: Qval                  ! Q-value
  real(sgl)        :: x(numenin)            ! help variable
  real(sgl)        :: xsadd(10*numenin)     ! cross section (to add)
!
! ************************* Make MF3 for isomers ***********************
!
  MF = 3
  do iMT = 1, nummt
    if (MTid(iMT) ==  -1) cycle
    if (iMT == 91) cycle
    if (iMT >= 600) cycle
    do idc = -1, idnum
      MT = iMT
      if (idchannel(idc) == MTid(iMT)) then
        iso = -1
        mtstring(MT) = reacstring(idc, 0)
        do nex = 0, nlevmax
!
! Fission
!
          if (idc ==  -1) then
            if (nex == 0) then
              mtstring(MT) = reacstring(idc, 0)
              goto 100
            else
              cycle
            endif
          endif
!
! Other channels
!
          flagiso = .true.
          if ( .not. isoexist(idc, nex)) flagiso = .false.
          if (Ethexcliso(idc, nex) >= eninccut) flagiso = .false.
          if (flagiso) iso = iso + 1
          if (iso > 2) flagiso = .false.
          if (iMT == MTinel .and. Lisomer == iso) flagiso = .false.
          if (flagiso) then
            MT = 300 * iso + iMT
            mtstring(MT) = reacstring(idc, nex)
            QM(MT) = Qexcl(idc) * 1.e6
            Qval = Qexcliso(idc, nex)
            QI(MT) = Qval * 1.e6
            E3(MT, 1) = max(Ethexcliso(idc, nex) * 1.e6, EmineV)
            if (eninc(1) == EminMeV) then
              xs(MT, 1) = xsexcliso(idc, nex, 1) * 1.e-3
            else
              xs(MT, 1) = 0.
            endif
            iE = 1
!
! Linearize exothermic partial cross sections
!
            Ehigh = 0.
            if (Qval > 0..and.iMT /= 102) then
              if (EHres > 0.) then
                Ehigh = EHres
              else
                Ehigh = 1000.
              endif
              call locate(eninc, 1, numinc, Ehigh * 1.e-6, nen)
              if (nen >= 1) then
                do i = 1, nen + 1
                  x(i) = xsexcliso(idc, nex, i)
                enddo
                Nlin = 100
                nen1 = min(nen + 1, numinc)
                call linear(eninc, x, nen1, Eadd, xsadd, Nlin)
                iE = 0
                do i = 1, Nlin
                  if (Eadd(i) <= Ehigh * 1.e-6) then
                    iE = iE + 1
                    E3(MT, iE) = Eadd(i) * 1.e6
                    xs(MT, iE) = xsadd(i) * 1.e-3
                  endif
                enddo
              endif
            endif
            if (NE3adopt(iMT) == 0) then
              do nin = 1, numcut
                if (eninc(nin) == EminMeV) cycle
                if (eninc(nin) <= Ethexcliso(idc, nex)) cycle
                Eev = eninc(nin) * 1.e6
                if (Eev < Ehigh) cycle
                if (xsexcliso(idc, nex, nin) == 0.) cycle
                iE = iE + 1
                E3(MT, iE) = eninc(nin) * 1.e6
                xs(MT, iE) = xsexcliso(idc, nex, nin) * 1.e-3
              enddo
            else
              do i = 1, NE3adopt(iMT)
                E3(MT, i) = E3adopt(iMT, i)
                ee = E3(MT, i) * 1.e-6
                call locate(eninc, 0, numcut, ee, nin)
                xs(MT, i) = xs3adopt(iMT, i) * branchiso(idc, nex, nin)
              enddo
              iE = NE3adopt(iMT)
            endif
            if (iE /= 1) then
!
! ENDF-6 parameters
!
              LFS3(MT) = iso
              NP(MF, MT) = iE
              NR(MF, MT) = 1
              NBT(MF, MT, 1) = iE
              INTER(MF, MT, 1) = 2
              LR3(MT) = 0
              mtexist(MF, MT) = .true.
            endif
          endif
!
! Make EAF-string for first line of MT-section
!
  100         eafstring(MT)(1:1) = nuclid(1:1)
          if (nuclid(2:2) /= ' ') eafstring(MT)(2:2) = achar(iachar(nuclid(2:2)) - 32)
          if (Lisomer == 0) then
            eafstring(MT)(3:3) = '-'
          else
            eafstring(MT)(3:3) = '*'
          endif
          write(eafstring(MT)(4:6), '(i3)') Atarget
          write(eafstring(MT)(7:8), '(a1, ",")') achar(iachar(parsym(k0)) - 32)
          str = '( ,'
          write(str(2:2), '(a1)') parsym(k0)
          ibeg = 0
          iend = 0
          do i = 1, 64
            if (mtstring(MT)(i:i + 2) == str) ibeg = i + 3
            if (mtstring(MT)(i:i) == ')') iend = i - 1
          enddo
          ilength = iend - ibeg + 1
          if (ibeg > 0 .and. iend > 0) then
            write(eafstring(MT)(9:9+ilength), '(a)') mtstring(MT)(ibeg:iend)
            if (MT > 300) write(eafstring(MT)(9+ilength:9+ilength), '("*")')
          endif
          do i = 9, 20
            if (eafstring(MT)(i:i) >= 'a' .and. eafstring(MT)(i:i) <= 'z') eafstring(MT)(i:i) = &
              achar(iachar(eafstring(MT)(i:i)) - 32)
          enddo
          write(eafstring(MT)(26:35), '("TALYS-2.2 ")')
        enddo
      endif
    enddo
  enddo
  return
end subroutine make3eaf
! Copyright A.J. Koning 2021
