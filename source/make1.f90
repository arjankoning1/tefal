subroutine make1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF1
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   nummf          ! number of MF numbers
!   nummt          ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt          ! logical for existence of MF information (per MT)
! Variables for input of ENDF structure
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   nomf           ! flag to exclude an entire MF
! Variables for input of ENDF MF1
!   author         ! author
!   diffweight     ! differential weight to be put in MF1
!   endftext       ! file with MF1 information
!   identifier     ! library identifier
!   lab            ! laboratory
! Variables for reaction initialization
!   includeres     ! flag to include resonance parameters
!   isochar        ! symbol for isomer
!   nuclid         ! nuclide
!   relmass        ! mass relative to neutron mass
! Variables for initialization of ENDF format
!   MAT            ! MAT number
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
! Variables for path names
!   month          ! month
!   year           ! year
! Variables for info from TALYS
!   Atarget        ! mass number of nucleus
!   k0             ! index of incident particle
!   targetE        ! energy of target
!   Ztarget        ! charge number of nucleus
! Constants
!   parA           ! mass number of particle
!   parname        ! name of particle
!   parZ           ! charge number of particle
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF1
!   ALAB           ! Lab
!   AUTH           ! author(s)
!   AWI            ! mass of projectile in neutron units
!   Cnubar         ! coefficients of the polynomial expansion
!   DDATE          ! date of distribution
!   E1             ! incident energy for MF1 (in ENDF - 6 format)
!   EDATE          ! date of evaluationon
!   ELIS           ! excitation energy of target nucleus
!   EMAX           ! upper limit of energy range for evaluation
!   ENDATE         ! date string
!   HSUB1          ! string
!   HSUB2          ! string
!   HSUB3          ! string
!   LDRV           ! special evaluation flag
!   LFI            ! flag for fission
!   LNU            ! polynomial(LNU = 1) or tabulated (LNU = 2) representation
!   LREL           ! release number
!   LRP            ! flag that indicates presence of file 2 (resonances)
!   MODN           ! modification indicator
!   NC             ! number of lines for MF / MT number
!   NCnubar        ! number of terms in polynomial expansion
!   NFOR           ! library format
!   NLIB           ! library identifier
!   NMOD           ! modification number of evaluation
!   NSUB           ! sub - library number
!   ntxt           ! number of text lines
!   nubar          ! average number of neutrons per fission
!   NVER           ! library version number
!   NWD            ! number of text lines
!   NXC            ! number of directory lines
!   RDATE          ! date of revision
!   REF            ! reference
!   STA            ! target stability flag
!   TEMP           ! target temperature
!   txt            ! text line
!   ZSYMAM         ! string
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical          :: flagfismf1     ! flag for fission
  logical          :: lexist         ! logical to determine existence
  character(len=1) :: isochar_up     !
  character(len=4) :: mffile         ! MF filename
  integer          :: i              ! counter
  integer          :: imf            ! MF counter
  integer          :: imt            ! MT counter
  integer          :: IPART          ! 1000*Z+A
  integer          :: istat          ! error code
  integer          :: ITYPE          ! type of data
  integer          :: j              ! counter
  integer          :: maxE           ! maximum energy point
  integer          :: MF             ! MF-number
  integer          :: MT             ! MT-number
  integer          :: nutype         ! help variable
  real(sgl)        :: factor         ! multiplication factor
  real(sgl)        :: nudel          ! delayed nubar
!
! ****************** Read MF1 information from file ********************
!
  MF = 1
  ntxt = 1
  inquire (file = endftext, exist = lexist)
  if (lexist) then
    open (unit = 1, file = endftext, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(endftext, istat)
    do
      read(1, '(a66)', iostat = istat) txt(ntxt)
      if (istat == -1) exit
      ntxt = ntxt + 1
    enddo
    close (unit = 1)
    ntxt = ntxt - 1
    if (diffweight >= 0.) then
      do i = 1, ntxt
        if (txt(i)(1:20) == 'Differential weight:') then
          write(txt(i)(22:33), '(es12.5)') diffweight
          exit
        endif
      enddo
    endif
  endif
!
! **************************** Make MF1 ********************************
!
  if (k0 == 1) then
    if (includeres) then
      LRP = 1
    else
      LRP = 0
    endif
  else
    LRP = -1
  endif
  flagfismf1 = (k0 <= 1 .and. flagfission .and. mtexist(1, 452) .and. .not. flagfis10)
  if (flagfismf1) then
    LFI = 1
  else
    LFI = 0
  endif
  NLIB = 17
  NMOD = 1
  ELIS = targetE * 1.e6
  if (targetE == 0.) then
    STA = 0.
  else
    STA = 1.
  endif
  NFOR = 6
!
! Overwrite AWI with slightly less precise values to keep BNL checking software happy.
!
  AWI = real(relmass(k0))
! if (k0 == 2) AWI = 0.998623
! if (k0 == 3) AWI = 1.996256
! if (k0 == 4) AWI = 2.989596
! if (k0 == 5) AWI = 2.989033
! if (k0 == 6) AWI = 3.967131
  LREL = 1
  IPART = 1000 * parZ(k0) + parA(k0)
  ITYPE = 0
  NSUB = 10 * IPART + ITYPE
  NVER = 2023
!
! The number of directory lines NXC is calculated at the end of the subroutine.
!
  TEMP = 0.
  LDRV = 0
  NWD = 5 + ntxt
  isochar_up = ' '
  if (isochar == 'm') isochar_up = 'M'
  if (isochar == 'n') isochar_up = 'N'
  if (isochar == 'o') isochar_up = 'O'
  write(ZSYMAM, '(i3, "-", a2, "-", i3, a1)') Ztarget, nuclid, Atarget, isochar_up
  ALAB = lab
  EDATE = 'EVAL-     '
  write(EDATE(6:8), '(a3)') month
  write(EDATE(9:10), '(i2.2)') year
  AUTH = author
  REF = identifier
  DDATE = 'DIST-     '
  RDATE = 'REV1-     '
  ENDATE = '      '
  write(HSUB1, '("----", a10, t23, "Material ", i4, t45, "REVISION 1")') identifier, MAT
  write(HSUB2, '("-----Incident ", a8, " data", 39x)') parname(k0)
  write(HSUB3, '("------ENDF-6 Format", 47x)')
!
! ******************* Make or adjust fission information ***************
!
! Unify energy grid in case of total nubar calculated by TAFIS
!
  if (flagfismf1) then
    if (adopt(MF, 452)) then
      NP(MF, 456) = NP(MF, 452)
      NR(MF, 456) = NR(MF, 452)
      NC(MF, 456) = NC(MF, 452)
      LNU(456) = LNU(452)
      do i = 1, NR(MF, 452)
        NBT(MF, 456, i) = NBT(MF, 452, i)
        INTER(MF, 456, i) = INTER(MF, 452, i)
      enddo
      do i = 1, NP(MF, 452)
        E1(3, i) = E1(1, i)
        do j = 1, NP(MF, 455) - 1
          if (E1(3, i) >= E1(2, j) .and. E1(3, i) <= E1(2, j + 1)) then
            factor = (E1(3, i) - E1(2, j)) / (E1(2, j + 1) - E1(2, j))
            nudel = nubar(2, j) + factor * (nubar(2, j + 1) - nubar(2, j))
            nubar(3, i) = nubar(1, i) - nudel
          endif
        enddo
      enddo
    endif
!
! Deal with photons and possible extension to high energies
!
    do nutype = 1, 3
      if (nutype == 1) MT = 452
      if (nutype == 2) MT = 455
      if (nutype == 3) MT = 456
      if (k0 == 0) then
        do i = 1, NP(MF, MT)
          E1(nutype, i) = E1(nutype, i) - QM(4)
        enddo
      endif
      maxE = NP(MF, MT)
      if (maxE == 0) cycle
      do i = 1, maxE
        if (EMAX < E1(nutype, i)) then
          NP(MF, MT) = i - 1
          NBT(MF, MT, NR(MF, MT)) = i - 1
          exit
        endif
      enddo
      if (EMAX > E1(nutype, NP(MF, MT)) + 10.) then
        if (LNU(MT) == 1) then
          if (NCNUBAR == 2) then
            NCnubar = 3
            Cnubar(3) = (9.9 - Cnubar(1) - Cnubar(2) * 2.e8) / (2.e8 **2)
          endif
        else
          NP(MF, MT) = NP(MF, MT) + 1
          if (mod(NP(MF, MT), 3) == 1) NC(MF, MT) = NC(MF, MT) + 1
          NBT(MF, MT, NR(MF, MT)) = NBT(MF, MT, NR(MF, MT)) + 1
          E1(nutype, NP(MF, MT)) = EMAX
          nubar(nutype, NP(MF, MT)) = nubar(nutype, NP(MF, MT) - 1)
        endif
      endif
    enddo
  endif
!
! **** Create the directory by reading all the individual MF-files *****
!
! Assess number of directory lines by checking the existence of each MF/MT number.
!
  NXC = 0
  mfexist(1) = .true.
  mtexist(1, 451) = .true.
  do imf = 1, nummf
    do imt = 1, nummt
      if (mtexist(imf, imt)) then
        NXC = NXC + 1
        MODN(imf, imt) = 1
      endif
    enddo
  enddo
  NC(1, 451) = 4 + NXC + NWD
!
! Read number of lines for each MF/MT section
!
  do imf = 2, nummf
    if (nomf(imf)) then
      do imt = 1, nummt
        mtexist(imf, imt) = .false.
      enddo
      cycle
    endif
    mffile = 'MF  '
    if (imf < 10) then
      write(mffile(3:3), '(i1)') imf
    else
      write(mffile(3:4), '(i2)') imf
    endif
    inquire (file = mffile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = mffile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(mffile, istat)
      imt = 0
      do
        read(1, '(72x, i3)', iostat = istat) MT
        if (istat == -1) exit
        if (MT /= 0) then
          imt = imt + 1
          NC(imf, MT) = imt
        else
          imt = 0
        endif
      enddo
      close (unit = 1)
    endif
  enddo
  return
end subroutine make1
! Copyright A.J. Koning 2021
