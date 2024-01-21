subroutine endfinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of ENDF-6 formats
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
! All global variables
!   nummt          ! number of MT numbers
! Variables for input of ENDF library type
!   endffile       ! name of ENDF file
!   flageaf        ! flag for EAF - formatted activation library
!   flaggpf        ! flag for general purpose library
! Variables for input of ENDF structure
!   flagexclude    ! flag to specifically exclude an MT number (per MT)
!   flaginclude    ! flag to specifically include an MT number (per MT)
!   flagmtall      ! flag to include all defined MT numbers from ENDF manual
!   flagmtextra    ! flag to include extra MT numbers up to MT200
!   flagmulti      ! flag to include multi - chance fission
! Variables for input of ENDF MF1
!   identifier     ! library identifier
! Variables for info from TALYS
!   Atarget        ! mass number of nucleus
!   k0             ! index of incident particle
!   Lisomer        ! isomeric number of target
!   Ltarget        ! excited level of target
!   Ztarget        ! charge number of nucleus
! Constants
!   light          ! mass number of lightest stable isotope
!   nuc            ! nuclide
!   parname        ! name of particle
!   parsym         ! symbol of particle
! Variables for reaction initialization
!   enincmax      ! maximal incident energy
!   isochar       ! symbol for isomer
!   massN         ! mass of nucleus in neutron units
!   nuclid        ! nuclide
! Variables for info from TALYS
!   eninc             ! incident energy
! Variables for initialization of ENDF format
!   AWR           ! standard mass parameter
!   blank1        ! blank string
!   blank2        ! blank string
!   CONT          ! ENDF - 6 format
!   FEND          ! ENDF - 6 format
!   HEAD          ! ENDF - 6 format
!   iclean        ! number of cleaned points
!   INT3          ! blank string
!   LIS           ! state number of target nucleus
!   LISO          ! isomeric state number
!   MAT           ! MAT number
!   MEND          ! ENDF - 6 format
!   MTid          ! channel identifier for MT - number
!   MTinel        ! MT - number for inelastic scattering
!   MTmax         ! highest MT number for exclusive channels
!   MTnum         ! channel identifier for MT - number (ENDF format)
!   rec0          ! first line of ENDF-6 file
!   SEND          ! ENDF - 6 format
!   TEND          ! ENDF - 6 format
!   TEXT          ! ENDF - 6 format
!   VALUE         ! ENDF - 6 format
!   ZA            ! standard charge parameter
! Variables for MF1
!   EMAX        ! upper limit of energy range for evaluation
!   ERN         ! total energy minus energy of neutrinos
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: Astring      ! string for mass number
  character(len=4)  :: ext          ! filename extension
  character(len=9)  :: energystring ! energy string
  character(len=16) :: filetype     ! type of ENDF-6 file
  integer           :: isp(4)       ! help variable for Es-Lr MAT numbers
  integer           :: MT           ! MT-number
  isp = (/9920, 9945, 9965, 9980/)
!
! ********************** Assignment of MT-numbers **********************
!
! 1. Initialization.
!
  iclean = 0
  MTnum = -1
  MTid = -1
  mtstring = ''
  eafstring = ''
  mfexist = .false.
  mtexist = .false.
  reacstring = ''
!
! Each MT-number corresponds to a certain number of particles in the outgoing channel.
!
! 1. Official ENDF-6 MT numbers.
!
  MTnum(4) = 100000
  MTnum(11) = 201000
  MTnum(16) = 200000
  MTnum(17) = 300000
  MTnum(18) = -99
  MTnum(19) = -98
  MTnum(20) = -97
  MTnum(21) = -96
  MTnum(22) = 100001
  MTnum(23) = 100003
  MTnum(24) = 200001
  MTnum(25) = 300001
  MTnum(28) = 110000
  MTnum(29) = 100002
  MTnum(30) = 200002
  MTnum(32) = 101000
  MTnum(33) = 100100
  MTnum(34) = 100010
  MTnum(35) = 101002
  MTnum(36) = 100102
  MTnum(37) = 400000
  MTnum(38) = -95
  MTnum(41) = 210000
  MTnum(42) = 310000
  MTnum(44) = 120000
  MTnum(45) = 110001
  MTnum(91) = 100000
  MTnum(102) = 000000
  MTnum(103) = 010000
  MTnum(104) = 001000
  MTnum(105) = 000100
  MTnum(106) = 000010
  MTnum(107) = 000001
  MTnum(108) = 000002
  MTnum(109) = 000003
  MTnum(111) = 020000
  MTnum(112) = 010001
  MTnum(113) = 000102
  MTnum(114) = 001002
  MTnum(115) = 011000
  MTnum(116) = 010100
  MTnum(117) = 001001
  MTnum(649) = 010000
  MTnum(699) = 001000
  MTnum(749) = 000100
  MTnum(799) = 000010
  MTnum(849) = 000001
  if (flagmtall) then
    MTid = MTnum
  else
    MTid(4) = MTnum(4)
    MTid(16) = MTnum(16)
    MTid(17) = MTnum(17)
    if (k0 <= 1) then
      MTid(18) = MTnum(18)
      if (flagmulti) then
        MTid(19) = MTnum(19)
        MTid(20) = MTnum(20)
        MTid(21) = MTnum(21)
        MTid(38) = MTnum(38)
      endif
    endif
    MTid(22) = MTnum(22)
    MTid(28) = MTnum(28)
    MTid(37) = MTnum(37)
    MTid(91) = MTnum(91)
    MTid(102) = MTnum(102)
    MTid(103) = MTnum(103)
    MTid(104) = MTnum(104)
    MTid(105) = MTnum(105)
    MTid(106) = MTnum(106)
    MTid(107) = MTnum(107)
    MTid(649) = MTnum(649)
    MTid(699) = MTnum(699)
    MTid(749) = MTnum(749)
    MTid(799) = MTnum(799)
    MTid(849) = MTnum(849)
  endif
  do MT = 1, nummt
    if (flaginclude(MT)) MTid(MT) = MTnum(MT)
    if (flagexclude(MT)) MTid(MT) = -1
  enddo
  if (MTid(103) ==  -1) MTid(649) = -1
  if (MTid(104) ==  -1) MTid(699) = -1
  if (MTid(105) ==  -1) MTid(749) = -1
  if (MTid(106) ==  -1) MTid(799) = -1
  if (MTid(107) ==  -1) MTid(849) = -1
!
! 2. Special high-energy MT numbers
!
  if (flagmtextra) then
    MTnum(152) = 500000
    MTnum(153) = 600000
    MTnum(154) = 200100
    MTnum(155) = 000101
    MTnum(156) = 410000
    MTnum(157) = 301000
    MTnum(158) = 101001
    MTnum(159) = 210001
    MTnum(160) = 700000
    MTnum(161) = 800000
    MTnum(162) = 510000
    MTnum(163) = 610000
    MTnum(164) = 710000
    MTnum(165) = 400001
    MTnum(166) = 500001
    MTnum(167) = 600001
    MTnum(168) = 700001
    MTnum(169) = 401000
    MTnum(170) = 501000
    MTnum(171) = 601000
    MTnum(172) = 300100
    MTnum(173) = 400100
    MTnum(174) = 500100
    MTnum(175) = 600100
    MTnum(176) = 200010
    MTnum(177) = 300010
    MTnum(178) = 400010
    MTnum(179) = 320000
    MTnum(180) = 300002
    MTnum(181) = 310001
    MTnum(182) = 001100
    MTnum(183) = 111000
    MTnum(184) = 110100
    MTnum(185) = 101100
    MTnum(186) = 110010
    MTnum(187) = 101010
    MTnum(188) = 100110
    MTnum(189) = 100101
    MTnum(190) = 220000
    MTnum(191) = 010010
    MTnum(192) = 001010
    MTnum(193) = 000011
    MTnum(194) = 420000
    MTnum(195) = 400002
    MTnum(196) = 410001
    MTnum(197) = 030000
    MTnum(198) = 130000
    MTnum(199) = 320001
    MTnum(200) = 520000
    do MT = 1, 200
      MTid(MT) = MTnum(MT)
    enddo
    MTmax = 200
  else
    MTmax = 117
  endif
!
! MT number of inelastic scattering
!
  if (k0 == 0) MTinel = 102
  if (k0 == 1) MTinel = 4
  if (k0 == 2) MTinel = 103
  if (k0 == 3) MTinel = 104
  if (k0 == 4) MTinel = 105
  if (k0 == 5) MTinel = 106
  if (k0 == 6) MTinel = 107
!
! ENDF format variables
!
  write(blank1(1:11), '(11x)')
  write(blank2(1:66), '(66x)')
  CONT = '(2a11, 4i11, i4, i2, i3, i5)'
  FEND = '(a66, i4, i2, i3, i5)'
  HEAD = '(2a11, 4i11, i4, i2, i3, i5)'
  INT3 = '(6i11, t67, i4, i2, i3, i5)'
  MEND = '(a66, i4, i2, i3, i5)'
  SEND = '(a66, i4, i2, i3, i5)'
  TEND = '(a66, i4, i2, i3, i5)'
  TEXT = '(a66, i4, i2, i3, i5)'
  VALUE = '(6a11, i4, i2, i3, i5)'
!
! *********************** Initialization of target *********************
!
  LIS = Ltarget
  LISO = Lisomer
  isochar = " "
  if (Lisomer == 1) isochar = "m"
  if (Lisomer == 2) isochar = "n"
  if (Lisomer == 3) isochar = "o"
!
! *********************** Mass and energy variables ********************
!
  ZA = real(1000 * Ztarget + Atarget)
  AWR = massN
  EMAX = enincmax * 1.e6
  ERN = 2.e8
!
! ***************************** MAT numbers ****************************
!
! In the EAF-format, the MAT numbers are different from the ENDF-6 format.
! (Not all possibilities for MAT-numbers are yet covered, e.g. Einsteinium and higher, to be done).
!
  if (flageaf .or. Ztarget > 103) then
    MAT = mod(Ztarget * 100, 10000) + mod(Atarget, 100)
    if (flageaf) then
      if (Lisomer == 1) MAT = MAT + 50
      if (Lisomer == 2) MAT = MAT + 70
      if (Lisomer == 3) MAT = MAT + 90
    endif
  else
    if (Ztarget < 99) MAT = Ztarget * 100 + 25 + 3 * (Atarget - light(Ztarget)) + Lisomer
    if (Ztarget == 99) MAT = Ztarget * 100 + Atarget - light(Ztarget) + 1 + 5 * Lisomer
    if (Ztarget > 99) MAT = isp(Ztarget - 99) + Atarget - light(Ztarget) + 1 + 5 * Lisomer
  endif
!
! *********************** Name of ENDF file ****************************
!
  if (endffile(1:1) == ' ') then
    Astring = '000'
    write(Astring(1:3), '(i3.3)') Atarget
    endffile = parsym(k0)//'-'//trim(nuc(Ztarget))//Astring//isochar
    if (flaggpf) then
      ext = '.gpf'
    else
      if (flageaf) then
        ext = '.eaf'
      else
        ext = '.acf'
      endif
    endif
    endffile = trim(endffile) //ext
  endif
!
! Write header of ENDF-6 file
!
  if (flaggpf) then
    filetype = 'gen. purp. file:'
  else
    filetype = 'activation file:'
  endif
  if (k0 == 1) then
    energystring = ' 1.e-5 eV'
  else
    energystring = '      MeV'
    write(energystring(2:5), '(i4)') int(eninc(1))
  endif
  write(rec0(1:37), '(" ", a10, " ", a16, " ", a8)') identifier, filetype, parname(k0)
  write(rec0(38:66), '(" + ", a2, "-", i3, a1, a9, " -", i4, " MeV")') nuclid, Atarget, isochar, energystring, int(enincmax)
  write(rec0(67:80), '("  99 0  0    0")')
  return
end subroutine endfinitial
! Copyright A.J. Koning 2021
