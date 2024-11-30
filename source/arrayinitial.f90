subroutine arrayinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of arrays
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
! ************ Initialization ******************************************
!
!
  xsnonel = 0.
  xselas = 0.
  xsgamexcl = 0.
  xsnonth = 0.
  xsgamdisctot = 0.
  xsgamtot = 0.
  xsgamcont = 0.
  yieldgam = 0.
  yieldtot = 0.
  yieldcont = 0.
  xsgam = 0.
  yielddisc = 0.
  xsnonel6 = 0.
  xselas6 = 0.
  xstot6 = 0.
  E1 = 0.
  nubar = 0.
  Cnubar = 0.
  isolevel = 0
  rpisolevel = 0
  Qexcliso = 0.
  Ethexcliso = 0.
  xsexcliso = 0.
  branchiso = 0.
  Estartdis = 0.
  Egammadis = 0.
  xsgamdis = 0.
  Ethexcl = 0.
  idchannel = 0
  xsexcl = 0.
  Qexcl = 0.
  Nisomer = 0
  Ngam = 0
  Eout = 0.
  specexcl = 0.
  Ehist = 0.
  f0ex = 0.
  Eocum = 0.
  preeqratio = 0.
  buratio = 0.
  xsemis = 0.
  Ehistcum = 0.
  f0cum = 0.
  Erec = 0.
  recexcl = 0.
  Ehistrec = 0.
  f0exrec = 0.
  Erecrp = 0.
  recrp = 0.
  Ehistcumrec = 0.
  f0cumrec = 0.
  edisc = 0.
  xsbin = 0.
  xsintdisc = 0.
  nbeg = 0
  nend = 0
  nbegrec = 0
  nendrec = 0
  nout = 0
  noutrec = 0
  Eparticles = 0.
  Erecav = 0.
  nbegcumrec = 0
  nendcumrec = 0
  xscont = 0.
  xsngn = 0.
  Nisorp = 0
  nucmass = 0.
  Qrp = 0.
  Ethrp = 0.
  Erpiso = 0.
  Qrpiso = 0.
  Ethrpiso = 0.
  xsngnsum = 0.
  noutrecrp = 0
  xsrp = 0.
  Yrp = 0.
  xsrpiso = 0.
  Yrpiso = 0.
  Eurr = 0.
  NLSurr = 0
  NJSurr = 0
  Jurr = 0.
  Djlurr = 0.
  GNjlurr = 0.
  GGjlurr = 0.
  GXjlurr = 0.
  GFjlurr = 0.
  Nurr = 0
  Ngamdis = 0
  Nbranch = 0
  yieldg = 0.
  levexist = .true.
  Qdisc = 0.
  jdis = 0.
  Ethdisc = 0.
  Egamdis = 0.
  yieldratio = 0.
  xsprod = 0.
  yieldp = 0.
  yieldany = 0.
  yieldngn = 0.
  ncumout = 0
  nbegcum = 0
  nendcum = 0
  ncleg = 0
  cleg0 = 0.
  cpang = 0.
  fcpang = 0.
  elasni = 0.
  xsdisc = 0.
  Egamma = 0.
  Estart = 0.
  rmu = 0.
  fang = 0.
  xsang = 0.
  numendisc = 0
  ndisc = 0
  NSt = 0
  Nisocov = 0
  NE3res = 0
  NE3adopt = 0
  xs = 0.
  E3 = 0.
  leg = 0.
  E3res = 0.
  xs3res = 0.
  E3adopt = 0.
  xs3adopt = 0.
  MTindex = 0
  MTindexiso = -1
  MTintindex = 0
  MTintindexiso = -1
  Rmt = 0.
  Ecov = 0.
  Eleg = 0.
  Ecov = 0.
  Eleg = 0.
  relerr = 0.
  xserr = 0.
  Rcov = 0.
  Rleg = 0.
  LI4 = 0
  LCT = 0
  LVT = 0
  LTT = 0
  NRh = 0
  NEh = 0
  NE = 0
  LR3 = 0
  NE3res = 0
  NE3adopt = 0
  LFS3 = 0
  QM = 0.
  QI = 0.
  EthMT = 0.
  MODN = 0
  NC = 0
  NP = 0
  NR = 0
  NK = 0
  NBT = 0
  INTER = 0
  NL = 0
  INTERh = 0
  XMF1 = 0
  XLFS1 = 0
  LF = 0
  NR5 = 0
  NP5 = 0
  NR5e = 0
  NE5e = 0
  U = 0.
  NE6ea = 0
  NBT5 = 0
  INTER5 = 0
  NBT5e = 0
  INTER5e = 0
  NBT5e2 = 0
  INTER5e2 = 0
  NR5e2 = 0
  NF = 0
  pE = 0.
  E5 = 0.
  E5p = 0.
  TM5 = 0.
  gE5 = 0.
  b6 = 0.
  kpart = 0
  b6gam = 0.
  flagrec = .false.
  ErpZAiso = 0.
  QMZA = 0.
  QIZA = 0.
  IZAP = 0
  LMF = 0
  LFS = 0
  b6rec = 0.
  Ey = 0.
  LRU = 0
  LRF = 0
  NRO = 0
  NAPS = 0
  NLS = 0
  LAD = 0
  NLSC = 0
  LSSF = 0
  EL = 0.
  EH = 0.
  SPI = 0.
  AP = 0.
  Lres = 0
  NRS = 0
  NJS = 0
  AWRI = 0.
  QX = 0.
  APL = 0.
  LRX = 0
  Er = 0.
  AJ = 0.
  GT = 0.
  GN = 0.
  GG = 0.
  GF = 0.
  GFA = 0.
  GFB = 0.
  Es = 0.
  D = 0.
  GX = 0.
  GN0 = 0.
  GGu = 0.
  GFu = 0.
  kINT = 0
  NEu = 0
  AJU = 0.
  AMUX = 0.
  AMUN = 0.
  AMUG = 0.
  AMUF = 0.
  E2 = 0.
  APE = 0.
  return
end subroutine arrayinitial
! Copyright A.J. Koning 2021
