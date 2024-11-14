module A0_tefal_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : General module with all global variables
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2024-11-14: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Definition of single and double precision variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
!
!-----------------------------------------------------------------------------------------------------------------------------------
! All global dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: memorypar=6                     ! memory parameter
  integer, parameter :: numpar=6                        ! number of particles
  integer, parameter :: numelem=124                     ! number of elements
  integer, parameter :: numenin=600                     ! number of incident energies
  integer, parameter :: numen4=4000                     ! number of incident energies for MF4
  integer, parameter :: numlines=1000                   ! number of input lines
  integer, parameter :: numtxt=1000                     ! number of text lines
  integer, parameter :: nummt=850                       ! number of MT numbers
  integer, parameter :: numang=100                      ! number of angles
  integer, parameter :: nummf=40                        ! number of MF numbers
  integer, parameter :: numchan=200                     ! maximum number of exclusive channels
  integer, parameter :: numlevels=40                    ! maximum number of discrete levels
  integer, parameter :: numlevin=110                    ! number of discrete levels
  integer, parameter :: numZ=10                         ! maximal number of protons from initial compound nucleus
  integer, parameter :: numN=20                         ! maximal number of neutrons from initial compound nucleus
  integer, parameter :: numen6=memorypar*2000           ! number of incident energies
  integer, parameter :: numen2=175+16*(memorypar-2)     ! number of emission energies
  integer, parameter :: numl=60                         ! number of l values
  integer, parameter :: numgam=100                      ! number of gamma lines
  integer, parameter :: numres=5                        ! number of resonance sections
  integer, parameter :: numnrs=6000                     ! number of resonances
  integer, parameter :: numrespar=5                     ! number of different resonance parameters
  integer, parameter :: numrescov=1000                  ! number of resonances with covariances
  integer, parameter :: numlres=5                       ! number of l-values
  integer, parameter :: numjres=10                      ! number of j-values
  integer, parameter :: numPP=10                        ! number of particle pairs
  integer, parameter :: numch7=10                       ! number of channels in LRF7 representation
  integer, parameter :: nummtres=107                    ! number of MT numbers with resonances
  integer, parameter :: numint=10                       ! number of interpolation sections
  integer, parameter :: numsec=300                      ! number of sections
  integer, parameter :: numsecea=8                      ! number of sections for energy-angle distributions
  integer, parameter :: numsecg=1                       ! number of sections for gamma cross sections
  integer, parameter :: numchancov=14                   ! number of channels with inter-channel correlations
  integer, parameter :: numencov=200                    ! number of incident energies for covariances
  integer, parameter :: numencovtot=1+numencov*numencov ! number of energies for covariances
  integer, parameter :: numenspec=80                    ! number of incident energies for spectra
  integer, parameter :: numenrec=25                     ! number of incident energies for recoils
  integer, parameter :: numengam=80                     ! number of incident energies for gamma cross sections
  integer, parameter :: numenang=80                     ! number of incident energies for angular distributions
  integer, parameter :: numddx=18                       ! number of DDX angles
  integer, parameter :: numrp=(numZ+1)*(numN+1)         ! nunber of residual products
  integer, parameter :: numiso=9                        ! number of isomers
  integer, parameter :: numleg=10                       ! number of Legendre coefficients
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for path names
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=3)   :: month   ! month
  character(len=132) :: path    ! directory containing files to be read
  integer            :: year    ! year
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Constants
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(-1:numpar) :: parsym    ! symbol of particle
  character(len=2), dimension(numelem)   :: nuc       ! nuclide
  character(len=8), dimension(-1:numpar) :: parname   ! name of particle
  integer                                :: fislim    ! mass above which nuclide fissions
  integer, dimension(numelem)            :: light     ! mass number of lightest stable isotope
  integer, dimension(0:numpar)           :: parA      ! mass number of particle
  integer, dimension(0:numpar)           :: parN      ! neutron number of particle
  integer, dimension(0:numpar)           :: parZ      ! charge number of particle
  real(dbl), dimension(0:numpar)         :: parmass   ! mass of particle in a.m.u.
  real(sgl), dimension(0:numpar)         :: parspin   ! spin of particle
  real(sgl)                              :: pi        ! pi
  real(sgl)                              :: xsepshigh ! upper limit for cross sections in millibarns
  real(sgl)                              :: xsepslow  ! lower limit for cross sections in millibarns
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for info from TALYS
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                         :: flagelastic  ! flag to give priority to elastic cross section for normalization
  logical                         :: flagtalysdet ! flag for detailed ENDF-6 information from TALYS
  logical                         :: flagtalysrec ! flag for recoil information from TALYS
  logical                         :: flagtalysurr ! flag for URR information from TALYS
  logical                         :: flagblock    ! flag to block spectra, angle and gamma files
  integer                         :: Ainit        ! mass number of initial compound nucleus
  integer                         :: Atarget      ! mass number of nucleus
  integer                         :: k0           ! index of incident particle
  integer                         :: Lisomer      ! isomeric number of target
  integer                         :: Ltarget      ! excited level of target
  integer                         :: NLmax        ! maximum number of included discrete levels
  integer                         :: numinc       ! number of incident energies
  integer                         :: Zinit        ! charge number of initial compound nucleus
  integer                         :: Ztarget      ! charge number of nucleus
  real(sgl), dimension(0:numenin) :: eninc        ! incident energy
  real(sgl)                       :: targetE      ! energy of target
  real(dbl)                       :: tarmass      ! mass of nucleus
  real(sgl)                       :: targetspin   ! spin of target
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for reaction initialization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                        :: includeres ! flag to include resonance parameters
  character(len=1)               :: isochar    ! symbol for isomer
  character(len=2)               :: nuclid     ! nuclide
  integer, dimension(numenang)   :: Eangindex  ! enegy index for angular distribution
  integer, dimension(numengam)   :: Egamindex  ! enegy index for gamma cross sections
  integer, dimension(numenspec)  :: Especindex ! enegy index for spectra
  integer                        :: Nenang     ! number of incident energies for angular distributions
  integer                        :: Nengam     ! number of incident energies for gamma production
  integer                        :: Nenspec    ! number of incident energies for spectra
  integer                        :: nlevmax    ! number of included discrete levels
  integer                        :: numcut     ! number of energies before high-energy format
  integer                        :: numcut4    ! number of energies before MF4 high-energy format
  real(sgl)                      :: EHres      ! upper energy in resonance range
  real(sgl)                      :: eninccut   ! last incident energy before high-energy format
  real(sgl)                      :: enincmax   ! maximal incident energy
  real(sgl)                      :: EmineV     ! minimal incident energy in eV
  real(sgl)                      :: EminMeV    ! minimal incident energy in MeV
  real(sgl)                      :: massN      ! mass of nucleus in neutron units
  real(dbl), dimension(0:numpar) :: relmass    ! mass relative to neutron mass
  real(sgl), dimension(0:numang) :: rmu        ! cosine of the angle
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for reading TEFAL input lines
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numlines) :: inline ! input line
  integer                                 :: nlines ! number of input lines
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for input of ENDF structure
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                   :: flaglinenum ! flag to write line numbers of ENDF file
  logical                   :: flagaddlow  ! flag to add low-energy Q>0 reactions to total cross section
  logical                   :: flagcapt6   ! flag to put MT102 gamma prod. in MF6 instead of MF12/14/15
  logical                   :: flagdisc6   ! flag for discrete angular distribution and gamma prod. in MF6
  logical, dimension(nummt) :: flagexclude ! flag to specifically exclude an MT number (per MT)
  logical                   :: flagfis10   ! flag to put (subactinide) fission cross sections in MF10
  logical                   :: flaggam13   ! flag to use MF13 for gamma prod. instead of MF12 (if not in MF6)
  logical                   :: flaggamdisc ! flag to store gamma prod. per level (MT51..) instead of in MF4
  logical                   :: flaggamspec ! flag to store gamma prod. only as a spectrum and not per level
  logical, dimension(nummt) :: flaginclude ! flag to specifically include an MT number (per MT)
  logical                   :: flagmtall   ! flag to include all defined MT numbers from ENDF manual
  logical                   :: flagmtextra ! flag to include extra MT numbers up to MT200
  logical                   :: flagmulti   ! flag to include multi-chance fission
  logical                   :: flagngn     ! flag to include (n,gamma n) data
  logical                   :: flagpara    ! flag to include partial cross sections for alphas
  logical                   :: flagparn    ! flag to include partial cross sections for neutrons
  logical                   :: flagpard    ! flag to include partial cross sections for deuterons
  logical                   :: flagparh    ! flag to include partial cross sections for helions
  logical                   :: flagparp    ! flag to include partial cross sections for protons
  logical                   :: flagpart    ! flag to include partial cross sections for tritons
  logical                   :: flagpart6   ! flag for gam. prod. for partial c.s. in MF6 not in MF12/14/15
  logical                   :: flagrenorm  ! flag for renormalization of spectra
  logical                   :: flagrp10    ! flag to put residual production cross sections in MF10
  logical                   :: flagrp6     ! flag to put residual production cross sections in MF6
  logical                   :: flagrecoil  ! flag to include recoil information
  logical                   :: flagres     ! flag to include resonance parameters
  logical                   :: flagsubfis  ! flag to include subactinide fission
  logical                   :: flagtabddx  ! flag to give explicit DDX in MF6
  logical, dimension(nummf) :: nomf        ! flag to exclude an entire MF
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for input of ENDF library type
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical            :: flageaf     ! flag for EAF-formatted activation library
  logical            :: flagendfdet ! flag for detailed ENDF-6 information per channel
  logical            :: flagclean   ! flag to clean up double points
  logical            :: flaggpf     ! flag for general purpose library
  logical            :: flagbreakup ! breakup flag
  logical            :: flaghigh    ! flag for high energies ( > 20 MeV)
  character(len=132) :: endffile    ! name of ENDF file
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for input of ENDF limits, switches and tolerances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer   :: maxrp    ! maximum number of residual products
  integer   :: NMTmax   ! maximum number of MT numbers
  real(sgl) :: cuteps   ! energy shift at MT5 cutoff energy (in eV)
  real(sgl) :: disclim  ! limit for specific MT numbers for discrete levels
  real(sgl) :: Eswitch  ! energy where ENDF-6 representation is switched (in MeV)
  real(sgl) :: Eswitch4 ! energy where MF4 representation is switched (in MeV)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for ENDF input for MF1
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=10)  :: identifier ! library identifier
  character(len=11)  :: lab        ! laboratory
  character(len=33)  :: author     ! author
  character(len=132) :: endftext   ! file with MF1 information
  real(sgl)          :: diffweight ! differential weight to be put in MF1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for input of specific ENDF data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(nummf,nummt)            :: adopt       ! logical for existence of MF information (per MT)
  character(len=132), dimension(nummf,nummt) :: adoptfile   ! name of library for MF information (per MT)
  character(len=132)                         :: background  ! file with background cross sections
  integer                                    :: lssfinp     ! 0: URR cross section from MF2, 1: URR cross section
  integer                                    :: urrcomp     ! mode for competition in the URR, 0: none, 1:MT4, 2:all
  integer                                    :: urrmode     ! 0: no URR, 1: URR from TALYS, 2: URR from data library
  real(sgl), dimension(nummf,nummt)          :: Eahigh      ! upper energy of MT values to be adopted
  real(sgl), dimension(nummf,nummt)          :: Ealow       ! lower energy of MT values to be adopted
  real(sgl)                                  :: urrenergy   ! upper energy of the URR in MeV
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for input of ENDF covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                   :: flagcovar   ! flag for covariances
  logical                   :: flagcovrp   ! flag for covariance of residual production cross sections
  logical, dimension(nummt) :: flagcross   ! flag for covariance cross-channel covariance data
  logical                   :: flagintercor! flag for inter-MT covariance data
  logical                   :: flagparcov  ! flag to include covariances for MT60
  integer                   :: covdiscrete ! number of discrete inelastic levels with covariances
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
  integer, dimension(nummf,nummt,numint) :: INTER ! interpolation scheme
  integer, dimension(nummf,nummt,numint) :: NBT   ! separation value for interpolation scheme
  integer, dimension(nummf,nummt)        :: NK    ! number of subsections
  integer, dimension(nummf,nummt)        :: NP    ! number of incident energies
  integer, dimension(nummf,nummt)        :: NR    ! number of interpolation ranges
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF1
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(nummt)      :: LNU     ! polynomial(LNU=1) or tabulated (LNU=2) representation
  integer                        :: NCnubar ! number of terms in polynomial expansion
  integer                        :: NNF     ! number of precursor families
  real(sgl), dimension(numen6)   :: Cnubar  ! coefficients of the polynomial expansion
  real(sgl)                      :: DEBDEL  ! uncertainty of total energy of delayed beta's
  real(sgl)                      :: DEFR    ! uncertainty of kinetic energy of fragments
  real(sgl)                      :: DEGD    ! uncertainty of total energy of delayed gamma rays
  real(sgl)                      :: DEGP    ! uncertainty of total energy of prompt gamma rays
  real(sgl)                      :: DEND    ! uncertainty of kinetic energy of delayed fission neutrons
  real(sgl)                      :: DENP    ! uncertainty of kinetic energy of prompt fission neutrons
  real(sgl)                      :: DENU    ! uncertainty of energy of neutrinos
  real(sgl)                      :: DERN    ! uncertainty of total energy minus energy of neutrinos
  real(sgl)                      :: DET     ! uncertainty of sum of all partial energies
  real(sgl), dimension(3,numen6) :: E1      ! incident energy for MF1 (in ENDF-6 format)
  real(sgl)                      :: EBDEL   ! total energy of delayed beta's
  real(sgl)                      :: EFR     ! kinetic energy of fragments
  real(sgl)                      :: EGD     ! total energy of delayed gamma rays
  real(sgl)                      :: EGP     ! total energy of prompt gamma rays
  real(sgl)                      :: END1    ! kinetic energy of delayed fission neutrons
  real(sgl)                      :: ENP     ! kinetic energy of prompt fission neutrons
  real(sgl)                      :: ENU     ! energy of neutrinos
  real(sgl)                      :: ERN     ! total energy minus energy of neutrinos
  real(sgl)                      :: ET      ! sum of all partial energies
  real(sgl), dimension(numen6)   :: lambda  ! decay constant for precursor
  real(sgl), dimension(3,numen6) :: nubar   ! average number of neutrons per fission
  real(sgl)                      :: nuspont ! nu for spontaneous fission
!
! make1
!
  character(len=6)                      :: ENDATE ! date string
  character(len=10)                     :: DDATE  ! date of distribution
  character(len=10)                     :: EDATE  ! date of evaluationon
  character(len=10)                     :: RDATE  ! date of revision
  character(len=11)                     :: ALAB   ! Lab
  character(len=11)                     :: ZSYMAM ! string
  character(len=21)                     :: REF    ! reference
  character(len=33)                     :: AUTH   ! author(s)
  character(len=66)                     :: HSUB1  ! string
  character(len=66)                     :: HSUB2  ! string
  character(len=66)                     :: HSUB3  ! string
  character(len=66),  dimension(numtxt) :: txt    ! text line
  integer                               :: LDRV   ! special evaluation flag
  integer                               :: LFI    ! flag for fission
  integer                               :: LREL   ! release number
  integer                               :: LRP    ! flag that indicates presence of file 2 (resonances)
  integer, dimension(nummf,nummt)       :: MODN   ! modification indicator
  integer, dimension(nummf,nummt)       :: NC     ! number of lines for MF/MT number
  integer                               :: NFOR   ! library format
  integer                               :: NMOD   ! modification number of evaluation
  integer                               :: NLIB   ! library identifier
  integer                               :: NSUB   ! sub-library number
  integer                               :: ntxt   ! number of text lines
  integer                               :: NVER   ! library version number
  integer                               :: NWD    ! number of text lines
  integer                               :: NXC    ! number of directory lines
  real(sgl)                             :: AWI    ! mass of projectile in neutron units
  real(sgl)                             :: ELIS   ! excitation energy of target nucleus
  real(sgl)                             :: EMAX   ! upper limit of energy range for evaluation
  real(sgl)                             :: STA    ! target stability flag
  real(sgl)                             :: TEMP   ! target temperature
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF2
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                             :: IFG7   ! flag for gamma width
  integer                                             :: INTERU ! URR interpolation scheme
  integer, dimension(numjres)                         :: KBK7   ! R-matrix background parameter
  integer, dimension(numres,numlres,numjres)          :: kINT   ! interpolation scheme
  integer, dimension(numjres)                         :: KPS7   ! non-hard-sphere phase shift parameter
  integer                                             :: KRL7   ! flag for relativistic kinematics
  integer                                             :: KRM7   ! flag of formula to be used (MLBW, RM, etc.) for R-matrix
  integer, dimension(numres)                          :: LAD    ! flag to indicate computation of angular distributions
  integer                                             :: LFW    ! flag for average fission width
  integer, dimension(numres,numlres)                  :: Lres   ! l-value
  integer, dimension(numres)                          :: LRF    ! representation indicator
  integer, dimension(numres)                          :: LRU    ! flag for resolved/unresolved
  integer, dimension(numres,numlres)                  :: LRX    ! flag to indicate competitive width
  integer, dimension(numres)                          :: LSSF   ! flag for interpretation of MF3 cross sections
  integer, dimension(numres)                          :: NAPS   ! flag for channel radius and scattering radius
  integer                                             :: NBTU   ! separation value for URR interpolation scheme
  integer, dimension(numjres)                         :: NCH7   ! number of channels per J
  integer                                             :: NER    ! number of resonance energy ranges
  integer, dimension(numres,numlres,numjres)          :: NEu    ! number of energy points
  integer                                             :: NIS    ! number of isotopes in the material
  integer, dimension(numres,numlres)                  :: NJS    ! number of J-states for a particular l-value
  integer                                             :: NJS7   ! number of J-states for a particular l-value
  integer, dimension(numres)                          :: NLS    ! number of l-values
  integer, dimension(numres)                          :: NLSC   ! number of l-values for convergence of angular distribution
  integer                                             :: NPP7   ! total number of particle pairs
  integer                                             :: NPU    ! number of URR incident energies
  integer, dimension(numres)                          :: NRO    ! flag for energy dependence of scattering radius
  integer, dimension(numres,numlres)                  :: NRS    ! number of resolved resonances per l-value
  integer, dimension(numjres)                         :: NRS7   ! number of resonances
  integer                                             :: NRU    ! number of URR interpolation ranges
  integer, dimension(numjres)                         :: NX7    ! number of lines required for all resonances
  real(sgl), dimension(numres,numlres,numnrs)         :: AJ     ! spin of the resonance
  real(sgl), dimension(numPP)                         :: AJ7    ! spin
  real(sgl), dimension(numres,numlres,numjres)        :: AJU    ! spin
  real(sgl), dimension(numres,numlres,numjres)        :: AMUF   ! number of degrees of freedom for fission width distributio
  real(sgl), dimension(numres,numlres,numjres)        :: AMUG   ! number of degrees of freedom for gamma width distribution
  real(sgl), dimension(numres,numlres,numjres)        :: AMUN   ! number of degrees of freedom for neutron width distributio
  real(sgl), dimension(numres,numlres,numjres)        :: AMUX   ! number of degrees of freedom for competitive width distrib
  real(sgl), dimension(numres)                        :: AP     ! scattering radius
  real(sgl), dimension(numres,numenin)                :: APE    ! scattering radius
  real(sgl), dimension(numjres,numch7)                :: APE7   ! effective channel radius
  real(sgl), dimension(numres,numlres)                :: APL    ! l-dependent scattering radius
  real(sgl), dimension(numjres,numch7)                :: APT7   ! true channel radius
  real(sgl), dimension(numres,numlres)                :: AWRI   ! ratio of isotope mass to neutron
  real(sgl), dimension(numjres,numch7)                :: BND7   ! boundary condition for this channel
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: D      ! average level spacing for resonances with spin J
  real(sgl), dimension(numres,numenin)                :: E2     ! incident energy for MF2 (in ENDF-6 format)
  real(sgl), dimension(numres)                        :: EH     ! boundary for resonance range
  real(sgl), dimension(numres)                        :: EL     ! boundary for resonance range
  real(sgl), dimension(numres,numlres,numnrs)         :: Er     ! resonance energy in LAB system
  real(sgl), dimension(numjres,numch7)                :: ER7    ! energy of resonance in eV
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: Es     ! energy of energy-dependent width
  real(sgl), dimension(numres,numlres)                :: QX     ! Q-value to be added to C.M. incident energy
  real(sgl), dimension(numres)                        :: SPI    ! target spin
  real(sgl), dimension(numjres,numnrs,numch7)         :: GAM7   ! channel width or reduced width
  real(sgl), dimension(numres,numlres,numnrs)         :: GF     ! fission width of the resonance
  real(sgl), dimension(numres,numlres,numnrs)         :: GFA    ! first partial fission width
  real(sgl), dimension(numres,numlres,numnrs)         :: GFB    ! second partial fission width
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: GFu    ! average fission width
  real(sgl), dimension(numres,numlres,numnrs)         :: GG     ! gamma width of the resonance
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: GGu    ! average radiation width
  real(sgl), dimension(numres,numlres,numnrs)         :: GN     ! neutron width of the resonance
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: GN0    ! average reduced neutron width
  real(sgl), dimension(numres,numlres,numnrs)         :: GT     ! total width of the resonance
  real(sgl), dimension(numres,numlres,numjres,numnrs) :: GX     ! average competitive reaction width
  real(sgl), dimension(numPP)                         :: IA7    ! spin of first particle in pair
  real(sgl), dimension(numPP)                         :: IB7    ! spin of second particle in pair
  real(sgl), dimension(numjres,numch7)                :: L7     ! orbital angular momentum
  real(sgl), dimension(numPP)                         :: MA7    ! mass of first particle in pair
  real(sgl), dimension(numPP)                         :: MB7    ! mass of second particle in pair
  real(sgl), dimension(numPP)                         :: MTP7   ! reaction type for this particle pair
  real(sgl), dimension(numPP)                         :: PA7    ! parity of first particle in pair
  real(sgl), dimension(numPP)                         :: PB7    ! parity of second particle in pair
  real(sgl), dimension(numPP)                         :: PJ7    ! parity
  real(sgl), dimension(numPP)                         :: PNT7   ! flag for penetrability
  real(sgl), dimension(numjres,numch7)                :: PPI7   ! particle-pair number for this channel
  real(sgl), dimension(numPP)                         :: QP7    ! Q-value for this particle pair
  real(sgl), dimension(numjres,numch7)                :: SCH7   ! channel spin
  real(sgl), dimension(numPP)                         :: SHF7   ! flag for shift factor
  real(sgl), dimension(numPP)                         :: ZA7    ! charge of first particle in pair
  real(sgl), dimension(numPP)                         :: ZB7    ! charge of second particle in pair
!
! make2
!
  real(sgl) :: ABN ! abundance
  real(sgl) :: ZAI ! (Z,A) designation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF3
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(nummt)               :: NE3adopt ! number of adopted incident energies
  integer, dimension(nummt)               :: NE3res   ! number of energies in resonance range
  real(sgl), dimension(nummtres,0:numen6) :: E3adopt  ! incident energy adopted from other library
  real(sgl), dimension(nummtres,0:numen6) :: E3res    ! energy in resonance range
  real(sgl), dimension(nummtres,0:numen6) :: xs3adopt ! cross section adopted from other library
  real(sgl), dimension(nummtres,0:numen6) :: xs3res   ! cross section in resonance range
!
! make3total
!
  integer, dimension(nummt)               :: LR3   ! breakup flag
  integer                                 :: NMT   ! total number of MT sections
  real(sgl), dimension(nummt,0:numen6)    :: E3    ! incident energy for MF3 (in ENDF-6 format)
  real(sgl), dimension(nummt)             :: EthMT ! threshold energy
  real(sgl), dimension(nummt)             :: QI    ! Q-value (in ENDF-6 format)
  real(sgl), dimension(nummt)             :: QM    ! Q-value (in ENDF-6 format)
  real(sgl), dimension(nummt,0:numen6)    :: xs    ! cross section
!
! make3eaf
!
  character(len=66), dimension(nummt) :: eafstring ! string with reaction information for EAF format
!
! make3partial
!
  integer, dimension(nummt)           :: LFS3     ! isomeric state number (EAF only)
  character(len=66), dimension(nummt) :: mtstring ! string with reaction information
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF4
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                 :: NE4r ! number of incident energies (MF4 only)
  integer                                 :: NEhr ! number of incident energies (MF4 only)
  integer, dimension(numen4)              :: NL4r ! number of Legendre coefficients
  integer, dimension(numen4+1)            :: NP4r ! number of angles (MF4 only)
  integer, dimension(numen4+1)            :: NR4r ! number of interpolation ranges (MF4 only)
  real(sgl), dimension(numen4+1)          :: E4hr ! incident energy for MF4 (in ENDF-6 format)
  real(sgl), dimension(numen4)            :: E4r  ! incident energy for MF4 (in ENDF-6 format)
  real(sgl), dimension(numen4+1,numang+3) :: f4r  ! angular distribution
  real(sgl), dimension(numen4,0:numl)     :: legr ! Legendre coefficients (in ENDF-6 format)
  real(sgl), dimension(numen4+1,numang+3) :: x4r  ! cosine of the angle
!
! make4elastic
!
  integer, dimension(numen4,numint)       :: INTER4 ! interpolation scheme
  integer, dimension(numint)              :: INTERh ! interpolation scheme
  integer                                 :: LCT    ! LAB/CM flag
  integer                                 :: LI4    ! isotropy flag
  integer                                 :: LTT    ! representation
  integer                                 :: LVT    ! specification of transformation matrix
  integer, dimension(numen4,numint)       :: NBT4   ! separation value for interpolation scheme
  integer, dimension(numint)              :: NBTh   ! separation value for interpolation scheme
  integer                                 :: NE     ! number of incident energies
  integer                                 :: NEh    ! number of incident energies (MF4 only)
  integer, dimension(6,nummt,numen4)      :: NL     ! Legendre order or number of cosines
  integer, dimension(numen4)              :: NP4    ! number of incident energies
  integer, dimension(numen4)              :: NR4    ! number of interpolation ranges
  integer                                 :: NRh    ! number of interpolation ranges
  real(sgl), dimension(numen4)            :: E4     ! incident energy for MF4 (in ENDF-6 format)
  real(sgl), dimension(numen4+1)          :: E4h    ! incident energy for MF4 (in ENDF-6 format)
  real(sgl), dimension(numen4+1,numang+3) :: f4     ! angular distribution
  real(sgl), dimension(numen4,0:numl)     :: leg    ! Legendre coefficients (in ENDF-6 format)
  real(sgl), dimension(numen4+1,numang+3) :: x4     ! cosine of the angle
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF5
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(numsecea,numint)              :: INTER5   ! interpolation scheme
  integer, dimension(numsecea,numint)              :: INTER5e  ! interpolation scheme
  integer, dimension(numsecea,numenin,numint)      :: INTER5e2 ! interpolation scheme
  integer, dimension(numsecea)                     :: LF       ! flag for energy distribution law
  integer, dimension(numsecea,numint)              :: NBT5     ! separation value for interpolation scheme
  integer, dimension(numsecea,numint)              :: NBT5e    ! separation value for interpolation scheme
  integer, dimension(numsecea,numenin,numint)      :: NBT5e2   ! separation value for interpolation scheme
  integer, dimension(numsecea)                     :: NE5e     ! number of incident energies for distribution
  integer, dimension(numsecea,numenin)             :: NF       ! number of secondary energy points
  integer, dimension(numsecea)                     :: NP5      ! number of incident energies
  integer, dimension(numsecea)                     :: NR5      ! number of interpolation ranges
  integer, dimension(numsecea)                     :: NR5e     ! number of interpolation ranges
  integer, dimension(numsecea,numenin)             :: NR5e2    ! number of interpolation ranges
  real(sgl), dimension(numsecea,numenin)           :: E5       ! incident energy for MF5 (in ENDF-6 format)
  real(sgl), dimension(numsecea,numenin)           :: E5p      ! incident energy for which tabulated distribution is given
  real(sgl)                                        :: EFH      ! constant used in the energy-dependent fission neutron spectrum
  real(sgl)                                        :: EFL      ! constant used in the energy-dependent fission neutron spectrum
  real(sgl), dimension(numsecea,numenin,10*numen2) :: gE5      ! energy-spectrum values
  real(sgl), dimension(numsecea,numenin)           :: pE       ! fractional part of cross section
  real(sgl), dimension(numsecea,10*numen2)         :: TM5      ! Effective (7) or maximum (12) temperature  FNS parameter
  real(sgl), dimension(numsecea)                   :: U        ! constant for upper energy limit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF6
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                :: LIDP   ! indicates that particles are identical
  integer                :: LTP    ! representation
  real(sgl)              :: SPIpar ! particle spin
!
! make6partial
!
  integer, dimension(numsec,numint)                    :: INTER6ea ! interpolation scheme
  integer, dimension(numsec,numint)                    :: INTER6y  ! interpolation scheme
  integer, dimension(numsec)                           :: LANG     ! flag for angular representation
  integer, dimension(numsec)                           :: LAW      ! flag for distribution function
  integer, dimension(numsec)                           :: LEP      ! interpolation scheme for secondary energy
  integer, dimension(numsec)                           :: LIP      ! product modifier flag
  integer, dimension(numsec,numenin)                   :: NA       ! number of angular parameters
  integer, dimension(numsec,numint)                    :: NBT6ea   ! separation value for interpolation scheme
  integer, dimension(numsec,numint)                    :: NBT6y    ! separation value for interpolation scheme
  integer, dimension(numsec,numenin)                   :: ND       ! number of discrete energies
  integer, dimension(numsec)                           :: NE6ea    ! number of incident energies for distribution
  integer, dimension(numsec,numenin)                   :: NEP      ! number of secondary energy points
  integer, dimension(numsec)                           :: NP6y     ! number of incident energies for yields
  integer, dimension(numsec)                           :: NR6ea    ! number of interpolation ranges for distribution
  integer, dimension(numsec)                           :: NR6y     ! number of interpolation ranges for yields
  integer, dimension(numsec,numenin)                   :: NW       ! number of words
  real(sgl), dimension(numsec)                         :: AWP      ! product mass
  real(sgl), dimension(numsecea,numenspec+1,40*numen2) :: b6       ! energy-angle values
  real(sgl), dimension(numsec,numenin)                 :: E6       ! incident energy (in ENDF-6 format) for distribution
  real(sgl), dimension(numsec,numenin)                 :: Ey       ! incident energy for yields (in ENDF-6 format)
  real(sgl), dimension(numsec,numenin)                 :: Y        ! product yield (in ENDF-6 format)
  real(sgl), dimension(numsec)                         :: ZAP      ! product identifier
!
! make6mt5
!
  logical, dimension(numsec)                          :: flagrec ! flag to state that b6 element concerns recoil grid
  integer                                             :: kpart   ! section number for particles
  real(sgl), dimension(numenspec+1,40*numen2)         :: b6gam   ! energy-angle values for photons
  real(sgl), dimension(numsec,numenspec+1,2*numenrec) :: b6rec   ! energy-angle values for recoils
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF8-10
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(nummt,numsec)   :: LMF  ! file number for information
  integer, dimension(nummt,numsec)   :: MATP ! material number for reaction product
  integer                            :: NDec ! number of decay branches
  integer                            :: NO   ! flag for decay information
  integer, dimension(nummt)          :: ZAPi ! designation of final nucleus
  real(sgl), dimension(nummt,numsec) :: ZAPr ! designation of final nucleus
!
! make10
!
  integer, dimension(nummt,numiso,numint)    :: INTER10  ! interpolation scheme
  integer, dimension(numsec,numint)          :: INTERZA  ! interpolation scheme
  integer, dimension(numsec)                 :: IZAP     ! second IZAP-number
  integer, dimension(nummt,numsec)           :: LFS      ! final state number
  integer, dimension(numsec)                 :: LFSZA    ! final state number
  integer, dimension(nummt,numiso,numint)    :: NBT10    ! separation value for interpolation scheme
  integer, dimension(numsec,numint)          :: NBTZA    ! separation value for interpolation scheme
  integer, dimension(nummt,numiso)           :: NP10     ! number of incident energies
  integer, dimension(numsec)                 :: NPZA     ! number of incident energies
  integer, dimension(nummt,numiso)           :: NR10     ! number of interpolation ranges
  integer, dimension(numsec)                 :: NRZA     ! number of interpolation ranges
  integer, dimension(nummt)                  :: NSt      ! number of final states
  integer                                    :: NZA      ! number of nuclides
  integer, dimension(numsec)                 :: XMFZA    ! second MF-number
  real(sgl), dimension(nummt,numiso,numenin) :: E10      ! incident energy (in ENDF-6 format)
  real(sgl), dimension(numsec,numenin)       :: E10ZA    ! incident energy (in ENDF-6 format)
  real(sgl), dimension(nummt,numsec)         :: ELFS     ! excitation energy of final state
  real(sgl), dimension(numsec)               :: ErpZAiso ! energy of isomer
  real(sgl), dimension(numsec)               :: EthZA    ! threshold energy
  real(sgl), dimension(nummt,numiso)         :: QIiso    ! Q-value for isomer (in ENDF-6 format)
  real(sgl), dimension(numsec)               :: QIZA     ! Q-value (in ENDF-6 format)
  real(sgl), dimension(numsec)               :: QMZA     ! Q-value (in ENDF-6 format)
  real(sgl), dimension(nummt,numiso,numenin) :: xsiso    ! cross section for isomer (in ENDF-6 format)
  real(sgl), dimension(numsec,numenin)       :: xsrpZA   ! cross section for residual production (in ENDF-6 format)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF12-15
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(nummt,numgam,numint)        :: INTERg      ! interpolation scheme
  integer, dimension(nummt)                      :: LG12        ! type setters
  integer, dimension(nummt)                      :: LO12        ! type setters
  integer, dimension(nummt,numgam,numint)        :: NBTg        ! separation value for interpolation scheme
  integer, dimension(nummt,numgam)               :: NPg         ! number of incident energies
  integer, dimension(nummt,numgam)               :: NRg         ! number of interpolation ranges
  integer, dimension(nummt)                      :: LP12        ! origin of photons
  integer, dimension(nummt,numgam)               :: LPg         ! primary photon flag
  integer, dimension(nummt,numgam)               :: LFg         ! photo energy distribution law
  integer, dimension(nummt)                      :: NS12        ! number of levels below the present one
  integer, dimension(nummt)                      :: NT12        ! number of transitions for which data are given
  real(sgl), dimension(0:numchan,numenin)        :: E12         ! incident energy (in ENDF-6 format)
  real(sgl), dimension(0:numchan,numgam,numenin) :: Eg          ! gamma energy
  real(sgl), dimension(0:numchan,numgam)         :: Egk         ! gamma energy
  real(sgl), dimension(nummt,numgam)             :: ES12        ! energy of level
  real(sgl), dimension(0:numchan,numgam)         :: Esk         ! starting level (in ENDF-6 format)
  real(sgl), dimension(nummt)                    :: ESNS        ! energy of mother level
  real(sgl), dimension(nummt,numgam)             :: TP12        ! probability of direct transition
  real(sgl), dimension(0:numchan,numenin)        :: xsgtotyield ! total discrete photon multiplicity
  real(sgl), dimension(0:numchan,numgam,numenin) :: xsgyield    ! gamma-ray multiplicity (in ENDF-6 format)
!
! make13
!
  real(sgl), dimension(0:numchan,numenin)        :: E13       ! incident energy (in ENDF-6 format)
  real(sgl), dimension(0:numchan,numgam,numenin) :: xsg       ! gamma-ray cross section (in ENDF-6 format)
  real(sgl), dimension(0:numchan,numenin)        :: xsgtot    ! total discrete photon production cross section
!
! make14
!
  integer, dimension(nummt)       :: LI14 ! isotropy flag
!
! make15
!
  integer, dimension(numsecg,numint)             :: INTER15   ! interpolation scheme
  integer, dimension(numsecg,numint)             :: INTER15g  ! interpolation scheme
  integer, dimension(numsecg,numenin,numint)     :: INTER15ge ! interpolation scheme
  integer, dimension(numsecg,numint)             :: NBT15     ! separation value for interpolation scheme
  integer, dimension(numsecg,numint)             :: NBT15g    ! separation value for interpolation scheme
  integer, dimension(numsecg,numenin,numint)     :: NBT15ge   ! separation value for interpolation scheme
  integer, dimension(numsecg)                    :: NE15g     ! number of incident energies for distribution
  integer, dimension(numsecg)                    :: NP15      ! number of incident energies
  integer, dimension(numsecg,numenin)            :: NP15ge    ! number of secondary energy point
  integer, dimension(numsecg)                    :: NR15      ! number of interpolation ranges
  integer, dimension(numsecg)                    :: NR15g     ! number of interpolation ranges
  integer, dimension(numsecg,numenin)            :: NR15ge    ! number of interpolation ranges
  real(sgl), dimension(numsecg,numenin)          :: E15       ! incident energy (in ENDF-6 format)
  real(sgl), dimension(numsecg,numenin,3*numen2) :: E15ge     ! secondary energy
  real(sgl), dimension(numsecg,numenin)          :: EPy       ! incident energy for probabilities (in ENDF-6 format)
  real(sgl), dimension(numsecg,numenin,3*numen2) :: ge        ! gamma distribution
  real(sgl), dimension(numsecg,numenin)          :: Pg        ! probability (in ENDF-6 format)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for MF31-40
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                              :: LB31 ! flag for meaning of numbers
  integer                              :: LS31 ! symmetry flag
  integer                              :: NC31 ! number of NC-type sub-subsections
  integer                              :: NE31 ! number of energies in energy array
  integer                              :: NI31 ! number of NI-type sub-subsections
  integer                              :: NL31 ! number of subsections
  integer                              :: NT31 ! total number of entries
  real(sgl), dimension(10*numencovtot) :: b31  ! covariance matrix element
!
! read32
!
  character(len=4), dimension(numrescov*numrescov/2,14) :: covdigit ! elements of compact covariance format
  integer, dimension(numrescov*numrescov/2,2)           :: covix32  ! covariance index
  integer, dimension(numres)                            :: ISR      ! flag for presence of scattering radius uncertainty
  integer, dimension(numres)                            :: LCOMP    ! compatibility flag
  integer                                               :: LRX32    ! flag to indicate competitive width
  integer, dimension(numres)                            :: MLS      ! number of DAP points
  integer                                               :: MPAR     ! number of parameters per resonance
  integer                                               :: MPARURR  ! number of parameters per resonance for URR
  integer                                               :: N32      ! number of values for MF32
  integer                                               :: N32URR   ! number of points for URR
  integer                                               :: NDIGIT   ! integer for compact covariance format
  integer, dimension(numres,numlres)                    :: NJS32    ! number of j-values
  integer                                               :: NLRS     ! number of subsections with long-range covariance
  integer, dimension(numres)                            :: NLS32    ! number of l-values
  integer                                               :: NM       ! integer for compact covariance format
  integer                                               :: NNN      ! integer for compact covariance format
  integer                                               :: NPARURR  ! number of parameters per resonance for URR
  integer                                               :: NRB      ! number of resonances
  integer                                               :: NSRS     ! number of subsections with covariances
  real(sgl), dimension(numres,numlres,numjres)          :: AJ32     ! spin of the resonance
  real(sgl)                                             :: APLQX    ! l-dependent scattering radius
  real(sgl), dimension(5*numjres*(numjres+1))           :: b32URR   ! covariance matrix elemen
  real(sgl), dimension(numres,numlres,numjres)          :: D32      ! MF32 resonance parameter
  real(sgl), dimension(numres)                          :: DAP      ! uncertainty in scattering radius (compact format)
  real(sgl), dimension(numres,numlres,numjres)          :: GF32     ! fission width of the resonance
  real(sgl), dimension(numres,numlres,numjres)          :: GG32     ! gamma width of the resonance
  real(sgl), dimension(numres,numlres,numjres)          :: GNO32    ! neutron width of the resonance
  real(sgl), dimension(numres,numlres,numjres)          :: GX32     ! competitive width of the resonance
  real(sgl), dimension(6*numrescov+numrescov*numrespar*(numrescov*numrespar+1)) :: b32      ! covariance matrix
!
! read33
!
  integer, dimension(numchan,numchan)                     :: LBread    ! flag for meaning of numbers
  integer                                                 :: LB8read   ! flag for meaning of numbers
  integer, dimension(numchan,numchan)                     :: LSread    ! symmetry flag
  integer, dimension(numchan,numchan)                     :: MAT1read  ! second MAT-number
  integer, dimension(numchan,numchan)                     :: MT33read  ! second MT-number
  integer, dimension(nummt)                               :: MTLread   ! lumped reaction identifier
  integer, dimension(numchan,numchan)                     :: NC33read  ! number of NC-type sub-subsections
  integer, dimension(numchan)                             :: NE8read   ! number of entries for LB=8 section
  integer, dimension(numchan,numchan)                     :: NE33read  ! number of energies in energy array
  integer, dimension(numchan,numchan)                     :: NI33read  ! number of NI-type sub-subsections
  integer, dimension(nummf,nummt)                         :: NL33read  ! number of subsections
  integer, dimension(numchan)                             :: NT8read   ! total number of entries
  integer, dimension(numchan,numchan)                     :: NT33read  ! total number of entries
  real(sgl), dimension(numchancov,numchancov,numencovtot) :: b33read   ! covariance matrix element
  real(sgl), dimension(numchan,numencovtot)               :: b33MTread ! covariance matrix element
  real(sgl), dimension(numchan,numencovtot)               :: b8read    ! covariance matrix element
  real(sgl), dimension(numchan,numchan)                   :: XLFS1read ! second state discrete level number
  real(sgl), dimension(numchan,numchan)                   :: XMF1read  ! second MF-number
!
! make33
!
  integer, dimension(numchan,numchan)                     :: IZAP1     ! second IZAP-number
  integer, dimension(numchan,numchan)                     :: LB        ! flag for meaning of numbers
  integer                                                 :: LB8       ! flag for meaning of numbers
  integer, dimension(numsec,numchan)                      :: LBZA      ! flag for meaning of numbers
  integer, dimension(numchan,numchan)                     :: LS        ! symmetry flag
  integer, dimension(numsec,numchan)                      :: LSZA      ! symmetry flag
  integer, dimension(numchan,numchan)                     :: MAT1      ! second MAT-number
  integer, dimension(numsec,numchan)                      :: MAT1ZA    ! second MAT-number
  integer, dimension(numchan,numchan)                     :: MT33      ! second MT-number
  integer, dimension(numsec,numchan)                      :: MT33ZA    ! second MT-number
  integer, dimension(nummt)                               :: MTL       ! lumped reaction identifier
  integer, dimension(numchan,numchan)                     :: NC33      ! number of NC-type sub-subsections
  integer, dimension(numsec,numchan)                      :: NC33ZA    ! number of NC-type sub-subsections
  integer, dimension(numchan,numchan)                     :: NE33      ! number of energies in energy array
  integer, dimension(numsec,numchan)                      :: NE33ZA    ! number of energies in energy array
  integer, dimension(numchan)                             :: NE8       ! number of entries for LB=8 section
  integer, dimension(numchan,numchan)                     :: NEC33     ! number of energies in energy array
  integer, dimension(numchan,numchan)                     :: NI33      ! number of NI-type sub-subsections
  integer, dimension(numsec,numchan)                      :: NI33ZA    ! number of NI-type sub-subsections
  integer, dimension(nummf,nummt)                         :: NL33      ! number of subsections
  integer, dimension(numchan,numchan)                     :: NT33      ! total number of entries
  integer, dimension(numsec,numchan)                      :: NT33ZA    ! total number of entries
  integer, dimension(numchan)                             :: NT8       ! total number of entries
  real(sgl), dimension(numchancov,numchancov,numencovtot) :: b33       ! covariance matrix element
  real(sgl), dimension(numchan,numencovtot)               :: b33MT     ! covariance matrix element
  real(sgl), dimension(numchan,numencovtot)               :: b33ZA     ! covariance matrix element
  real(sgl), dimension(numchan,numencovtot)               :: b8        ! covariance matrix element
  real(sgl)                                               :: Emincov   ! minimum energy for covariance information
  real(sgl), dimension(numchan,numchan)                   :: XLFS1     ! second state discrete level number
  real(sgl), dimension(numchan,numchan)                   :: XMF1      ! second MF-number
!
! make34
!
  integer                                         :: LB34  ! flag for meaning of numbers
  integer                                         :: LS34  ! symmetry flag
  integer                                         :: LTT34 ! representation
  integer                                         :: LVT34 ! specification of transformation matrix
  integer                                         :: MAT34 ! material number for MF34
  integer                                         :: MT34  ! MT number for MF34
  integer, dimension(numleg,numleg)               :: NE34  ! number of energies in energy array
  integer, dimension(numleg,numleg)               :: NI34  ! number of NI-type sub-subsections
  integer                                         :: NL34  ! number of Legendre coefficients with covariances
  integer                                         :: NL341 ! number of Legendre coefficients with covariances
  integer                                         :: NMT34 ! total number of MT sections
  integer, dimension(numleg,numleg)               :: NT34  ! total number of entries
  real(sgl), dimension(numleg,numleg,numencovtot) :: b34   ! covariance matrix element
!
! read35
!
  integer, parameter                            :: numencov35=750 ! number of incident energies for FNS covariances
  integer                                       :: LB35           ! flag for meaning of numbers
  integer                                       :: LS35           ! symmetry flag
  integer, dimension(numencov35)                :: NE35           ! number of energies in energy array
  integer, dimension(numencov35)                :: NT35           ! total number of entries
  real(sgl), dimension(numencov35, numencovtot) :: b35            ! covariance matrix element
  real(sgl), dimension(numencov35)              :: E35b           ! start energy of block
  real(sgl), dimension(numencov35)              :: E35e           ! end energy of block
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for initialization of ENDF variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=11)               :: blank1  ! blank string
  character(len=66)               :: blank2  ! blank string
  character(len=30)               :: CONT    ! ENDF-6 format
  character(len=30)               :: FEND    ! ENDF-6 format
  character(len=30)               :: HEAD    ! ENDF-6 format
  character(len=30)               :: INT3    ! ENDF-6 format
  character(len=30)               :: MEND    ! ENDF-6 format
  character(len=30)               :: SEND    ! ENDF-6 format
  character(len=30)               :: TEND    ! ENDF-6 format
  character(len=30)               :: TEXT    ! ENDF-6 format
  character(len=30)               :: VALUE   ! ENDF-6 format
  character(len=80)               :: rec0    ! first line of ENDF-6 file
  integer                         :: iclean  ! number of cleaned points
  integer                         :: idnum   ! number of different exclusive cross sections
  integer                         :: LIS     ! state number of target nucleus
  integer                         :: LISO    ! isomeric state number
  integer                         :: MAT     ! MAT number
  logical, dimension(nummf)       :: mfexist ! flag for existence of MF-number
  integer, dimension(nummt)       :: MTid    ! channel identifier for MT-number
  integer                         :: MTinel  ! MT-number for inelastic scattering
  integer, dimension(nummt)       :: MTnum   ! channel identifier for MT-number (ENDF format)
  logical, dimension(nummf,nummt) :: mtexist ! flag for existence of MT-number
  integer                         :: MTmax   ! highest MT number for exclusive channels
  real(sgl)                       :: AWR     ! standard mass parameter
  real(sgl)                       :: ZA      ! standard charge parameter
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for total cross sections in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                         :: numcut6  ! number of energies before high-energy format for total cross sections
  integer                         :: numinc6  ! number of incident energies for total cross sections
  real(sgl), dimension(0:numen6)  :: eninc6   ! incident energy for total cross section
  real(sgl)                       :: Rprime   ! potential scattering radius
  real(dbl), dimension(0:numenin) :: xsany    ! (x,anything) cross section (MT5)
  real(dbl), dimension(0:numen6)  :: xsany6   ! (x,anything) cross section (MT5)
  real(dbl), dimension(0:numen6)  :: xselas   ! total elastic cross section
  real(dbl), dimension(0:numen6)  :: xselas6  ! total elastic cross section
  real(dbl), dimension(numen6)    :: xsnonel  ! nonelastic cross section
  real(dbl), dimension(0:numen6)  :: xsnonel6 ! nonelastic cross section
  real(dbl), dimension(numen6)    :: xstot6   ! total cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for partial cross sections in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numchan,0:numlevels)                        :: isoexist    ! flag for existence of isomer
  logical                                                          :: flagfission ! flag for fission
  character(len=54), dimension(-1:numchan,0:numlevels)             :: reacstring  ! string with reaction information
  integer, dimension(-5:numchan)                                   :: idchannel   ! identifier for channel
  integer, dimension(0:numchan,0:numpar,numenspec)                 :: nbeg        ! first outgoing energy
  integer, dimension(0:numchan,numenspec)                          :: nbegrec     ! first outgoing energy
  integer, dimension(0:numchan,0:numpar,numenspec)                 :: nend        ! last outgoing energy
  integer, dimension(0:numchan,numenspec)                          :: nendrec     ! last outgoing energy
  integer, dimension(0:numchan)                                    :: Nisomer     ! number of isomers
  integer, dimension(0:numchan,numenspec)                          :: nout        ! number of emission energies
  integer, dimension(0:numchan,numenspec)                          :: noutrec     ! number of recoil energies
  real(sgl), dimension(0:numchan,0:numlevels,numenin)              :: branchiso   ! branching ratio for isomer
  real(sgl), dimension(0:numchan,0:numlevels,0:numlevels)          :: Egammadis   ! gamma energy
  real(sgl), dimension(0:numchan,numenspec,0:numpar,0:numen2)      :: Ehist       ! histogram emission energy
  real(sgl), dimension(0:numchan,numenspec,0:numenrec)             :: Ehistrec    ! histogram recoil energy
  real(sgl), dimension(0:numchan,numenspec,0:numen2)               :: Eout        ! emission energy
  real(sgl), dimension(0:numchan,numenspec)                        :: Eparticles  ! total energy carried away by particles
  real(sgl), dimension(0:numchan,numenspec,0:numenrec)             :: Erec        ! recoil energy
  real(sgl), dimension(0:numchan,numenspec)                        :: Erecav      ! average recoil energy
  real(sgl), dimension(0:numchan,0:numlevels,0:numlevels)          :: Estartdis   ! starting level
  real(sgl), dimension(-5:numchan)                                 :: Ethexcl     ! threshold energy
  real(sgl), dimension(0:numchan,0:numlevels)                      :: Ethexcliso  ! threshold energy for isomer
  real(sgl), dimension(0:numchan,numenspec,0:numpar,0:numen2)      :: f0ex        ! energy distribution for exclusive c
  real(sgl), dimension(0:numchan,numenspec,0:numenrec)             :: f0exrec     ! energy distribution for recoil
  real(sgl), dimension(0:numchan)                                  :: Qexcl       ! Q-value
  real(sgl), dimension(0:numchan,0:numlevels)                      :: Qexcliso    ! Q-value for isomer
  real(sgl), dimension(0:numchan,numenspec,0:numenrec)             :: recexcl     ! exclusive recoils
  real(sgl), dimension(0:numchan,numenspec,0:numpar,0:numen2)      :: specexcl    ! exclusive spectra
  real(sgl), dimension(-5:numchan,0:numenin)                       :: xsexcl      ! exclusive cross section
  real(sgl), dimension(0:numchan,0:numlevels,numenin)              :: xsexcliso   ! exclusive cross section for isomer
  real(sgl), dimension(0:numchan,numengam,0:numlevels,0:numlevels) :: xsgamdis    ! exclusive discrete gamma-ray c
  real(sgl), dimension(0:numchan,numenin)                          :: xsgamexcl   ! exclusive gamma cross section
  real(sgl), dimension(numenin)                                    :: xsnonth     ! sum of all non-thr. reac. except (n,g) and (n,f)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for discrete cross sections in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numpar,0:numlevels)                      :: levexist  ! flag for existence of discrete level
  integer, dimension(0:numpar,0:numlevin,numenang)              :: ncleg     ! number of Legendre coefficients
  integer, dimension(0:numpar)                                  :: ndisc     ! number of discrete levels
  integer, dimension(0:numpar)                                  :: numendisc ! number of incident energies including discrete states
  real(sgl), dimension(0:numpar,0:numlevin,numenang,0:numl)     :: cleg0     ! Legendre coefficients
  real(sgl), dimension(0:numpar,0:numenin+150)                  :: edisc     ! incident energy for discrete level cross sections
  real(sgl), dimension(0:numpar,0:numlevels)                    :: Ethdisc   ! threshold energy
  real(sgl), dimension(0:numpar,0:numlevels)                    :: jdis      ! spin of level
  real(dbl), dimension(0:numpar,0:numlevels)                    :: Qdisc     ! Q-value
  real(sgl), dimension(0:numpar,0:numlevin,0:numenang,0:numang) :: xsang     ! differential cross section
  real(sgl), dimension(0:numpar,0:numenin+150)                  :: xsbin     ! binary cross section
  real(sgl), dimension(0:numpar,0:numenin)                      :: xscont    ! continuum cross section
  real(sgl), dimension(0:numpar,0:numlevels,0:numenin)          :: xsdisc    ! discrete state cross section
  real(sgl), dimension(0:numpar,0:numlevels,0:numenin+150)      :: xsintdisc ! interpolated discrete cross section
  real(sgl), dimension(0:numpar,0:numenin)                      :: xsngn     ! (projectile,gamma-ejectile) cross section
  real(sgl), dimension(numenin)                                 :: xsngnsum  ! total (projectile,gamma-ejectile) cross section
  real(sgl), dimension(0:numpar,numenin)                        :: yieldngn  ! (projectile,gamma-ejectile) particle yield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for URR in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                         :: NLSurr  ! number of l-values
  integer, dimension(numlres)                     :: NJSurr  ! number of j-values
  integer                                         :: Nurr    ! number of incident energies with URR data
  real(sgl), dimension(0:numenin,numlres,numjres) :: Djlurr  ! average level spacing for resonances for j,l value
  real(sgl), dimension(numenin)                   :: Eurr    ! incident energy with URR data
  real(sgl), dimension(numlres,numjres)           :: Jurr    ! spin value
  real(sgl), dimension(0:numenin,numlres,numjres) :: GFjlurr ! average fission width for j,l value
  real(sgl), dimension(0:numenin,numlres,numjres) :: GGjlurr ! average radiative width for j,l value
  real(sgl), dimension(0:numenin,numlres,numjres) :: GNjlurr ! average reduced neutron width for j,l value
  real(sgl), dimension(0:numenin,numlres,numjres) :: GXjlurr ! average fission width for j,l value
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for angular distribution in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                                       :: limang ! smallest angle for charged-particle elastic scattering
  real(sgl), dimension(0:numenin,0:numang)                      :: cpang  ! differential cross section
  real(sgl), dimension(0:numpar,0:numlevin,0:numenang,0:numang) :: fang   ! scattering angular distribution
  real(sgl), dimension(0:numenin,0:numang)                      :: fcpang ! scattering angular dist. for charged-particle elastic
  real(sgl), dimension(0:numenin,0:numang)                      :: elasni ! nuclear + interference term
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for spectra in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numpar,numenspec,numddx)            :: ncumddx    ! number of emission energies for DDX spectra
  integer, dimension(0:numpar,numenspec)                   :: ncumout    ! number of emission energies for total production spectra
  integer                                                  :: Nddx       ! number of angles
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: buratio    ! break-up ratio
  real(sgl), dimension(0:numpar,numenspec,numddx,0:numen2) :: ddxemis    ! double-differential emission spectra
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: Eocum      ! emission energies for total production spectra
  real(sgl), dimension(0:numpar,numenspec,numddx,0:numen2) :: Eoddx      ! emission energies for double-differential sp
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: preeqratio ! pre-equilibrium ratio
  real(sgl), dimension(numddx)                             :: rmuddx     ! cosine of angle
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: xsemis     ! total production emission spectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for yields in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numpar,numenspec)                   :: nbegcum  ! first outgoing energy
  integer, dimension(0:numpar,numenspec)                   :: nendcum  ! last outgoing energy
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: Ehistcum ! histogram emission energy for total production spec
  real(sgl), dimension(0:numpar,numenspec,0:numen2)        :: f0cum    ! energy distribution for cumulative spectra
  real(sgl), dimension(0:numpar,numenspec,numddx,0:numen2) :: f0ddx    ! energy distribution for DDX
  real(sgl), dimension(0:numpar,numenin)                   :: xsprod   ! particle production cross section
  real(sgl), dimension(0:numpar,numenin)                   :: yieldany ! yield for (n,anything) channel
  real(sgl), dimension(0:numpar,numenin)                   :: yieldp   ! particle production yield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for residual cross sections in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ,0:numN,0:numlevels)            :: isorpexist  ! flag for existence of isomer of residual nucleus
  integer, dimension(0:numZ,0:numN,numenspec)              :: nbegcumrec  ! first outgoing energy
  integer, dimension(0:numZ,0:numN,numenspec)              :: nendcumrec  ! last outgoing energy
  integer, dimension(0:numZ,0:numN)                        :: Nisorp      ! number of isomers
  integer, dimension(0:numZ,0:numN,numenspec)              :: noutrecrp   ! number of recoil energies for residual nucleus
  logical, dimension(0:numZ,0:numN)                        :: rpexist     ! flag for existence of residual nuclide
  real(sgl), dimension(0:numZ,0:numN,numenin,0:numenrec)   :: Erecrp      ! recoil energy of residual nucleus
  real(sgl), dimension(0:numZ,0:numN,numenspec,0:numenrec) :: Ehistcumrec ! histogram emission energy for recoil spectr
  real(sgl), dimension(0:numZ,0:numN,0:numlevels)          :: Erpiso      ! energy of isomer
  real(sgl), dimension(0:numZ,0:numN)                      :: Ethrp       ! threshold energy for residual product
  real(sgl), dimension(0:numZ,0:numN,0:numlevels)          :: Ethrpiso    ! threshold energy for isomer of residual product
  real(sgl), dimension(0:numZ,0:numN,numenspec,0:numenrec) :: f0cumrec    ! energy distribution for recoil spectra
  real(sgl), dimension(0:numZ,0:numN)                      :: nucmass     ! mass of nucleus
  real(sgl), dimension(0:numZ,0:numN)                      :: Qrp         ! Q-value for residual product
  real(sgl), dimension(0:numZ,0:numN,0:numlevels)          :: Qrpiso      ! Q-value for isomer of residual product
  real(sgl), dimension(0:numZ,0:numN,numenin,0:numenrec)   :: recrp       ! recoils or residual nucleus
  real(sgl), dimension(0:numZ,0:numN,numenin)              :: xsrp        ! residual production cross section
  real(sgl), dimension(0:numZ,0:numN,0:numlevels,numenin)  :: xsrpiso     ! residual production cross section for isomer
  real(sgl), dimension(0:numZ,0:numN,numenin)              :: Yrp         ! residual production yield
  real(sgl), dimension(0:numZ,0:numN,0:numlevels,numenin)  :: Yrpiso      ! residual production yield for isomer
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for photon production in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numpar,0:numlevels,0:numlevels)   :: branchlevel  ! level to which branching takes place
  integer, dimension(0:numpar,0:numlevels)               :: Nbranch      ! number of branches for level
  integer, dimension(0:numchan)                          :: Ngam         ! number of gamma transitions for this nucleus
  integer, dimension(0:numpar,0:numlevels)               :: Ngamdis      ! number of gamma ray lines per level
  integer, dimension(0:numpar)                           :: nlev         ! number of excited levels for nucleus
  real(sgl), dimension(0:numpar,0:numlevels,0:numlevels) :: branchratio  ! branch ratio
  real(sgl), dimension(0:numpar,0:numlevels,numgam)      :: Egamdis      ! energy of gamma ray
  real(sgl), dimension(0:numchan,numgam)                 :: Egamma       ! gamma energy
  real(sgl), dimension(0:numchan,numgam)                 :: Estart       ! starting level
  real(sgl), dimension(0:numchan,numgam,numengam)        :: xsgam        ! exclusive discrete gamma-ray cross section
  real(sgl), dimension(0:numchan,numengam)               :: xsgamcont    ! gamma production cross section per residual nucleus
  real(sgl), dimension(0:numchan,numengam)               :: xsgamdisctot ! total discrete gamma-ray cross section
  real(sgl), dimension(0:numchan,numengam)               :: xsgamtot     ! total (discrete+continuum) gamma-ray cross section
  real(sgl), dimension(0:numchan,numengam)               :: yieldcont    ! continuum gamma-ray yield
  real(sgl), dimension(0:numchan,numgam,numengam)        :: yielddisc    ! discrete gamma-ray yield
  real(sgl), dimension(0:numpar,0:numlevels)             :: yieldg       ! total discrete gamma yield per level
  real(sgl), dimension(0:numchan,numengam)               :: yieldgam     ! gamma yield
  real(sgl), dimension(0:numpar,0:numlevels,numgam)      :: yieldratio   ! yield ratio for level
  real(sgl), dimension(0:numchan,numengam)               :: yieldtot     ! total gamma-ray yield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variables for covariances in ENDF format
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                                       :: flagcovleg    ! flag for covariances for Legendre coefficients
  logical, dimension(numchan)                                   :: flagMTint     ! flag for channel with inter-MT covariance data
  integer, dimension(numchan)                                   :: Arpcov        ! mass number of residual product
  integer, dimension(numencov)                                  :: Ecovindex     ! energy index for main energy grid
  integer, dimension(numchan)                                   :: isorpcov      ! isomer number of residual product
  integer, dimension(numchan)                                   :: MTindex       ! index for MT number
  integer, dimension(numchan)                                   :: MTindexiso    ! index for isomer of MT number
  integer, dimension(numchan)                                   :: MTintindex    ! index for MT number with inter-MT covariance data
  integer, dimension(numchan)                                   :: MTintindexiso ! index for MT number with inter-MT covariance data
  integer                                                       :: Nchancov      ! number of channels with covariance data
  integer                                                       :: Nchancovint   ! number of channels with inter-MT covariance data
  integer                                                       :: Nchanleg      ! number of energies with Legendre covariance data
  integer                                                       :: Ncovrp        ! neutron number of residual product
  integer                                                       :: Nencov        ! number of energies with covariance data
  integer, dimension(nummt)                                     :: Nisocov       ! number of isomers with covariances per MT number
  integer                                                       :: Nleg34        ! number of Legendre coeff. with covariance data
  integer, dimension(numchan)                                   :: Zrpcov        ! charge number of residual product
  real(sgl), dimension(0:numencov+1)                            :: Ecov          ! energy grid for covariances
  real(sgl), dimension(0:numencov)                              :: Eleg          ! energy for Legendre covariance data
  real(sgl), dimension(numchancov,numencov,numchancov,numencov) :: Rcov          ! relative covariance matrix
  real(sgl), dimension(numchan,numencov)                        :: relerr        ! relative cross section uncertainty
  real(sgl), dimension(0:numencov,0:numleg,0:numencov,0:numleg) :: Rleg          ! covariance matrix element for Legend
  real(sgl), dimension(numchan,numencov,numencov)               :: Rmt           ! relative covariance matrix within same MT number
  real(sgl), dimension(numchan,numencov,numencov)               :: Rmtres        ! relative covariance matrix within same MT number
  real(sgl), dimension(numchan,numencov,numencov)               :: Rrp           ! covariance element for residual cross sections
  real(sgl), dimension(numchan,numencov)                        :: xserr         ! cross section uncertainty
end module A0_tefal_mod
