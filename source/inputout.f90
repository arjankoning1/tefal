subroutine inputout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write input parameters
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
!   nummf           ! number of MF numbers
!   nummt           ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt           ! logical for existence of MF information (per MT)
!   adoptfile       ! name of library for MF information (per MT)
!   background      ! file with background cross sections
!   Eahigh          ! upper energy of MT values to be adopted
!   Ealow           ! lower energy of MT values to be adopted
!   lssfinp         ! 0: URR cross section from MF2, 1: URR cross section
!   urrcomp         ! mode for competition in the URR, 0: none, 1:MT4, 2:all
!   urrenergy       ! upper energy of the URR in MeV
!   urrmode         ! 0: no URR, 1: URR from TALYS, 2: URR from data library
! Variables for input of ENDF library type
!   flagendfdet     ! flag for detailed ENDF - 6 information per channel
!   endffile        ! name of ENDF file
!   flageaf         ! flag for EAF - formatted activation library
!   flagclean       ! flag to clean up double points
!   flaggpf         ! flag for general purpose library
!   flaghigh        ! flag for high energies ( > 20 MeV)
! Variables for ENDF limits, switches and tolerances
!   cuteps          ! energy shift at MT5 cutoff energy (in eV)
!   disclim         ! limit for specific MT numbers for discrete levels
!   Eswitch         ! energy where ENDF - 6 representation is switched (in MeV)
!   Eswitch4        ! energy where MF4 representation is switched (in MeV)
!   NMTmax          ! maximum number of MT numbers
! Variables for input of ENDF structure
!   flagaddlow      ! flag to add low - energy Q>0 reactions to total cross section
!   flagcapt6       ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagdisc6       ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
!   flagfis10       ! flag to put (subactinide) fission cross sections in MF10
!   flaggam13       ! flag to use MF13 for gamma prod. instead of MF12 (if not in MF6)
!   flaggamdisc     ! flag to store gamma prod. per level (MT51..) instead of MT4
!   flaggamspec     ! flag to store gamma prod. only as a spectrum and not per level
!   flagmtall       ! flag to include all defined MT numbers from ENDF manual
!   flagmtextra     ! flag to include extra MT numbers up to MT200
!   flagmulti       ! flag to include multi - chance fission
!   flagngn         ! flag to include (n,gamma n) data
!   flagpara        ! flag to include partial cross sections for alphs
!   flagpard        ! flag to include partial cross sections for deuterons
!   flagparh        ! flag to include partial cross sections for helions
!   flagparp        ! flag to include partial cross sections for protons
!   flagpart        ! flag to include partial cross sections for tritons
!   flagpart6       ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
!   flagrecoil      ! flag to include recoil information
!   flagrenorm      ! flag for renormalization of spectra
!   flagres         ! flag to include resonance parameters
!   flagrp10        ! flag to put residual production cross sections in MF10
!   flagrp6         ! flag to put residual production cross sections in MF6
!   flagsubfis      ! flag to include subactinide fission
!   flagtabddx      ! flag to give explicit DDX in MF6
! Variables for input of ENDF MF1
!   author          ! author
!   endftext        ! file with MF1 information
!   identifier      ! library identifier
!   lab             ! laboratory
! Variables for ENDF covariance input
!   covdiscrete     ! number of disc. inelastic levels with covariances
!   flagcovar       ! flag for covariances
!   flagintercor    ! flag for inter - MT covariance data
!   flagparcov      ! flag to include covariances for MT600 - 849
! Variables for initialization of ENDF format
!   rec0            ! first line of ENDF - 6 file
! Variables for reading TEFAL input lines
!   inline          ! input line
!   nlines          ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  character(len=1) :: yesno     ! function to assign y or n to logical value
  integer          :: i         ! counter
  integer          :: imf       ! MF counter
  integer          :: imt       ! MT counter
!
! ******************** Info about created data library *****************
!
  write(*, '(/" ########### DATA FILE ##########")')
  write(*, '(/, a66)') rec0(1:66)
!
! ************************** User input file ***************************
!
  write(*, '(/" ########## USER INPUT ##########")')
  write(*, '(/" USER INPUT FILE"/)')
  do i = 1, nlines
    write(*, '(1x, a)') trim(inline(i))
  enddo
!
! ********* All possible input parameters including defaults ***********
!
  write(*, '(/" USER INPUT FILE + DEFAULTS"/)')
  write(*, '(" Keyword           Value   Variable     Explanation"/)')
!
! 1. Type of data library
!
! yesno    : function to assign y or n to logical value
!
  write(*, '(" #"/" # Type of data library"/" #")')
  write(*, '(" gpf                 ", a1, "     flaggpf       flag for general purpose library")') yesno(flaggpf)
  write(*, '(" eaf                 ", a1, "     flageaf       flag for EAF-formatted activation library")') yesno(flageaf)
!
! 2. Specific structure of the data file
!
  write(*, '(" #"/" # Specific structure of the data file")')
  write(*, '(" #")')
  write(*, '(" linenumbers         ", a1, "     flaglinenum   flag to write line numbers of ENDF file")') yesno(flaglinenum)
  write(*, '(" high                ", a1, "     flaghigh      flag for high energies")') yesno(flaghigh)
  write(*, '(" mtall               ", a1, "     flagmtall     flag to include all defined MT numbers from ENDF manual")') &
 &  yesno(flagmtall)
  write(*, '(" mtextra             ", a1, "     flagmtextra   flag to include extra numbers up to MT200")') yesno(flagmtextra)
  write(*, '(" Eswitch      ", f8.3, "     Eswitch       energy where ENDF-6 representation is switched (in MeV)")') Eswitch
  write(*, '(" Eswitch4     ", f8.3, "     Eswitch4      energy where MF4 representation is switched (in MeV)")') Eswitch4
  write(*, '(" cuteps       ", f8.3, "     cuteps        energy shift  at MT5 cutoff energy (in eV)")') cuteps
  write(*, '(" resonance           ", a1, "     flagres       flag to include resonance parameters")') yesno(flagres)
  write(*, '(" capt6               ", a1, "     flagcapt6     flag to put all MT102 gamma production ", &
 &  "in MF6 instead of MF12/14/15")') yesno(flagcapt6)
  write(*, '(" disc6               ", a1, "     flagdisc6     flag to put all discrete angular distributions and ", &
 &  "gamma production in MF6 instead of MF4/12/14")') yesno(flagdisc6)
  write(*, '(" part6               ", a1, "     flagpart6     flag to put all gamma production for partial cross ", &
 &  "sections in MF6 instead of MF12/14/15")') yesno(flagpart6)
  write(*, '(" rp6                 ", a1, "     flagrp6       flag to put all residual production cross ", &
 &   "sections in MF6")') yesno(flagrp6)
  write(*, '(" rp10                ", a1, "     flagrp10      flag to put all residual production cross ", &
 &   "sections in MF10")') yesno(flagrp10)
  write(*, '(" fis10               ", a1, "     flagfis10     flag to put (subactinide) fission cross ", &
 &   "sections in MF10")') yesno(flagfis10)
  write(*, '(" gam13               ", a1, "     flaggam13     flag to use MF13 for gamma production instead of ", &
 &   "MF12 (if not in MF6)")') yesno(flaggam13)
  write(*, '(" gamdisc             ", a1, "     flaggamdisc   flag to store gamma production per discrete level ", &
 &   "cross section (MT51..) instead of in the total of MT4")') yesno(flaggamdisc)
  write(*, '(" gamspectrum         ", a1, "     flaggamspec   flag to store gamma production only as a spectrum ", &
 &   "and not per discrete level")') yesno(flaggamspec)
  write(*, '(" ngn                 ", a1, "     flagngn       flag to include (n,gamma n) data")') yesno(flagngn)
  write(*, '(" partialn            ", a1, "     flagparn      flag to include partial cross sections for neutrons")') &
 &  yesno(flagparn)
  write(*, '(" partialp            ", a1, "     flagparp      flag to include partial cross sections for protons")') &
 &  yesno(flagparp)
  write(*, '(" partiald            ", a1, "     flagpard      flag to include partial cross sections for deuterons")') &
 &  yesno(flagpard)
  write(*, '(" partialt            ", a1, "     flagpart      flag to include partial cross sections for tritons")') &
 &  yesno(flagpart)
  write(*, '(" partialh            ", a1, "     flagparh      flag to include partial cross sections for helions")') &
 &  yesno(flagparh)
  write(*, '(" partiala            ", a1, "     flagpara      flag to include partial cross sections for alphas")') &
 &  yesno(flagpara)
  write(*, '(" endfdetail          ", a1, "     flagendfdet   flag for detailed ENDF-6 information per channel")') &
 &  yesno(flagendfdet)
  write(*, '(" addlow              ", a1, "     flagaddlow    flag to add low-energy Q>0 reactions to total", &
 &  " cross section")') yesno(flagaddlow)
  write(*, '(" clean               ", a1, "     flagclean     flag for cleaning up double points")') yesno(flagclean)
  write(*, '(" renorm              ", a1, "     flagrenorm    flag for renormalization of spectra")') yesno(flagrenorm)
  write(*, '(" multichance         ", a1, "     flagmulti     flag to include multi-chance fission")') yesno(flagmulti)
  write(*, '(" recoil              ", a1, "     flagrecoil    flag to include recoil information")') yesno(flagrecoil)
  write(*, '(" subfission          ", a1, "     flagsubfis    flag to include subactinide fission")') yesno(flagsubfis)
  write(*, '(" urrmode             ", i1, "     urrmode       0: no URR, 1: URR from TALYS, 2: URR from data library")') urrmode
  write(*, '(" lssf               ", i2, "     lssfinp       0: URR cross section from MF2, 1: URR cross section", &
 &  " from URR (input option)")') lssfinp
  write(*, '(" urrcomp             ", i1, "     urrcomp       mode for competition in the URR, 0:none, 1:MT4, 2:all")') urrcomp
  write(*, '(" urrenergy        ", f8.5, " urrenergy     upper energy of the URR in MeV")') urrenergy
  write(*, '(" NMTmax               ", i4, " NMTmax        maximum number of MT numbers")') NMTmax
  write(*, '(" disclim             ", f5.2, " disclim       limit for specific MT numbers for discrete levels")')  disclim
  write(*, '(" tabddx              ", a1, "     flagtabddx    flag to give explicit DDX in MF6")') yesno(flagtabddx)
!
! 3. Use of specific data
!
  write(*, '(" #"/" # Use of specific data"/" #")')
  do imf = 1, nummf
    if (adopt(imf, nummt)) then
      write(*, '(" adopt    ", i4, 1x, a, 2es12.5)') imf, trim(adoptfile(imf, nummt)), Ealow(imf, nummt), Eahigh(imf, nummt)
    else
      do imt = 1, nummt
        if (adopt(imf, imt)) write(*, '(" adopt    ", 2i4, 1x, a, 2es12.5)') &
 &        imf, imt, trim(adoptfile(imf, imt)), Ealow(imf, imt), Eahigh(imf, imt)
      enddo
    endif
  enddo
  if (background(1:1) /= ' ') write(*, '(" background ", a)') trim(background)
!
! 4. Information for MF1
!
  write(*, '(" #"/" # Information for MF1"/" #")')
  write(*, '(" author             ", a33)') author
  write(*, '(" lab                ", a11)') lab
  write(*, '(" identifier         ", a10)') identifier
  write(*, '(" endftext           ", a)') trim(endftext)
!
! 5. Covariance data
!
  write(*, '(" #"/" # Covariance data"/" #")')
  write(*, '(" covariance          ", a1, "     flagcovar     flag for covariances")') yesno(flagcovar)
  write(*, '(" partialcov          ", a1, "     flagparcov    flag to include covariances for MT600-849")') yesno(flagparcov)
  write(*, '(" covdiscrete        ", i2, "     covdiscrete   number of discrete inelastic levels with covariances")') &
 &  covdiscrete
  write(*, '(" intercor            ", a1, "     flagintercor  flag for inter-MT covariance data")') yesno(flagintercor)
!
! 6. Miscellaneous
!
  write(*, '(" #"/" # Miscellaneous"/" #")')
  write(*, '(" endffile ", a, " endffile      name of ENDF file")') trim(endffile)
  write( * , '()' )
  return
end subroutine inputout
! Copyright A.J. Koning 2021
