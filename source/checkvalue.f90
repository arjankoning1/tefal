subroutine checkvalue
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in values
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
! All global variables
!   numelem         ! number of elements
!   numenin         ! number of incident energies
!   numlevin        ! number of discrete levels
!   nummf           ! number of MF numbers
!   nummt           ! number of MT numbers
!   numrp           ! nunber of residual products
! Variables for input of ENDF structure
!   flagrecoil      ! flag to include recoil information
!   nomf            ! flag to exclude an entire MF
! Variables for input of specific ENDF data
!   adopt           ! logical of existence of MF information (per MT)
!   Eahigh          ! upper energy of MT values to be adopted
!   Ealow           ! lower energy of MT values to be adopted
!   lssfinp         ! 0: URR cross section from MF2, 1: URR cross section
!   urrcomp         ! mode for competition in the URR, 0: none, 1:MT4, 2:all
!   urrenergy       ! upper energy of the URR in MeV
!   urrmode         ! 0: no URR, 1: URR from TALYS, 2: URR from data library
! Variables for input of ENDF library type
!   flageaf         ! flag for EAF - formatted activation library
!   flaggpf         ! flag for general purpose library
! Variables for ENDF limits, switches and tolerances
!   cuteps          ! energy shift at MT5 cutoff energy (in eV)
!   disclim         ! limit for specific MT numbers for discrete levels
!   Eswitch         ! energy where ENDF - 6 representation is switched (in MeV)
!   Eswitch4        ! energy where MF4 representation is switched (in MeV)
!   maxrp           ! maximum number of residual products
!   NMTmax          ! maximum number of MT numbers
! Variables for ENDF covariance input
!   covdiscrete     ! number of disc. inelastic levels with covariances
! Variables for info from TALYS
!   flagtalysdet    ! flag for detailed ENDF - 6 information from TALYS
!   flagtalysrec    ! flag for recoil information from TALYS
!   flagtalysurr    ! flag for URR information from TALYS
!   k0              ! index of incident particle
!   Ltarget         ! excited level of target
!   NLmax           ! maximum number of included discrete levels
!   numinc          ! number of incident energies
!   Ztarget         ! charge number of nucleus
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   range_real_error       ! Test if real variable is out of range
!
! *** Declaration of local data
!
  implicit none
  integer :: imf              ! MF counter
  integer :: imt              ! MT counter
!
! ******************* Check for input variables from TALYS *************
!
  call range_integer_error('k0', k0, 0, 6)
  call range_integer_error('numinc from TALYS', numinc, 0, numenin)
  call range_integer_error('Ztarget', Ztarget, 3, numelem)
  call range_integer_error('nlevmax', NLmax, 0, numlevin)
  call range_integer_error('Ltarget', Ltarget, 0, numlevin)
!
! All parameters need to fall within certain ranges.
! These ranges are specified in this subroutine and in the manual.
!
! ******************* Check for wrong input variables ******************
!
  do imf = 1, nummf
    do imt = 1, nummt
      call range_real_error('Elow', Ealow(imf, imt), 1.e-5, 2.e8, unit = 'eV', index1 = imf, name1 = 'MF', &
 &      index2 = imt, name2 = 'MT')
      call range_real_error('Ehigh', Eahigh(imf, imt), 1.e-5, 2.e8, unit = 'eV', index1 = imf, name1 = 'MF', &
 &      index2 = imt, name2 = 'MT')
    enddo
  enddo
  call range_real_error('Eswitch', Eswitch, 0., 300., unit = 'MeV')
  call range_real_error('Eswitch4', Eswitch4, 0., 300., unit = 'MeV')
  call range_real_error('cuteps', cuteps, 0., 100., unit = 'eV')
  call range_real_error('disclim', disclim, 0., 10000., unit = 'mb')
  call range_real_error('urrenergy', urrenergy, 0., 1., default = -1., unit = 'MeV')
  call range_integer_error('covdiscrete', covdiscrete, 0, 40)
  call range_integer_error('NMTmax', NMTmax, 3, 1000)
  call range_integer_error('urrmode', urrmode, 0, 2)
  call range_integer_error('urrcomp', urrcomp, 0, 2)
  call range_integer_error('lssf', lssfinp, 0, 1, default = -1)
  call range_integer_error('maxrp', maxrp, 0, numrp)
!
! ***************** Check for conflicting keywords *********************
!
! Transport data are not possible in the European Activation File format, only in the ENDF-6 format.
!
  if (flaggpf .and. flageaf) then
    write(*, '(" TEFAL-error: Transport file must be in ", "ENDF-6 format, not in EAF format")')
    stop
  endif
!
! Detailed ENDF-6 information from TALYS should be available if Eswitch > 0.
!
  if ( .not. flagtalysdet .and. (Eswitch > 0..or.Eswitch4 > 0.)) then
    write(*, '(" TEFAL-error: Detailed ENDF-6 information ", "from TALYS should be available if Eswitch > 0.")')
    stop
  endif
!
! Check for presence of recoil data from TALYS
!
  if (flagrecoil .and. .not. flagtalysrec) then
    write(*, '(" TEFAL-error: TALYS recoil information must ", "be available for recoil y. Rerun TALYS with recoil y or ", &
 &    "set recoil n in TEFAL")')
    stop
  endif
!
! Check for presence of URR data from TALYS
!
  if (urrmode == 1 .and. .not. flagtalysurr) then
    write(*, '(" TEFAL-error: TALYS URR information must ", "be available for urrmode 1. Rerun TALYS with urr y or ", &
 &    "set urrmode equal to 0 or 2")')
    stop
  endif
!
! Exclude MF32 if MF33 for resonance channels are given
!
  if (adopt(33,1)) nomf(32) = .true.
  if (adopt(33,2)) nomf(32) = .true.
  if (adopt(33,18)) nomf(32) = .true.
  if (adopt(33,102)) nomf(32) = .true.
  return
end subroutine checkvalue
! Copyright A.J. Koning 2021
