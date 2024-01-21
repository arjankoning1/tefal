subroutine reacinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of nuclear reaction info
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
!   numenspec      ! number of incident energies for spectra
! Variables for ENDF limits, switches and tolerances
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
!   Eswitch4       ! energy where MF4 representation is switched (in MeV)
! Variables for input of ENDF library type
!   flageaf        ! flag for EAF - formatted activation library
!   flaggpf        ! flag for general purpose library
!   flaghigh       ! flag for high energies ( > 20 MeV)
! Variables for input of ENDF structure
!   flagmtextra    ! flag to include extra MT numbers up to MT200
!   flagrecoil     ! flag to include recoil information
! Variables for reaction initialization
!   Eangindex      ! enegy index for angular distributions
!   Egamindex      ! enegy index for gamma cross sections
!   EHres          ! upper energy in resonance range
!   eninccut       ! last incident energy before high - energy format
!   enincmax       ! maximal incident energy
!   Especindex     ! enegy index for spectra
!   includeres     ! flag to include resonance parameters
!   massN          ! mass of nucleus in neutron units
!   Nenang         ! number of incident energies for ang. dist.
!   Nengam         ! number of incident energies for gamna c.s.
!   Nenspec        ! number of incident energies for spectra
!   nlevmax        ! number of included discrete levels
!   nuclid         ! nuclide
!   numcut         ! number of energies before high - energy format
!   numcut4        ! number of energies before MF4 high - energy format
!   relmass        ! mass relative to neutron mass
!   rmu            ! cosine of the angle
! Variables for info from TALYS
!   eninc          ! incident energy
!   NLmax          ! maximum number of included discrete levels
!   numinc         ! number of incident energies
!   tarmass        ! mass of nucleus
!   Ztarget        ! charge number of nucleus
! Constants
!   nuc            ! nuclide
!   parmass        ! mass of particle in a.m.u.
!   pi             ! pi
!
! *** Declaration of local data
!
  implicit none
  integer           :: iang                  ! running variable for angle
  integer           :: ii                    ! counter
  integer           :: ispec                 ! counter for spectra
  integer           :: ispecprev             ! counter for spectra
  integer           :: nen                   ! energy counter
  integer           :: nin                   ! counter for incident energy
  integer           :: type                  ! particle type
  real(sgl)         :: angle                 ! angle in degrees
  real(sgl)         :: Edist                 ! help variable
  real(sgl)         :: Edisttry              ! help variable
  real(sgl)         :: Espec(0:numenspec)    ! energy of spectrum
!
! ************ Derived variables and overruling of input ***************
!
! "Special" high energy MT format for FISPACT
!
  if (flageaf .and. flaghigh) then
    Eswitch = 60.
    flagmtextra = .true.
  endif
!
! Reset energy for switch from Legendre coeff. to angular distribution
!
  Eswitch4 = min(Eswitch4, Eswitch)
!
! Check for Eswitch energies to be in the energy grid
!
  do nin = 1, numinc - 1
    if (eninc(nin) < Eswitch .and. eninc(nin + 1) > Eswitch) then
      write(*, '(" TEFAL-error: Eswitch=", f10.5, " MeV should be in energy grid")') Eswitch
      stop
    endif
  enddo
  do nin = 1, numinc - 1
    if (eninc(nin) < Eswitch4 .and. eninc(nin + 1) > Eswitch4) then
      write(*, '(" TEFAL-error: Eswitch4=", f10.5, " MeV should be in energy grid")') Eswitch4
      stop
    endif
  enddo
!
! Recoil information is only included in general purpose files
!
  if ( .not. flaggpf) flagrecoil = .false.
  nlevmax = min(NLmax, 40)
!
! ************** Initialization of nuclear reaction info ***************
!
! Masses relative to the neutron mass and standard ENDF-6 quantities are determined.
!
  do type = 0, 6
    relmass(type) = parmass(type) / parmass(1)
  enddo
  nuclid = nuc(Ztarget)
  massN = real(tarmass / parmass(1))
  includeres = .false.
  EHres = 0.
!
! Determine minimal and maximal energy and energy where the transition from the low energy to the high energy transition occurs.
!
  if (k0 == 1) then
    EminMeV = 1.e-11
  else
    EminMeV = 0.5 * eninc(1)
  endif
  EmineV = EminMeV * 1.e6
  numcut = numinc
  do nin = 1, numinc
    if (eninc(nin) > Eswitch) then
      numcut = nin - 1
      exit
    endif
  enddo
  numcut = max(numcut, 1)
  if (numcut == numinc .and. .not. flageaf) flaghigh = .false.
  eninccut = eninc(numcut)
  if (flaghigh) then
    enincmax = eninc(numinc)
  else
    enincmax = eninccut
  endif
  numcut4 = numinc
  do nin = 1, numinc
    if (eninc(nin) > Eswitch4) then
      numcut4 = nin - 1
      exit
    endif
  enddo
!
! Determine grid for angular distributions.
!
  if (flaggpf) then
    do iang = 0, 90
      angle = 2 * iang * pi / 180.
      rmu(iang) = cos(angle)
    enddo
  endif
!
! Set secondary energy grids (spectra, angles, gamma's).
! These grids will be less dense as that of the cross section.
!
  Espec(0) = 0.
  Espec(1) = 1.e-11
  Espec(2) = 2.53e-8
  Espec(3) = 1.e-6
  Espec(4) = 1.e-5
  Espec(5) = 1.e-4
  Espec(6) = 0.001
  Espec(7) = 0.002
  Espec(8) = 0.005
  Espec(9) = 0.01
  Espec(10) = 0.02
  Espec(11) = 0.05
  Espec(12) = 0.1
  Espec(13) = 0.2
  Espec(14) = 0.4
  Espec(15) = 0.6
  Espec(16) = 0.8
  Espec(17) = 1.
  Espec(18) = 1.4
  Espec(19) = 1.8
  Espec(20) = 2.2
  Espec(21) = 2.6
  Espec(22) = 3.
  Espec(23) = 3.4
  Espec(24) = 3.8
  Espec(25) = 4.2
  Espec(26) = 4.6
  Espec(27) = 5.
  Espec(28) = 5.5
  Espec(29) = 6.
  Espec(30) = 6.5
  Espec(31) = 7.
  Espec(32) = 7.5
  Espec(33) = 8.
  Espec(34) = 8.5
  Espec(35) = 9.
  Espec(36) = 10.
  Espec(37) = 12.
  Espec(38) = 14.
  Espec(39) = 16.
  Espec(40) = 18.
  Espec(41) = 20.
  Espec(42) = 22.
  Espec(43) = 24.
  Espec(44) = 26.
  Espec(45) = 28.
  Espec(46) = 30.
  Espec(47) = 34.
  Espec(48) = 38.
  Espec(49) = 42.
  Espec(50) = 46.
  Espec(51) = 50.
  Espec(52) = 55.
  Espec(53) = 60.
  Espec(54) = 70.
  Espec(55) = 80.
  Espec(56) = 100.
  Espec(57) = 120.
  Espec(58) = 140.
  Espec(59) = 160.
  Espec(60) = 200.
  do nen = 1, 20
    Espec(60 + nen) = 200. + 40. * nen
  enddo
!
! Set grid
!
  ii = 0
  ispec = 0
  ispecprev = 0
  do nen = 1, numenspec
    Edist = 1000.
    do nin = 1, numinc
      Edisttry = abs(Espec(nen) - eninc(nin))
      if (Edisttry < Edist) then
        Edist = Edisttry
        ispec = nin
      endif
    enddo
    if (ispec /= ispecprev) then
      ii = ii + 1
      Especindex(ii) = ispec
      ispecprev = ispec
    endif
  enddo
  Nenspec = ii
  Nengam = Nenspec
  Nenang = Nenspec
  do nen = 1, Nenspec
    Egamindex(nen) = Especindex(nen)
    Eangindex(nen) = Especindex(nen)
  enddo
  return
end subroutine reacinitial
! Copyright A.J. Koning 2021
