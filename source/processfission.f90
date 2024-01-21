subroutine processfission
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process fission data
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
!   numN           ! maximal number of neutrons from initial compound nucleus
!   numZ           ! maximal number of protons from initial compound nucleus
! Variables for ENDF limits, switches and tolerances
!   Eswitch        ! energy where ENDF - 6 representation is switched (in MeV)
! Variables for input of ENDF structure
!   flagsubfis     ! flag to include subactinide fission
! Variables for info from TALYS
!   Atarget        ! mass number of nucleus
!   eninc          ! incident energy
!   numinc         ! number of incident energies
! Constants
!   fislim         ! mass above which nuclide fissions
!   xsepslow       ! lower limit for cross sections in millibarns
! Variables for total cross sections in ENDF format
!   xsnonel        ! nonelastic cross section
! Variables for residual production cross sections in ENDF format
!   xsrp           ! residual production cross section
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
!   xsexcl         ! exclusive cross section
!
! *** Declaration of local data
!
  implicit none
  logical   :: limit              ! help variable
  integer   :: id                 ! counter for deuterons
  integer   :: nin                ! counter for incident energy
  integer   :: Nix                ! neutron number index for residual nucleus
  integer   :: Zix                ! charge number index for residual nucleus
  real(sgl) :: deltaF             ! help variable
  real(sgl) :: factor             ! multiplication factor
  real(sgl) :: xsrpsum            ! sum over residual production cross sections
!
! ************* Determine limits for fission cross sections ************
!
  limit = .false.
  do nin = 1, numinc
    if (xsexcl(-1, nin) > 0.) then
      xsexcl(-1, nin) = max(xsepslow, xsexcl(-1, nin))
      limit = .true.
    endif
  enddo
!
! For the moment, we set fission false for subactinides
!
  if (Atarget <= fislim .and. .not. flagsubfis) limit = .false.
  if (limit) then
    do id = -2, -1
      do nin = 1, numinc
        xsexcl(id, nin) = max(xsexcl(id, nin), xsepslow)
      enddo
    enddo
  else
    do id = -5, -1
      do nin = 1, numinc
        xsexcl(id, nin) = 0.
      enddo
    enddo
    flagfission = .false.
  endif
!
! Due to numerical instabilities, the sum fission + residual production cross section is sometimes slightly larger than
! the non-elastic cross section. This is normalized.
!
  do nin = 1, numinc
    if (eninc(nin) <= Eswitch) cycle
    xsrpsum = 0.
    do Zix = 0, numZ
      do Nix = 0, numN
        xsrpsum = xsrpsum + xsrp(Zix, Nix, nin)
    enddo
  enddo
    deltaF = xsrpsum + xsexcl(-1, nin) - real(xsnonel(nin))
    if (deltaF > 0.) then
      factor = real(xsnonel(nin)) / (xsrpsum + xsexcl(-1, nin))
      do id = -5, -1
        xsexcl(id, nin) = factor * xsexcl(id, nin)
      enddo
    endif
  enddo
  return
end subroutine processfission
! Copyright A.J. Koning 2021
