subroutine processcpelastic
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process charged-particle elastic scattering data
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
!   sgl        ! single precision kind
! All global variables
!   numang     ! number of angles
! Constants
!   pi         ! pi
! Variables for reaction initialization
!   rmu        ! cosine of the angle
! Variables for info from TALYS
!   numinc     ! number of incident energies
! Variables for total cross sections in ENDF format
!   xselas     ! total elastic cross section
! Variables for angular distributions in ENDF format
!   cpang      ! differential cross section
!   fcpang     ! scattering angular distribution for charged particles
!   limang     ! smallest angle for charged - particle elastic scattering
!   elasni     ! nuclear + interference elastic angular distribution
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang                        ! running variable for angle
  integer   :: nin                         ! counter for incident energy
  integer   :: ninprev                     ! previous incident energy
  real(sgl) :: da                          ! help variable
  real(sgl) :: xsel                        ! total elastic cross section
!
! **** Angular distributions for charged-particle elastic scattering ***
!
! The nuclear term is removed from charged-particle elastic scattering.
!
  limang = 3
  ninprev = 1
  do nin = numinc, 1, -1
    xsel = 0.
    do iang = 90, limang + 1, -1
      da = 0.5 * (rmu(iang - 1) - rmu(iang))
      xsel = xsel + da * (elasni(nin, iang) + elasni(nin, iang - 1))
    enddo
    if (nin < numinc .and. xsel == 0.) xsel = real(xselas(nin + 1))
    xselas(nin) = xsel
    if (xsel == xselas(ninprev)) then
      do iang = limang, 90
        fcpang(nin, iang) = fcpang(ninprev, iang)
      enddo
    else
      do iang = limang, 90
        if (xsel /= 0.) then
          fcpang(nin, iang) = elasni(nin, iang) / xsel
        else
          fcpang(nin, iang) = 1. / (2. * pi)
        endif
      enddo
    endif
    ninprev = nin
  enddo
  return
end subroutine processcpelastic
! Copyright A.J. Koning 2021
