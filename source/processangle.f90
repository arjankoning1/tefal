subroutine processangle
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process discrete angular data from TALYS
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
!   sgl         ! single precision kind
! Variables for input of ENDF library type
!   flageaf     ! flag for EAF - formatted activation library
! Variables for reaction initialization
!   Nenang      ! number of incident energies for ang. dist.
!   nlevmax     ! number of included discrete levels
!   rmu         ! cosine of the angle
! Variables for info from TALYS
!   k0          ! index of incident particle
! Variables for discrete state cross sections in ENDF format
!   cleg0       ! Legendre coefficients
!   jdis        ! spin of level
!   ncleg       ! number of Legendre coefficients
!   xsang       ! differential cross section
! Variables for angular distributions in ENDF format
!   fang        ! scattering angular distribution
!
! *** Declaration of local data
!
  implicit none
  integer   :: deltaj             ! spin difference
  integer   :: iang               ! running variable for angle
  integer   :: j                  ! counter
  integer   :: j0                 ! ground state spin
  integer   :: L                  ! counter for Legendre coefficients
  integer   :: nen                ! energy counter
  integer   :: nex                ! discrete level
  integer   :: nex0               ! base number for discrete level
  integer   :: nexinel            ! leven number
  integer   :: nleg               ! number of Legendre coefficients
  integer   :: type               ! particle type
  real(sgl) :: da                 ! help variable
  real(sgl) :: xsi                ! help variable
!
! ********** Legendre coefficients and angular distributions ***********
!
! The highest order of the Legendre coefficients is determined.
! In the ENDF-6 format, the Legendre coefficients are defined relative to the first Legendre coefficient,
! which is always equal to 1.
!
  do type = k0, k0
    do nex = 0, nlevmax
Loop1: do nen = 1, Nenang
        nleg = ncleg(type, nex, nen)
        do L = 0, nleg + 1
          if (mod(L, 2) == 0 .and. abs(cleg0(type, nex, nen, L)) < 1.e-10) then
            ncleg(type, nex, nen) = max(L - 2, 2)
            cycle Loop1
          endif
        enddo
      enddo Loop1
    enddo
  enddo
!
! For (n,p),.....(n,a) scattering, the Legendre coefficients of inelastic scattering are adopted.
! For each (n,p)...(n,a) level, we search for inelastic level with the same (or nearest) spin, and
! adopt the corresponding Legendre coefficients.
! This is done with a preference for the weakest coupled (i.e. highest) inelastic levels.
!
  nexinel = 0
  do type = 2, 6
    do nex = 0, nlevmax
      j = int(jdis(type, nex))
      deltaj = 100
      do nex0 = nlevmax, 1, -1
        j0 = int(jdis(k0, nex0))
        if (abs(j - j0) < deltaj) then
          deltaj = abs(j - j0)
          nexinel = nex0
        endif
        if (deltaj == 0) exit
      enddo
      do nen = 1, Nenang
        ncleg(type, nex, nen) = ncleg(k0, nexinel, nen)
        do L = 0, ncleg(type, nex, nen)
          cleg0(type, nex, nen, L) = cleg0(k0, nexinel, nen, L)
        enddo
      enddo
    enddo
  enddo
!
! For high energies: Renormalize angular distributions with cross section
!
  if ( .not. flageaf) then
    do type = 1, 6
      do nex = 0, nlevmax
        do nen = 1, Nenang
          xsi = 0.
          do iang = 0, 89
            da = 0.5 * (rmu(iang) - rmu(iang + 1))
            xsi = xsi + da * (xsang(type, nex, nen, iang) + xsang(type, nex, nen, iang + 1))
          enddo
          if (xsi > 0.) then
            do iang = 0, 90
              fang(type, nex, nen, iang) = xsang(type, nex, nen, iang) / xsi
            enddo
          else
            do iang = 0, 90
              fang(type, nex, nen, iang) = fang(type, nex, nen - 1, iang)
            enddo
          endif
        enddo
      enddo
    enddo
  endif
  return
end subroutine processangle
! Copyright A.J. Koning 2021
