subroutine write9_10
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF9 and MF10
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
!   numenin    ! number of incident energies
!   numint     ! number of interpolation sections
!   nummt      ! number of MT numbers
! Variables for initialization of ENDF variables
!   AWR        ! standard mass parameter
!   blank2     ! blank string
!   mtexist    ! flag for existence of MT - number
!   FEND       ! ENDF - 6 format
!   LIS        ! discrete state number
!   MAT        ! MAT number
!   SEND       ! ENDF - 6 format
!   ZA         ! standard charge parameter
! Variables for MF3
!   QM         ! Q - value (in ENDF - 6 format)
! Variables for MF8_10
!   E10        ! incident energy (in ENDF - 6 format)
!   E10ZA      ! incident energy (in ENDF - 6 format)
!   INTER10    ! interpolation scheme
!   INTERZA    ! interpolation scheme
!   IZAP       ! second IZAP - number
!   LFS        ! final state number
!   LFSZA      ! final state number
!   NBT10      ! separation value for interpolation scheme
!   NBTZA      ! separation value for interpolation scheme
!   NP10       ! number of incident energies
!   NPZA       ! number of incident energies
!   NR10       ! number of interpolation ranges
!   NRZA       ! number of interpolation ranges
!   NSt        ! number of final states
!   NZA        ! number of nuclides
!   QIiso      ! Q - value for isomer (in ENDF - 6 format)
!   QIZA       ! Q - value (in ENDF - 6 format)
!   QMZA       ! Q - value (in ENDF - 6 format)
!   xsiso      ! cross section for isomer (in ENDF - 6 format)
!   xsrpZA     ! cross section for residual production (in ENDF - 6 format)
!   ZAPi       ! designation of final nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: ii                        ! counter
  integer   :: iso                       ! counter for isomer
  integer   :: iza                       ! counter for Z,A combinations
  integer   :: MF                        ! MF-number
  integer   :: MT                        ! MT-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(2*numenin)              ! help variable
!
! ***************************** Write MF10 *****************************
!
! hrwrite: subroutine to write header with real values
!
  do MF = 9, 10
    NS = 0
    if (MF == 9) then
      open (unit = 2, file = 'MF9', status = 'replace')
    else
      open (unit = 2, file = 'MF10', status = 'replace')
    endif
    do MT = 1, nummt
      if ( .not. mtexist(MF, MT)) cycle
!
! 1. Residual production cross sections
!
      if (MT == 5) then
        call hrwrite(ZA, AWR, LIS, 0, NZA, 0, MAT, MF, MT, NS)
        do iza = 1, NZA
          call hrwrite(QMZA(iza), QIZA(iza), IZAP(iza), LFSZA(iza), NRZA(iza), NPZA(iza), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
! kwrite    : subroutine to write integer value block
!
          do i = 1, NRZA(iza)
            ii = 2 * i - 1
            Nval(ii) = NBTZA(iza, i)
            Nval(ii + 1) = INTERZA(iza, i)
          enddo
          N = 2 * NRZA(iza)
          call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write cross sections
!
! xwrite: subroutine to write real value block
!
          do i = 1, NPZA(iza)
            ii = 2 * i - 1
            x(ii) = E10ZA(iza, i)
            x(ii + 1) = xsrpZA(iza, i)
          enddo
          N = 2 * NPZA(iza)
          call xwrite(N, x, MAT, MF, MT, NS)
        enddo
      else
!
! 2. Isomeric cross sections
!
        call hrwrite(ZA, AWR, LIS, 0, NSt(MT), 0, MAT, MF, MT, NS)
        do iso = 1, NSt(MT)
          call hrwrite(QM(MT), QIiso(MT, iso), ZAPi(MT), LFS(MT, iso), NR10(MT, iso), NP10(MT, iso), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
! kwrite    : subroutine to write integer value block
!
          do i = 1, NR10(MT, iso)
            ii = 2 * i - 1
            Nval(ii) = NBT10(MT, iso, i)
            Nval(ii + 1) = INTER10(MT, iso, i)
          enddo
          N = 2 * NR10(MT, iso)
          call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write cross sections
!
! xwrite: subroutine to write real value block
!
          do i = 1, NP10(MT, iso)
            ii = 2 * i - 1
            x(ii) = E10(MT, iso, i)
            x(ii + 1) = xsiso(MT, iso, i)
          enddo
          N = 2 * NP10(MT, iso)
          call xwrite(N, x, MAT, MF, MT, NS)
        enddo
      endif
      write(2, fmt = SEND) blank2, MAT, MF, 0, NS
    enddo
    write(2, fmt = FEND) blank2, MAT, 0, 0, NS
    close (unit = 2)
  enddo
  return
end subroutine write9_10
! Copyright A.J. Koning 2021
