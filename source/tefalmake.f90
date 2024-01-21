subroutine tefalmake
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Construction of ENDF-6 data
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
! Variables for input of ENDF library type
!   flageaf        ! flag for EAF - formatted activation library
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
! Variables for input of ENDF structure
!   flagcapt6      ! flag to put MT102 gamma prod. in MF6 instead of MF12 / 14 / 15
!   flagdisc6      ! flag for disc. ang. distr. and gam prod. in MF6 not MF4 / 12 / 14
!   flagfis10      ! flag to put (subactinide) fission cross sections in MF10
!   flaggam13      ! flag to use MF13 for gamma prod. instead of MF12 (if not in MF6)
!   flagpart6      ! flag for gam. prod. for partial c.s. in MF6 not in MF12 / 14 / 15
! Variables for ENDF covariance input
!   flagcovar      ! flag for covariances
! Variables for TALYS info
!   k0             ! index of incident particle
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
! Variables for covariances in ENDF format
!   flagcovleg     ! flag for covariances for Legendre coefficients
!
! *** Declaration of local data
!
  implicit none
!
! *************************** Make MF files ****************************
!
! read1      : subroutine to read MF1 (fission neutrons) from existin ENDF-6 data library
!
  if (flaggpf .and. k0 <= 1) then
    if (flagfission .and. .not. flagfis10) call read1
!
! For general purpose libraries, resonance parameters can be read from existing data libraries.
! The corresponding background cross sections are read from MF3.
!
! make2 : subroutine to make MF2
! write2: subroutine to write MF2
!
! k0=1 means neutron-induced data.
!
    if (k0 == 1) then
      call make2
      call write2
    endif
  endif
!
! Transport and activation libraries (In the EAF-format, all data is stored in MF3).
!
! make3  : subroutine to make MF3
! write3 : subroutine to write MF3
!
  call make3
  call write3
  if (.not. flageaf) then
!
! Neutron angular distributions.
!
! make4      : subroutine to make MF4
!
    if (k0 == 1 .and. flagendfdet .and. flaggpf) call make4
!
! Secondary energy distributions
!
! make5    : subroutine to make MF5 (contains write5)
!
    if (k0 <= 1 .and. flagfission .and. .not. flagfis10) call make5
!
! Yields and secondary distributions
!
! make6: subroutine to make MF6 (contains write6)
!
    call make6
!
! Isomeric cross sections
!
! make9_10 : subroutine to make MF9 and MF10
! make8    : subroutine to make MF8
! write9_10: subroutine to write MF9 and MF10
! write8   : subroutine to write MF8
!
    if (flagendfdet) then
      call make9_10
      call make8
      call write9_10
      call write8
!
! Photon production data
!
! make12   : subroutine to make MF12 (contains write12)
! make13   : subroutine to make MF13 (contains write13)
! make14   : subroutine to make MF14 (contains write14)
! make15   : subroutine to make MF15 (contains write15)
!
      if ( .not. flagcapt6 .or. .not. flagdisc6 .or. .not. flagpart6) then
        call make12
        if (flaggam13) call make13
        call make14
        call make15
      endif
    endif
  endif
!
! Covariance data
!
! make33    : subroutine to make MF33 and MF40
! write33   : subroutine to write MF33 and MF40
! make31    : subroutine to make MF31
! write31   : subroutine to write MF31
! make32    : subroutine to make MF32
! write32   : subroutine to write MF32
! make34    : subroutine to make MF34
! write34   : subroutine to write MF34
! make35    : subroutine to make MF35
! write35   : subroutine to write MF35
!
  if (flagcovar) then
    call make33(33)
    call write33(33)
    if (.not. flageaf) then
      call make33(40)
      call write33(40)
      if (k0 == 1 .and. flagendfdet .and. flaggpf) then
        if (flagfission .and. .not. flagfis10) then
          call read31
          call write31
          call make35
          call write35
        endif
        call read32
        call write32
        if (flagcovleg) then
          call make34
          call write34
        endif
      endif
    endif
  endif
!
! MF1
!
! make1 : subroutine to make MF1
! write1: subroutine to write MF1
!
  if (.not. flageaf) then
    call make1
    call write1
  endif
!
! Write final ENDF-6 file
!
! finalwrite: subroutine to write final ENDF-6 file
! timer     : subroutine for output of execution time
!
  call finalwrite
  call timer
  return
end subroutine tefalmake
! Copyright A.J. Koning 2021
