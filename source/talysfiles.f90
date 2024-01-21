subroutine talysfiles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read data files from TALYS
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
! Variables for input of specific ENDF data
!   urrmode        ! 0: no URR, 1: URR from TALYS, 2: URR from data library
! Variables for ENDF covariance input
!   flagcovar      ! flag for covariances
! Variables for TALYS info
!   k0             ! index of incident particle
!
! *** Declaration of local data
!
  implicit none
!
! talyscovar    : subroutine to read covariance data from TASMAN
! talystotal    : subroutine to read total cross sections from TALYS
! talyschannels : subroutine to read exclusive channel data from TALYS
! talysresidual : subroutine to read residual production data from TALYS
! talysfission  : subroutine to read fission data from TALYS
! talysyields   : subroutine to read yields from TALYS
! talysspectra  : subroutine to read spectra and pre-equilibrium ratios from TALYS
! talysdiscrete : subroutine to read discrete level and continuum data from TALYS
! talysurr      : subroutine to read URR data from TALYS
! talysphoton   : subroutine to read gamma decay scheme from TALYS
! talyscpelastic: subroutine to read charged-particle elastic scattering data from TALYS
!
! Exclusive reaction channels, discrete level and photon production data per reaction channel are only stored for incident neutrons.
!
  if (flagcovar) then
    call talyscovar
  endif
  call talystotal
  if (flagendfdet .or. flageaf) call talyschannels
  call talysresidual
  call talysfission
  if (flaggpf) then
    call talysyields
    call talysspectra
    if (flagendfdet) then
      call talysdiscrete
      if (urrmode == 1) then
        call talysurr
      endif
      call talysphoton
    endif
    if (k0 > 1) call talyscpelastic
  endif
  return
end subroutine talysfiles
! Copyright A.J. Koning 2021
