subroutine process
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process data from TALYS
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
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flaggpf        ! flag for general purpose library
! Variables for TALYS info
!   k0             ! index of incident particle
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
!
! *** Declaration of local data
!
  implicit none
!
! ************************ Process data from TALYS *********************
!
! processchannels : subroutine to process exclusive channel data
! processresidual : subroutine to process residual production data
! processfission  : subroutine to process fission data
! processdiscrete : subroutine to process discrete level and continuum data
! processyields   : subroutine to process yields and spectra
! processphoton   : subroutine to process gamma data
! processangle    : subroutine to process discrete angular data from TALYS
! processcpelastic: subroutine to process charged-particle elastic scattering data
! processtotal    : subroutine to process total cross sections
!
  if (flagendfdet) call processchannels
  call processresidual
  if (flagfission) call processfission
  call processdiscrete
  if (flaggpf) then
    call processyields
    if (flagendfdet) call processphoton
    if (k0 <= 1) then
      call processangle
    else
      call processcpelastic
    endif
  endif
  call processtotal
  return
end subroutine process
! Copyright A.J. Koning 2021
