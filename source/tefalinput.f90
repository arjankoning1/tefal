subroutine tefalinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Input
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! talysinfo   : subroutine for basic reaction information from TALYS
! readinput   : subroutine to read user input
! input       : subroutine to read input for keywords
! checkkeyword: subroutine to check for errors in keywords
! checkvalue  : subroutine to check for errors in values
! checkfiles  : subroutine to check for errors in external data files
!
  call talysinfo
  call readinput
  call input_endfstruc
  call input_endfMF1
  call input_endftype
  call input_endflimits
  call input_endfspec
  call input_endfcovar
  call checkkeyword
  call checkvalue
  call checkfiles
  return
end subroutine tefalinput
! Copyright A.J. Koning 2021
