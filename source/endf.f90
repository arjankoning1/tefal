function endf(x)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Create a number in ENDF-6 format
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  character(len=11) :: endf      ! function for a number in ENDF-6 format
  character(len=13) :: string    ! line with parameter value
  real(sgl)         :: x         ! help variable
!
! *********************** ENDF-6 format function ***********************
!
! This function removes the 'E' from the exponential
!
! endf: function for a number in ENDF-6 format
!
  write(string, '(es13.6)') x
  if (x == 0 .or. (abs(x) <= 1.e10 .and. abs(x) >= 1.e-9)) then
    endf = string(1:9) //string(11:11) //string(13:13)
  else
    endf = string(1:8) //string(11:13)
  endif
  return
end function endf
! Copyright A.J. Koning 2021
