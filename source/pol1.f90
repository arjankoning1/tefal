subroutine pol1(x1, x2, y1, y2, x, y)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolation of first order
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
  real(sgl) :: fac            ! factor
  real(sgl) :: x              ! help variable
  real(sgl) :: x1             ! coordinates of intersection points inside the bin
  real(sgl) :: x2             ! coordinates of the 2nd summit of the triangle
  real(sgl) :: y              ! coordinates of the point to test
  real(sgl) :: y1             ! coordinates of the 1st summit of the triangle
  real(sgl) :: y2             ! coordinates of the 2nd summit of the triangle
!
! ***************************** Interpolation **************************
!
  fac = (x-x1)/(x2-x1)
  y = y1 + fac * (y2 - y1)
  return
end subroutine pol1
! Copyright A.J. Koning 2021
