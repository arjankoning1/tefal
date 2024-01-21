subroutine linear(x, y, N, xlin, ylin, Nlin)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Linearize cross section with exponential shape (1/v)
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
  integer   :: i                     ! counter
  integer   :: j                     ! counter
  integer   :: N                     ! neutron number of residual nucleus
  integer   :: Nlin                  ! number of linearized points
  real(sgl) :: logx                  ! log of x value
  real(sgl) :: logx1                 ! log of x value
  real(sgl) :: logx2                 ! log of x value
  real(sgl) :: logxb                 ! log of begin x value
  real(sgl) :: logxe                 ! log of end x value
  real(sgl) :: logy1                 ! log of y value
  real(sgl) :: logy2                 ! log of y value
  real(sgl) :: x(0:N)                ! help variable
  real(sgl) :: xb                    ! begin x value
  real(sgl) :: xe                    ! end x value
  real(sgl) :: xlin(Nlin)            ! x value
  real(sgl) :: y(N)                  ! coordinates of the point to test
  real(sgl) :: ylin(Nlin)            ! y value
!
! ******************************* Linearization ************************
!
  xb = x(1)
  xe = x(N)
  xlin(1) = xb
  ylin(1) = y(1)
  if (N == 1) return
  xlin(Nlin) = xe
  ylin(Nlin) = y(N)
  logxb = log(xb)
  logxe = log(xe)
  do i = 2, Nlin - 1
    xlin(i) = exp(logxb + real(i) / Nlin * (logxe - logxb))
  enddo
  do j = 2, N - 1
    do i = 1, Nlin - 1
      if (x(j) >= xlin(i) .and. x(j) < xlin(i + 1)) then
        if (abs(xlin(i) - x(j)) < abs(xlin(i + 1) - x(j))) then
          xlin(i) = x(j)
        else
          xlin(i + 1) = x(j)
        endif
      endif
    enddo
  enddo
  do i = 2, Nlin - 1
    call locate(x, 1, N, xlin(i), j)
    logx = log(xlin(i))
    logx1 = log(x(j))
    logx2 = log(x(j + 1))
    if (y(j) > 0..and.y(j + 1) > 0.) then
      logy1 = log(y(j))
      logy2 = log(y(j + 1))
      ylin(i) = exp(logy1 + (logx - logx1) / (logx2 - logx1) * (logy2 - logy1))
    else
      ylin(i) = y(1) + (xlin(i) - x(j)) / (x(j + 1) - x(j)) * (y(j + 1) - y(j))
    endif
  enddo
  return
end subroutine linear
! Copyright A.J. Koning 2021
