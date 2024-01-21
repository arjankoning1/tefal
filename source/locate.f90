subroutine locate(xx, ib, ie, x, j)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Find value in ordered table
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
  logical   :: ascend              ! logical for ascending values
  integer   :: ib                  ! counter
  integer   :: ie                  ! counter
  integer   :: j                   ! counter
  integer   :: jl                  ! lower value
  integer   :: jm                  ! middle value
  integer   :: ju                  ! higher value
  real(sgl) :: x                   ! help variable
  real(sgl) :: xx(0:ie)            ! x value
!
! ******************************* Search *******************************
!
! Find j such that xx(j) <= x < xx(j+1) or xx(j) > x >= xx(j+1)
!
  jl = ib-1
  ju = ie + 1
  ascend = xx(ie) >= xx(ib)
  10  if(ju - jl > 1)then
    jm = (ju + jl) / 2
    if(ascend.eqv.(x >= xx(jm)))then
      jl = jm
    else
      ju = jm
    endif
    goto 10
  endif
  if(x == xx(ib))then
    j = ib
  else if(x == xx(ie))then
    j = ie - 1
  else
    j = jl
  endif
  return
end subroutine locate
! Copyright A.J. Koning 2021
