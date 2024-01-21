subroutine xwrite(Nx, x, MATnum, MF, MT, NS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write real value block
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
!   sgl       ! single precision kind
! Variables for initialization of ENDF variables
!   blank1    ! blank string
!   VALUE     ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  character(len=11) :: endf      ! function for a number in ENDF-6 format
  character(len=11) :: xx(6)     ! x value
  integer           :: i         ! counter
  integer           :: ii        ! counter
  integer           :: j         ! counter
  integer           :: MATnum    ! material number
  integer           :: MF        ! MF-number
  integer           :: NS        ! line number
  integer           :: MT        ! MT-number
  integer           :: nlin      ! number of lines
  integer           :: nrest     ! help variable
  integer           :: Nx        ! number of values
  real(sgl)         :: x(Nx)     ! help variable
!
! ***************************** Write VALUE block **********************
!
! endf  : function for a number in ENDF-6 format
!
! Write values
!
  nlin = Nx/6
  do i = 1, nlin
    do j = 1, 6
      ii = 6 * (i - 1) + j
      xx(j) = endf(x(ii))
    enddo
    write(2, fmt = VALUE) (xx(j), j = 1, 6), MATnum, MF, MT, NS
  enddo
  nrest = Nx - nlin * 6
  do j = 1, nrest
    ii = 6 * nlin + j
    xx(j) = endf(x(ii))
  enddo
  if (nrest /= 0) write(2, fmt = VALUE) (xx(j), j = 1, nrest), (blank1, j = 1, 6 - nrest), MATnum, MF, MT, NS
  return
end subroutine xwrite
! Copyright A.J. Koning 2021
