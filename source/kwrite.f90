subroutine kwrite(Nx, Nval, MATNUM, MF, MT, NS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write integer value block
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
! Variables for initialization of ENDF format
!   INT3    ! ENDF - 6 format
!
! *** Declaration of local data
!
  implicit none
  integer :: i                     ! counter
  integer :: j                     ! counter
  integer :: MATnum                ! material number
  integer :: MF                    ! MF-number
  integer :: NS                    ! line number
  integer :: MT                    ! MT-number
  integer :: nlin                  ! number of lines
  integer :: nrest                 ! help variable
  integer :: Nx                    ! number of values
  integer :: Nval(Nx)              ! value
!
! ***************************** Write VALUE block **********************
!
! Write values
!
  do i = 1, Nx-5, 6
    write(2, fmt = INT3) (Nval(j), j = 1, 6), MATnum, MF, MT, NS
  enddo
  nlin = Nx / 6
  nrest = Nx - nlin * 6
  if (nrest == 1) write(2, '(i11, t67, i4, i2, i3, i5)') Nval(nlin * 6 + 1), MATnum, MF, MT, NS
  if (nrest == 2) write(2, '(2i11, t67, i4, i2, i3, i5)') (Nval(nlin * 6 + j), j = 1, 2), MATnum, MF, MT, NS
  if (nrest == 3) write(2, '(3i11, t67, i4, i2, i3, i5)') (Nval(nlin * 6 + j), j = 1, 3), MATnum, MF, MT, NS
  if (nrest == 4) write(2, '(4i11, t67, i4, i2, i3, i5)') (Nval(nlin * 6 + j), j = 1, 4), MATnum, MF, MT, NS
  if (nrest == 5) write(2, '(5i11, t67, i4, i2, i3, i5)') (Nval(nlin * 6 + j), j = 1, 5), MATnum, MF, MT, NS
  return
end subroutine kwrite
! Copyright A.J. Koning 2021
