subroutine allocate_mf31
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF31
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  integer :: nb31
!
! read31 reads data in blocks of 6 values. Allow for the final
! partially filled block.
!
  nb31 = 6 * ((10 * numencovtot + 5) / 6)

  allocate(b31(nb31))

  b31 = 0.

  return
end subroutine allocate_mf31
