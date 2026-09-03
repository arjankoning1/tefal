subroutine allocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
  allocate(Ehist(0:idnum,1:numenspec,0:numpar,0:numen2))
  Ehist = 0.
  allocate(f0ex(0:idnum,1:numenspec,0:numpar,0:numen2))
  f0ex = 0.
  allocate(specexcl(0:idnum,1:numenspec,0:numpar,0:numen2))
  specexcl = 0.
  return
end subroutine allocate_arrays
