subroutine deallocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
  if (allocated(Ehist)) deallocate(Ehist)
  if (allocated(f0ex)) deallocate(f0ex)
  if (allocated(specexcl)) deallocate(specexcl)
  if (allocated(xsang)) deallocate(xsang)
  if (allocated(xsgamdis)) deallocate(xsgamdis)
  if (allocated(Rmt)) deallocate(Rmt)
  if (allocated(relerr)) deallocate(relerr)
  if (allocated(xserr)) deallocate(xserr)
  if (allocated(Rcov)) deallocate(Rcov)
  if (allocated(Rleg)) deallocate(Rleg)
  if (allocated(Rrp)) deallocate(Rrp)
  if (allocated(xsexcliso)) deallocate(xsexcliso)
  if (allocated(branchiso)) deallocate(branchiso)
  if (allocated(Eout)) deallocate(Eout)
  if (allocated(ncleg)) deallocate(ncleg)
  if (allocated(cleg0)) deallocate(cleg0)
  if (allocated(xsang)) deallocate(xsang)
  if (allocated(fang)) deallocate(fang)
  if (allocated(xsrpiso)) deallocate(xsrpiso)
  if (allocated(Yrpiso)) deallocate(Yrpiso)
!
  call deallocate_recoil
  call deallocate_ddx
  return
end subroutine deallocate_arrays
