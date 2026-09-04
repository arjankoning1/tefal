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
  if (allocated(xsgamdis)) deallocate(xsgamdis)
  if (allocated(b33read)) deallocate(b33read)
  if (allocated(b33)) deallocate(b33)
  if (allocated(b33MT)) deallocate(b33MT)
  if (allocated(b33ZA)) deallocate(b33ZA)
  if (allocated(b33MTread)) deallocate(b33MTread)
  if (allocated(b8)) deallocate(b8)
  if (allocated(b8read)) deallocate(b8read)
  if (allocated(Rmt)) deallocate(Rmt)
  if (allocated(relerr)) deallocate(relerr)
  if (allocated(xserr)) deallocate(xserr)
  if (allocated(Rcov)) deallocate(Rcov)
  if (allocated(Rleg)) deallocate(Rleg)
  if (allocated(Rrp)) deallocate(Rrp)
!
  return
end subroutine deallocate_arrays
