subroutine deallocate_recoil
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate arrays for recoil information
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  if (allocated(nbegrec)) deallocate(nbegrec)
  if (allocated(nendrec)) deallocate(nendrec)
  if (allocated(noutrec)) deallocate(noutrec)

  if (allocated(Erec)) deallocate(Erec)
  if (allocated(recexcl)) deallocate(recexcl)
  if (allocated(Ehistrec)) deallocate(Ehistrec)
  if (allocated(f0exrec)) deallocate(f0exrec)

  if (allocated(nbegcumrec)) deallocate(nbegcumrec)
  if (allocated(nendcumrec)) deallocate(nendcumrec)
  if (allocated(noutrecrp)) deallocate(noutrecrp)

  if (allocated(Erecrp)) deallocate(Erecrp)
  if (allocated(recrp)) deallocate(recrp)
  if (allocated(Ehistcumrec)) deallocate(Ehistcumrec)
  if (allocated(f0cumrec)) deallocate(f0cumrec)

  return
end subroutine deallocate_recoil
