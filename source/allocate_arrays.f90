subroutine allocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
  integer:: Nencovtot
!
  allocate(Ehist(0:idnum,1:numenspec,0:numpar,0:numen2))
  Ehist = 0.
  allocate(f0ex(0:idnum,1:numenspec,0:numpar,0:numen2))
  f0ex = 0.
  allocate(specexcl(0:idnum,1:numenspec,0:numpar,0:numen2))
  specexcl = 0.
  if (flagcovar) then
    Nencovtot = 1 + Nencov * Nencov
    allocate(b33read(Nchancov,Nchancov,Nencovtot))
    allocate(b33(Nchancov,Nchancov,Nencovtot))
    allocate(b33MT(idnum,Nencovtot))
    allocate(b33ZA(idnum,Nencovtot))
    allocate(b8(idnum,Nencovtot))
    allocate(b33MTread(idnum,Nencovtot))
    allocate(b8read(idnum,Nencovtot))
  endif
  return
end subroutine allocate_arrays
