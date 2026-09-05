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
  if (flaggpf) then
    allocate(Ehist(0:idnum,1:Nenspec,0:numpar,0:numen2))
    Ehist = 0.
    allocate(f0ex(0:idnum,1:Nenspec,0:numpar,0:numen2))
    f0ex = 0.
    allocate(specexcl(0:idnum,1:Nenspec,0:numpar,0:numen2))
    specexcl = 0.
  endif
  if (flaggpf .and. flagendfdet) then
    allocate(xsgamdis(0:idnum,1:Nengam,0:numlevels,0:numlevels))
    xsgamdis = 0.
  endif
  if (flagcovar) then
    Nencovtot = 1 + Nencov * Nencov
    allocate(b33read(Nchancov,Nchancov,Nencovtot))
    allocate(b33(Nchancov,Nchancov,Nencovtot))
    allocate(b33MT(Nchancov,Nencovtot))
    allocate(b33ZA(Ncovrp,Nencovtot))
    allocate(b8(Nchancov,Nencovtot))
    allocate(b33MTread(Nchancov,Nencovtot))
    allocate(b8read(Nchancov,Nencovtot))
    allocate(Rmt(Nchancov,Nencov,Nencov))
    Rmt = 0.
    allocate(relerr(Nchancov,Nencov))
    relerr = 0.
    allocate(xserr(Nchancov,Nencov))
    xserr = 0.
    allocate(Rcov(Nchancovint,Nencov,Nchancovint,Nencov))
    Rcov = 0.
    if (flagcovleg) then
      allocate(Rleg(0:Nchanleg,0:Nleg34,0:Nchanleg,0:Nleg34))
      Rleg = 0.
    endif
    if (flagcovrp .and. Ncovrp > 0) then
      allocate(Rrp(Ncovrp,Nencov,Nencov))
      Rrp = 0.
    endif
  endif
  return
end subroutine allocate_arrays
