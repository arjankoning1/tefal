subroutine allocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
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
  if (flagendfdet .or. flageaf) then
    allocate(xsexcliso(0:idnum,0:nlevmax,1:numinc))
    allocate(branchiso(0:idnum,0:nlevmax,1:numinc))

    xsexcliso = 0.
    branchiso = 0.
  endif
  if (flaggpf) then
    allocate(Eout(0:idnum,1:Nenspec,0:numen2))
    Eout = 0.
  endif
  ! Angular arrays are also used by processangle and MF4 when endfdetail is disabled.
  if (flaggpf) then
    allocate(ncleg(0:numpar,0:numlevin,Nenang))
    allocate(cleg0(0:numpar,0:numlevin,Nenang,0:numl))
    allocate(xsang(0:numpar,0:numlevin,0:Nenang,0:numang))
    allocate(fang(0:numpar,0:numlevin,0:Nenang,0:numang))

    ncleg = 0
    cleg0 = 0.
    xsang = 0.
    fang = 0.
  endif
  if (flagcovar) then
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
  call allocate_recoil
  return
end subroutine allocate_arrays
