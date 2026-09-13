subroutine allocate_mf2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate large arrays for MF2
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none

  allocate(AJ(numres,numlres,numnrs))
  allocate(Er(numres,numlres,numnrs))
  allocate(ER7(numjres,numnrs))

  allocate(GF(numres,numlres,numnrs))
  allocate(GFA(numres,numlres,numnrs))
  allocate(GFB(numres,numlres,numnrs))
  allocate(GG(numres,numlres,numnrs))
  allocate(GN(numres,numlres,numnrs))
  allocate(GT(numres,numlres,numnrs))

  allocate(D(numres,numlres,numjres,numnrs))
  allocate(Es(numres,numlres,numjres,numnrs))
  allocate(GFu(numres,numlres,numjres,numnrs))
  allocate(GGu(numres,numlres,numjres,numnrs))
  allocate(GN0(numres,numlres,numjres,numnrs))
  allocate(GX(numres,numlres,numjres,numnrs))

  allocate(GAM7(numjres,numnrs,numch7))

  AJ = 0.
  Er = 0.
  ER7 = 0.

  GF = 0.
  GFA = 0.
  GFB = 0.
  GG = 0.
  GN = 0.
  GT = 0.

  D = 0.
  Es = 0.
  GFu = 0.
  GGu = 0.
  GN0 = 0.
  GX = 0.

  GAM7 = 0.

  return
end subroutine allocate_mf2
