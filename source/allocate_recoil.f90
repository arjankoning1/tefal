subroutine allocate_recoil
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for recoil information
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
  if (.not. flagrecoil) return
!
! Exclusive-channel recoil information
!
  allocate(nbegrec(0:idnum,1:Nenspec))
  allocate(nendrec(0:idnum,1:Nenspec))
  allocate(noutrec(0:idnum,1:Nenspec))

  allocate(Erec(0:idnum,1:Nenspec,0:numenrec))
  allocate(recexcl(0:idnum,1:Nenspec,0:numenrec))
  allocate(Ehistrec(0:idnum,1:Nenspec,0:numenrec))
  allocate(f0exrec(0:idnum,1:Nenspec,0:numenrec))
!
! Residual-production recoil information
!
  allocate(nbegcumrec(0:numZ,0:numN,1:Nenspec))
  allocate(nendcumrec(0:numZ,0:numN,1:Nenspec))
  allocate(noutrecrp(0:numZ,0:numN,1:Nenspec))

  allocate(Erecrp(0:numZ,0:numN,1:Nenspec,0:numenrec))
  allocate(recrp(0:numZ,0:numN,1:Nenspec,0:numenrec))
  allocate(Ehistcumrec(0:numZ,0:numN,1:Nenspec,0:numenrec))
  allocate(f0cumrec(0:numZ,0:numN,1:Nenspec,0:numenrec))
!
! Initialization
!
  nbegrec = 0
  nendrec = 0
  noutrec = 0

  nbegcumrec = 0
  nendcumrec = 0
  noutrecrp = 0

  Erec = 0.
  recexcl = 0.
  Ehistrec = 0.
  f0exrec = 0.

  Erecrp = 0.
  recrp = 0.
  Ehistcumrec = 0.
  f0cumrec = 0.

  return
end subroutine allocate_recoil
