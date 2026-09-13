subroutine allocate_mf12_15
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate arrays for MF12-15
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tefal_mod
  implicit none
!
! MF12
!
  allocate(INTERg(nummt,numgam,numint))
  allocate(LG12(nummt))
  allocate(LO12(nummt))
  allocate(NBTg(nummt,numgam,numint))
  allocate(NPg(nummt,numgam))
  allocate(NRg(nummt,numgam))
  allocate(LP12(nummt))
  allocate(LPg(nummt,numgam))
  allocate(LFg(nummt,numgam))
  allocate(NS12(nummt))
  allocate(NT12(nummt))

  allocate(E12(0:idnum,numenin))
  allocate(Eg(0:idnum,numgam,numenin))
  allocate(Egk(0:idnum,numgam))
  allocate(ES12(nummt,numgam))
  allocate(Esk(0:idnum,numgam))
  allocate(ESNS(nummt))
  allocate(TP12(nummt,numgam))
  allocate(xsgtotyield(0:idnum,numenin))
  allocate(xsgyield(0:idnum,numgam,numenin))
!
! MF14
!
  allocate(LI14(nummt))
!
! MF15
!
  allocate(INTER15(numsecg,numint))
  allocate(INTER15g(numsecg,numint))
  allocate(INTER15ge(numsecg,numenin,numint))
  allocate(NBT15(numsecg,numint))
  allocate(NBT15g(numsecg,numint))
  allocate(NBT15ge(numsecg,numenin,numint))
  allocate(NE15g(numsecg))
  allocate(NP15(numsecg))
  allocate(NP15ge(numsecg,numenin))
  allocate(NR15(numsecg))
  allocate(NR15g(numsecg))
  allocate(NR15ge(numsecg,numenin))

  allocate(E15(numsecg,numenin))
  allocate(E15ge(numsecg,numenin,3*numen2))
  allocate(EPy(numsecg,numenin))
  allocate(ge(numsecg,numenin,3*numen2))
  allocate(Pg(numsecg,numenin))
!
! Initialization
!
  INTERg = 0
  LG12 = 0
  LO12 = 0
  NBTg = 0
  NPg = 0
  NRg = 0
  LP12 = 0
  LPg = 0
  LFg = 0
  NS12 = 0
  NT12 = 0

  E12 = 0.
  Eg = 0.
  Egk = 0.
  ES12 = 0.
  Esk = 0.
  ESNS = 0.
  TP12 = 0.
  xsgtotyield = 0.
  xsgyield = 0.

!
! MF13
!
  if (flaggam13) then
    allocate(E13(0:idnum,numenin))
    allocate(xsg(0:idnum,numgam,numenin))
    allocate(xsgtot(0:idnum,numenin))

    E13 = 0.
    xsg = 0.
    xsgtot = 0.
  endif

  LI14 = 0

  INTER15 = 0
  INTER15g = 0
  INTER15ge = 0
  NBT15 = 0
  NBT15g = 0
  NBT15ge = 0
  NE15g = 0
  NP15 = 0
  NP15ge = 0
  NR15 = 0
  NR15g = 0
  NR15ge = 0

  E15 = 0.
  E15ge = 0.
  EPy = 0.
  ge = 0.
  Pg = 0.
!
  return
end subroutine allocate_mf12_15
