subroutine write32
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF32
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! Variables for initialization of ENDF format
!   AWR         ! standard mass parameter
!   blank2      ! blank string
!   FEND        ! ENDF - 6 format
!   MAT         ! MAT number
!   mfexist     ! flag for existence of MF - number
!   SEND        ! ENDF - 6 format
!   ZA          ! standard charge parameter
! Variables for MF2
!   ABN         ! abundance
!   AWRI        ! ratio of isotope mass to neutron
!   AP          ! scattering radius
!   EH          ! boundary for resonance range
!   EL          ! boundary for resonance range
!   LFW         ! flag for average fission width
!   Lres        ! l - value
!   LRF         ! representation indicator
!   LRU         ! flag for resolved / unresolved
!   NAPS        ! flag for channel radius and scattering radius
!   NER         ! number of resonance energy ranges
!   NIS         ! number of isotopes in the material
!   NRO         ! flag for energy dependence of scattering radius
!   SPI         ! target spin
!   ZAI         ! (Z, A) designation
! Variables for MF31_40
!   AJ32        ! spin of the resonance
!   APLQX       ! l - dependent scattering radius
!   b32         ! covariance matrix element
!   b32URR      ! covariance matrix element
!   covdigit    ! elements of compact covariance format
!   covix32     ! covariance index
!   D32         ! MF32 resonance parameter
!   DAP         ! uncertainty in scattering radius (compact format only
!   GF32        ! fission width of the resonance
!   GG32        ! gamma width of the resonance
!   GNO32       ! neutron width of the resonance
!   GX32        ! competitive width of the resonance
!   ISR         ! flag for presence of scattering radius uncertainty
!   LCOMP       ! compatibility flag
!   LRX32       ! flag to indicate competitive width
!   MLS         ! number of DAP points
!   MPAR        ! number of parameters per resonance
!   MPARURR     ! number of parameters per resonance for URR
!   N32         ! number of values for MF32
!   N32URR      ! number of points for URR
!   NDIGIT      ! integer for compact covariance format
!   NJS32       ! number of j - values
!   NLRS        ! number of subsections with long - range covariances
!   NLS32       ! number of l - values
!   NM          ! integer for compact covariance format
!   NNN         ! integer for compact covariance format
!   NPARURR     ! number of parameters per resonance for URR
!   NRB         ! number of resonances
!   NSRS        ! number of subsections with covariances
!
! *** Declaration of local data
!
  implicit none
  integer :: i               ! counter
  integer :: j               ! counter
  integer :: l               ! counter
  integer :: MF              ! MF-number
  integer :: MT              ! MT-number
  integer :: n               ! counter
  integer :: NS              ! line number
!
! ***************************** Write MF2 ******************************
!
! hrwrite: subroutine to write header with real values
!
  MF = 32
  if ( .not. mfexist(MF)) return
  MT = 151
  NS = 0
  open (unit = 2, file = 'MF32', status = 'unknown')
  call hrwrite(ZA, AWR, 0, 0, NIS, 0, MAT, MF, MT, NS)
  call hrwrite(ZAI, ABN, 0, LFW, NER, 0, MAT, MF, MT, NS)
  do n = 1, NER
    call hrwrite(EL(n), EH(n), LRU(n), LRF(n), NRO(n), NAPS(n), MAT, MF, MT, NS)
!
! RRR
!
! xwrite   : subroutine to write real value block
!
    if (LRU(n) == 1) then
      call hrwrite(SPI(n), AP(n), 0, LCOMP(n), NLS32(n), ISR(n), MAT, MF, MT, NS)
!
! Uncertainty of scattering radius
!
      if (ISR(n) == 1) then
        if (LRF(n) <= 2) then
          call hrwrite(0., DAP(1), 0, 0, 0, 0, MAT, MF, MT, NS)
        else
          call hrwrite(0., 0., 0, 0, MLS(n), 1, MAT, MF, MT, NS)
          call xwrite(MLS(n), DAP, MAT, MF, MT, NS)
        endif
      endif
!
! Conventional format (non-compact)
!
      if (LCOMP(n) <= 1) then
        call hrwrite(AWRI(n, 1), 0., 0, 0, NSRS, NLRS, MAT, MF, MT, NS)
        if (LRF(n) <= 2) call hrwrite(0., 0., MPAR, 0, N32, NRB, MAT, MF, MT, NS)
        call xwrite(N32, b32, MAT, MF, MT, NS)
      else
!
! Compact format
!
        call hrwrite(AWRI(n, 1), APLQX, 0, LRX32, NSRS, NLRS, MAT, MF, MT, NS)
        call xwrite(N32, b32, MAT, MF, MT, NS)
        call hrwrite(0., 0., NDIGIT, NNN, NM, 0, MAT, MF, MT, NS)
        do i = 1, NM
          write(2, '(2i5, 14a4, i4, i2, i3, i5)') (covix32(i, j), j = 1, 2), (covdigit(i, j), j = 1, 14), MAT, MF, MT, NS
        enddo
      endif
    else
!
! URR
!
      call hrwrite(SPI(n), AP(n), 0, 0, NLS32(n), 0, MAT, MF, MT, NS)
      do l = 1, NLS32(n)
        call hrwrite(AWRI(n, l), 0., Lres(n, l), 0, 6 * NJS32(n, l), NJS32(n, l), MAT, MF, MT, NS)
        do j = 1, NJS32(n, l)
          call rwrite(D32(n, l, j), AJ32(n, l, j), GNO32(n, l, j), GG32(n, l, j), GF32(n, l, j), GX32(n, l, j), MAT, MF, MT, NS)
        enddo
      enddo
      call hrwrite(0., 0., MPARURR, 0, N32URR, NPARURR, MAT, MF, MT, NS)
      call xwrite(N32URR, b32URR, MAT, MF, MT, NS)
    endif
  enddo
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write32
! Copyright A.J. Koning 2021
