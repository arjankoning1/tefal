subroutine read32
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF32 from existing ENDF-6 data library
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
  use A1_error_handling_mod
!
! Variables for input of specific ENDF data
!   adopt        ! logical for existence of MF information (per MT)
!   adoptfile    ! name of library for MF information (per MT)
! Variables for initialization of ENDF format
!   mfexist      ! flag for existence of MF - number
!   mtexist      ! flag for existence of MT - number
! Variables for MF2
!   AP           ! scattering radius
!   AWRI         ! ratio of isotope mass to neutron
!   EH           ! boundary for resonance range
!   EL           ! boundary for resonance range
!   LFW          ! flag for average fission width
!   Lres         ! l - value
!   LRF          ! representation indicator
!   LRU          ! flag for resolved / unresolved
!   NAPS         ! flag for channel radius and scattering radius
!   NER          ! number of resonance energy ranges
!   NIS          ! number of isotopes in the material
!   NRO          ! flag for energy dependence of scattering radius
!   SPI          ! target spin
! Variables for MF31_40
!   AJ32         ! spin of the resonance
!   APLQX        ! l - dependent scattering radius
!   b32          ! covariance matrix element
!   b32URR       ! covariance matrix element
!   covdigit     ! elements of compact covariance format
!   covix32      ! covariance index
!   D32          ! MF32 resonance parameter
!   DAP          ! uncertainty in scattering radius (compact format only
!   GF32         ! fission width of the resonance
!   GG32         ! gamma width of the resonance
!   GNO32        ! neutron width of the resonance
!   GX32         ! competitive width of the resonance
!   ISR          ! flag for presence of scattering radius uncertainty
!   LCOMP        ! compatibility flag
!   LRX32        ! flag to indicate competitive width
!   MLS          ! number of DAP points
!   MPAR         ! number of parameters per resonance
!   MPARURR      ! number of parameters per resonance for URR
!   N32          ! number of values for MF32
!   N32URR       ! number of points for URR
!   NDIGIT       ! integer for compact covariance format
!   NJS32        ! number of j - values
!   NLRS         ! number of subsections with long - range covariances
!   NLS32        ! number of l - values
!   NM           ! integer for compact covariance format
!   NNN          ! integer for compact covariance format
!   NPARURR      ! number of parameters per resonance for URR
!   NRB          ! number of resonances
!   NSRS         ! number of subsections with covariances
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: ii        ! counter
  integer           :: is        ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: jj        ! counter
  integer           :: ll        ! angular momentum
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: n         ! counter
  integer           :: nlin      ! number of lines
!
! ****************** Read covariance parameters from MF32 **************
!
! No covariance data for potential scattering radius option.
!
  MF = 32
  MT = 151
  if (LRU(1) == 0) return
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
  MTstr = '     '
  write(MTstr(1:2), '(i2)') MF
  write(MTstr(3:5), '(i3)') MT
  do
    read(3, '(a80)', iostat = istat) string
    if (istat /= 0) then
      close (unit = 3)
      return
    endif
    if (string(71:75) == MTstr) exit
  enddo
  do is = 1, NIS
    read(3, '(33x, 2i11)', iostat = istat) LFW, NER
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    do n = 1, NER
      read(3, '(2e11.6, 4i11)', iostat = istat) EL(n), EH(n), LRU(n), LRF(n), NRO(n), NAPS(n)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
!
! RRR
!
      if (LRU(n) == 1) then
        read(3, '(33x, 3i11)', iostat = istat) LCOMP(n), NLS32(n), ISR(n)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
!
! Uncertainty of scattering radius
!
        if (ISR(n) == 1) then
          if (LRF(n) <= 2) then
            read(3, '(11x, e11.6)', iostat = istat) DAP(1)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          else
            read(3, '(44x, i11)', iostat = istat) MLS(n)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
            read(3, '(e11.6)', iostat = istat) DAP(1)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          endif
        endif
!
! Conventional format (non-compact)
!
        if (LCOMP(n) <= 1) then
          read(3, '(44x, 2i11)', iostat = istat) NSRS, NLRS
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          if (LRF(n) <= 2) then
            read(3, '(22x, i11, 11x, 2i11)', iostat = istat) MPAR, N32, NRB
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          endif
          nlin = 1 + (N32 - 1) / 6
          do i = 1, nlin
            ii = 6 * (i - 1)
            read(3, '(6e11.6)', iostat = istat) (b32(j), j = ii+1, ii+6)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          enddo
!
! Compact covariance format
!
        else
          read(3, '(11x, e11.6, 11x, 3i11)', iostat = istat) APLQX, LRX32, NSRS, NLRS
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          N32 = NSRS
          nlin = 1 + (N32 - 1) / 6
          do i = 1, nlin
            ii = 6 * (i - 1)
            read(3, '(6e11.6)', iostat = istat) (b32(j), j = ii+1, ii+6)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          enddo
          read(3, '(22x, 3i11)', iostat = istat) NDIGIT, NNN, NM
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          do i = 1, NM
            read(3, '(2i5, 14a4)', iostat = istat) (covix32(i, j), j = 1, 2), (covdigit(i, j), j = 1, 14)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          enddo
        endif
      else
!
! SLBW (URR)
!
        read(3, '(2e11.6, 22x, i11)', iostat = istat) SPI(n), AP(n), NLS32(n)
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        do ll = 1, NLS32(n)
          read(3, '(e11.6, 11x, i11, 22x, i11)', iostat = istat) AWRI(n, ll), Lres(n, ll), NJS32(n, ll)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          do jj = 1, NJS32(n, ll)
            read(3, '(6e11.6)', iostat = istat) D32(n, ll, jj), AJ32(n, ll, jj), GNO32(n, ll, jj), GG32(n, ll, jj),  &
 &            GF32(n, ll, jj), GX32(n, ll, jj)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          enddo
        enddo
        read(3, '(22x, i11, 22x, i11)', iostat = istat) MPARURR, NPARURR
        if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        N32URR = NPARURR * (NPARURR + 1) / 2
        nlin = 1 + (N32URR - 1) / 6
        do i = 1, nlin
          ii = 6 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (b32URR(j), j = ii+1, ii+6)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        enddo
      endif
    enddo
  enddo
  close (unit = 3)
  mfexist(MF) = .true.
  mtexist(MF, MT) = .true.
  return
end subroutine read32
! Copyright A.J. Koning 2021
