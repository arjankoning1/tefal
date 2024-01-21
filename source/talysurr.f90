subroutine talysurr
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read URR data from TALYS
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
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Variables for URR in ENDF format
!   Djlurr     ! average level spacing for resonances for j, l value
!   Eurr       ! incident energy with URR data
!   GFjlurr    ! average fission width for j, l value
!   GGjlurr    ! average radiative width for j, l value
!   GNjlurr    ! average reduced neutron width for j, l value
!   GXjlurr    ! average fission width for j, l value
!   Jurr       ! spin value
!   NJSurr     ! number of j - values
!   NLSurr     ! number of l - values
!   Nurr       ! number of incident energies with URR data
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical          :: lexist      ! logical to determine existence
  character(len=7) :: urrfile     ! name of URR file
  character(len=132) :: line
  character(len=132) :: key
  integer          :: istat       ! logical for file access
  integer          :: jix         ! index for j
  integer          :: k           ! counter
  integer          :: l           ! counter
  integer          :: N           ! counter
  integer          :: keyix
  integer          :: lix         ! index for l
  integer          :: nin         ! counter for incident energy
  real(sgl)        :: Djl         ! average level spacing for resonances for j,l value
  real(sgl)        :: Dl          ! mean s-wave resonance spacing per l value
  real(sgl)        :: GFjl        ! average fission width for j,l value
  real(sgl)        :: GGjl        ! average radiative width for j,l value
  real(sgl)        :: GNjl        ! average reduced neutron width for j,l value
  real(sgl)        :: GXjl        ! average competitive width for j,l value
  real(sgl)        :: rj          ! help variable
  real(sgl)        :: p           ! parity
  real(sgl)        :: Sjl         ! strenght function per j,l
  real(sgl)        :: Sl          ! strenght function per l
!
! **************************** Read URR data **************************
!
  NLSurr = 0
  urrfile = 'urr.dat'
  inquire (file = urrfile, exist = lexist)
  if (lexist) then
    open (unit = 1, file = urrfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(urrfile, istat)
    nin = 1
    do
      read(1,'(a)',iostat = istat) line
      if (istat == -1) exit
      key='E-incident [MeV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Eurr(nin)
      if (istat /= 0) call read_error(urrfile, istat)
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) N
        lix = -1
        jix = 0
        read(1,'(/)')
        do k=1, N
          read(1, * , iostat = istat) l, rj, p, Dl, Djl, Sl, Sjl, GXjl, GNjl, GGjl, GFjl
          if (istat /= 0) exit
          if (l + 1 /= lix) then
            jix = 1
            lix = l + 1
            NLSurr = max(lix, NLSurr)
          else
            jix = jix + 1
          endif
          NJSurr(lix) = jix
          Jurr(lix, jix) = rj
          Djlurr(nin, lix, jix) = Djl
          GXjlurr(nin, lix, jix) = GXjl
          GNjlurr(nin, lix, jix) = GNjl
          GGjlurr(nin, lix, jix) = GGjl
          GFjlurr(nin, lix, jix) = GFjl
        enddo
        nin = nin + 1
      endif
    enddo
    Nurr = nin - 1
    close (unit = 1)
  endif
  return
end subroutine talysurr
! Copyright A.J. Koning 2021
