subroutine talyscovar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read covariance data from TASMAN
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
! Definition of single and double precision variables
!   sgl              ! single precision kind
! All global variables
!   numchan          ! maximum number of exclusive channels
!   numchancov       ! number of channels with inter - channel correlations
! Variables for TALYS info
!   eninc            ! incident energy
! Variables for ENDF covariance input
!   flagcovrp        ! flag for covariance of residual production c.s.
!   flagcross        ! flag for covariance cross - channel correlation
!   flagintercor     ! flag for inter - MT covariance data
! Variables for input of ENDF library type
!   flageaf          ! flag for EAF - formatted activation library
! Variables for input of ENDF structure
!   flagsubfis       ! flag to include subactinide fission
! Variables for covariances in ENDF format
!   Arpcov           ! mass number of residual product
!   Ecov             ! energy grid for covariances
!   Ecovindex        ! energy index for main energy grid
!   Eleg             ! energy for Legendre covariance data
!   flagcovleg       ! flag for covariances for Legendre coefficients
!   flagMTint        ! flag for channel with inter - MT covariance data
!   MTindex          ! index for MT number
!   MTindexiso       ! index for isomer of MT number
!   MTintindex       ! index for MT number with inter - MT covariance data
!   MTintindexiso    ! index for MT number with inter - MT covariance data
!   isorpcov         ! isomer number of residual product
!   Nisocov          ! number of isomers with covariances per MT number
!   Nchancov         ! number of channels with covariance data
!   Nchancovint      ! number of channels with inter - MT covariance data
!   Nchanleg         ! number of energies with Legendre covariance dat
!   Nencov           ! number of energies with covariance data
!   Ncovrp           ! neutron number of residual product
!   Nleg34           ! number of Legendre coefficients with covariances
!   Rcov             ! relative covariance matrix
!   relerr           ! relative cross section uncertainty
!   Rleg             ! covariance matrix element for Legendre coeff.
!   Rmt              ! relative covariance matrix within same MT number
!   Rrp              ! covariance element for residual cross section
!   xserr            ! cross section uncertainty
!   Zrpcov           ! charge number of residual product
!
! *** Declaration of local data
!
  implicit none
  logical   :: lexist                       ! logical to determine existence
  character(len=132) :: line        !
  character(len=132) :: key        !
  integer   :: ichan                        ! counter for channels
  integer   :: ichan2                       ! counter for channels
  integer   :: istat                        ! logical for file access
  integer   :: j                            ! counter
  integer   :: keyix
  integer   :: l                            ! counter
  integer   :: mt                           ! MT number
  integer   :: MTint(numchan)               ! MT-numbers with inter-MT covariance data
  integer   :: MTintiso(numchan)            ! isomer of MT-number with inter-MT covariance data
  integer   :: nen                          ! energy counter
  integer   :: nen2                         ! energy counter
  integer   :: nin                          ! counter for incident energy
  real(sgl) :: ee                           ! energy
  real(sgl) :: ee1                          ! help variable
!
! ****************** Read covariance info ******************************
!
  inquire (file = 'covariance.inf', exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TEFAL-error: Non-existent covariance file")')
    stop
  endif
  open (unit = 1, file = 'covariance.inf', status = 'old')
  read(1, * )
  read(1, * ) Nencov
  do j = 1, Nencov
    read(1, * ) Ecovindex(j)
  enddo
  read(1, * ) Nchancov
  do ichan = 1, Nchancov
    read(1, * ) MTindex(ichan), MTindexiso(ichan)
    if (MTindexiso(ichan) >= 0) Nisocov(MTindex(ichan)) = Nisocov(MTindex(ichan)) + 1
  enddo
  read(1, * ) Nchancovint
  do ichan = 1, Nchancovint
    read(1, * ) MTint(ichan), MTintiso(ichan)
  enddo
  if (Nchancovint > numchancov) then
    write(*, '(" TEFAL-error: Increase numchancov dimension")')
    stop
  endif
  read(1, * ) Nchanleg, Nleg34
  do ichan = 1, Nchanleg
    read(1, * ) Eleg(ichan)
  enddo
  read(1, * ) Ncovrp
  do ichan = 1, Ncovrp
    read(1, * ) Zrpcov(ichan), Arpcov(ichan), isorpcov(ichan)
  enddo
  close (unit = 1)
  if (Nchanleg > 0) then
    flagcovleg = .true.
  else
    flagcovleg = .false.
  endif
  do ichan = 1, Nchancov
    flagMTint(ichan) = .false.
    do ichan2 = 1, Nchancovint
      if (MTint(ichan2) == MTindex(ichan) .and. MTintiso(ichan2) == MTindexiso(ichan) .and. flagcross(MTindex(ichan))) then
        MTintindex(ichan) = ichan2
        MTintindexiso(ichan) = MTintiso(ichan2)
        if (flagintercor) flagMTint(ichan) = .true.
      endif
    enddo
  enddo
  do ichan = 1, Nchancovint
    if (flagsubfis .and. MTint(ichan) == 18) flagMTint(ichan) = .false.
  enddo
!
! ************ Read covariance data within same MT number **************
!
  open (unit = 2, file = 'cov_intra.ave', status = 'old')
Loop1:  do ichan = 1, Nchancov
    do
      read(2,'(a)',iostat = istat) line
      if (istat == -1) exit Loop1
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(2,'(/)')
        do nen = 1, Nencov
          do nen2 = 1, Nencov
!           read(2, * , iostat = istat) mt, ie, ie, ee, mt, ie, ie, ee, Rmt(ichan, nen, nen2)
            read(2, * , iostat = istat) ee, ee, Rmt(ichan, nen, nen2)
            if (istat > 0) cycle
          enddo
        enddo
        exit
      endif
    enddo
  enddo Loop1
  close (unit = 2)
!
! ************* Read covariance data between MT numbers ****************
!
  open (unit = 3, file = 'cov_inter.ave', status = 'old')
Loop2:  do ichan = 1, Nchancovint
    do ichan2 = ichan, Nchancovint
      do
        read(3,'(a)',iostat = istat) line
        if (istat == -1) exit Loop2
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(3,'(/)')
          do nen = 1, Nencov
            do nen2 = 1, Nencov
              read(3, * , iostat = istat) ee, ee, Rcov(ichan, nen, ichan2, nen2)
              if (istat > 0) cycle
            enddo
          enddo
          exit
        endif
      enddo
    enddo
  enddo Loop2
  close (unit = 3)
!
! ************************* Read variance data *************************
!
  if (flageaf) then
    do ichan = 1, Nchancov
      if (MTindexiso(ichan) == 1) MTindex(ichan) = MTindex(ichan) + 300
      if (MTindexiso(ichan) == 2) MTindex(ichan) = MTindex(ichan) + 600
    enddo
  endif
  open (unit = 4, file = 'variance.ave', status = 'old')
Loop3:  do ichan = 1, Nchancov
    do
      read(4,'(a)',iostat = istat) line
      if (istat == -1) exit Loop3
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(4,'(/)')
        do nen = 1, Nencov
          read(4, * , iostat = istat) ee1, ee, relerr(ichan, nen), ee, xserr(ichan, nen)
          if (istat > 0) cycle
        enddo
        exit
      endif
    enddo
  enddo Loop3
  close (unit = 4)
!
! ******************* Read Legendre covariance data ********************
!
  if (flagcovleg) then
    open (unit = 3, file = 'cov_legendre.ave', status = 'old')
Loop4: do ichan = 1, Nchanleg
      do ichan2 = 1, Nchanleg
        do
          read(3,'(a)',iostat = istat) line
          if (istat == -1) exit Loop4
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(3,'(/)')
            do j = 0, Nleg34
              do l = 0, Nleg34
                read(3, * , iostat = istat) mt, mt, Rleg(ichan, j, ichan2, l)
                if (istat > 0) cycle
              enddo
            enddo
            exit
          endif
        enddo
      enddo
    enddo Loop4
    close (unit = 3)
!
! Assign covariance data of first energy to 1.e-5 eV point
!
    Eleg(0) = EminMeV
    do j = 0, Nleg34
      do ichan2 = 0, Nchanleg
        do l = 0, Nleg34
          Rleg(0, j, ichan2, l) = Rleg(1, j, ichan2, l)
        enddo
      enddo
    enddo
  endif
!
! ******* Read residual production cross section covariance data *******
!
  if (flagcovrp) then
    inquire (file = 'cov_residual.ave', exist = lexist)
    if (lexist) then
      open (unit = 3, file = 'cov_residual.ave', status = 'old')
Loop5:do ichan = 1, Ncovrp
        do
          read(3,'(a)',iostat = istat) line
          if (istat == -1) exit Loop5
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(3,'(/)')
            do nen = 1, Nencov
              do nen2 = 1, Nencov
                read(3, * , iostat = istat)  ee, ee, Rrp(ichan, nen, nen2)
                if (istat > 0) cycle
              enddo
            enddo
            exit
          endif
        enddo
      enddo Loop5
      close (unit = 3)
    else
      write(*, '(" TEFAL-warning: No covariance matrix for ", "residual production cross sections available, continuing without")')
      flagcovrp = .false.
    endif
  endif
!
! ****************** Determine covariance energy grid ******************
!
  do nen = 1, Nencov
    nin = Ecovindex(nen)
    Ecov(nen) = eninc(nin)
  enddo
  return
end subroutine talyscovar
! Copyright A.J. Koning 2021
