subroutine write33(MF)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF33 and MF40
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
!   sgl            ! single precision kind
! All global variables
!   numencovtot    ! number of energies for covariances
!   nummt          ! number of MT numbers
! Variables for input of ENDF library type
!   flageaf      ! flag for EAF - formatted activation library
! Variables for covariances in ENDF format
!   flagMTint      ! flag for channel with inter - MT covariance data
!   MTindex        ! index for MT number
!   MTindexiso     ! index for isomer of MT number
!   MTintindex     ! index for MT number with inter - MT covariance data
!   Nchancov       ! number of channels with covariance data
! Variables for initialization of ENDF format
!   AWR            ! standard mass parameter
!   blank2         ! blank string
!   FEND           ! ENDF - 6 format
!   LIS            ! state number of target nucleus
!   MAT            ! MAT number
!   mfexist        ! flag for existence of MF - number
!   mtexist        ! flag for existence of MT - number
!   SEND           ! ENDF - 6 format
!   TEXT           ! ENDF - 6 format
!   ZA             ! standard charge parameter
! Variables for MF3
!   eafstring      ! string with reaction information for EAF format
!   QM             ! Q - value (in ENDF - 6 format)
! Variables for MF8_10
!   IZAP           ! second IZAP - number
!   LFS            ! final state number
!   LFSZA          ! final state number
!   NSt            ! number of final states
!   NZA            ! number of nuclides
!   QIiso          ! Q - value for isomer (in ENDF - 6 format)
!   QIZA           ! Q - value (in ENDF - 6 format)
!   QMZA           ! Q - value (in ENDF - 6 format)
!   XMFZA          ! second MF - number
!   ZAPi           ! designation of final nucleus
! Variables for MF31_40
!   b33            ! covariance matrix element
!   b33MT          ! covariance matrix element
!   b33ZA          ! covariance matrix element
!   b8             ! covariance matrix element
!   LB             ! flag for meaning of numbers
!   LB8            ! flag for meaning of numbers
!   LBZA           ! flag for meaning of numbers
!   LS             ! symmetry flag
!   LSZA           ! symmetry flag
!   MAT1           ! second MAT - number
!   MAT1ZA         ! second MAT - number
!   MT33           ! second MT - number
!   MT33ZA         ! second MT - number
!   MTL            ! lumped reaction identifier
!   NC33           ! number of NC - type sub - subsections
!   NC33ZA         ! number of NC - type sub - subsections
!   NE33           ! number of energies in energy array
!   NE33ZA         ! number of energies in energy array
!   NE8            ! number of entries for LB = 8 section
!   NI33           ! number of NI - type sub - subsections
!   NI33ZA         ! number of NI - type sub - subsections
!   NL33           ! number of subsections
!   NT33           ! total number of entries
!   NT33ZA         ! total number of entries
!   NT8            ! total number of entries
!   XLFS1          ! second state discrete level number
!   XMF1           ! second MF - number
!
! *** Declaration of local data
!
  implicit none
  character(len=4) :: cfile              ! covariance file
  integer          :: i                  ! counter
  integer          :: i1                 ! value
  integer          :: ichan              ! counter for channels
  integer          :: iso                ! counter for isomer
  integer          :: iza                ! counter for Z,A combinations
  integer          :: k                  ! counter
  integer          :: MF                 ! MF-number
  integer          :: MT                 ! MT-number
  integer          :: MT2                ! MT number
  integer          :: N                  ! neutron number of residual nucleus
  integer          :: Niso               ! number of isotopes produced after irradiation
  integer          :: NS                 ! line number
  real(sgl)        :: x(numencovtot)     ! help variable
!
! ***************************** Write MF33 *****************************
!
! hrwrite   : subroutine to write header with real values
!
  NS = 0
  cfile = 'MF  '
  write(cfile(3:4), '(i2)') MF
  open (unit = 2, file = cfile, status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(MF, MT)) cycle
    if (flageaf) write(2, fmt = TEXT) eafstring(MT), MAT, MF, MT, NS
    if (MF == 33) then
      call hrwrite(ZA, AWR, 0, MTL(MT), 0, NL33(MF, MT), MAT, MF, MT, NS)
      Niso = 1
    else
      if (MT == 5) then
        call hrwrite(ZA, AWR, LIS, 0, NZA, 0, MAT, MF, MT, NS)
      else
        call hrwrite(ZA, AWR, LIS, 0, NSt(MT), 0, MAT, MF, MT, NS)
        Niso = NSt(MT)
      endif
    endif
!
! 1. Residual production covariances
!
    if (MF == 40 .and. MT == 5) then
      do iza = 1, NZA
        call hrwrite(QMZA(iza), QIZA(iza), IZAP(iza), LFSZA(iza), 0, NL33(MF, MT), MAT, MF, MT, NS)
        call hrwrite(real(XMFZA(iza)), real(LFSZA(iza)), MAT1ZA(iza, 1), MT33ZA(iza, 1), NC33ZA(iza, 1), NI33ZA(iza, 1), &
 &        MAT, MF, MT, NS)
        call hrwrite(0., 0., LSZA(iza, 1), LBZA(iza, 1), NT33ZA(iza, 1), NE33ZA(iza, 1), MAT, MF, MT, NS)
        N = NT33ZA(iza, 1)
        do i = 1, N
          x(i) = b33ZA(iza, i)
        enddo
        call xwrite(N, x, MAT, MF, MT, NS)
      enddo
    else
!
! 2. Cross section and isomeric covariances
!
      do iso = 1, Niso
        do ichan = 1, Nchancov
          if (MT /= MTindex(ichan)) cycle
          if (MF == 33 .and. MTindexiso(ichan) /=  -1) cycle
          if (MF == 40 .and. MT /= 18 .and. iso /= MTindexiso(ichan) + 1) cycle
!
! 1. Covariance data per MT number
!
! xwrite    : subroutine to write real value block
!
! Intra-MT correlations
!
          if (MF == 40) call hrwrite(QM(MT), QIiso(MT, iso), ZAPi(MT), LFS(MT, iso), 0, NL33(MF, MT), MAT, MF, MT, NS)
          call hrwrite(real(XMF1(ichan, 1)), real(LFS(MT, iso)), &
 &          MAT1(ichan, 1), MT33(ichan, 1), NC33(ichan, 1), NI33(ichan, 1), MAT, MF, MT, NS)
          if ( .not. flageaf) then
            call hrwrite(0., 0., LS(ichan, 1), LB(ichan, 1), NT33(ichan, 1), NE33(ichan, 1), MAT, MF, MT, NS)
            N = NT33(ichan, 1)
            do i = 1, N
              x(i) = b33MT(ichan, i)
            enddo
            call xwrite(N, x, MAT, MF, MT, NS)
          endif
!
! Sub-subsection LB=8 for variances
!
          if (MF == 33) then
            call hrwrite(0., 0., 0, LB8, NT8(ichan), NE8(ichan), MAT, MF, MT, NS)
            N = NT8(ichan)
            do i = 1, N
              x(i) = b8(ichan, i)
            enddo
            call xwrite(N, x, MAT, MF, MT, NS)
          endif
          if (flageaf) cycle
!
! Inter-MT correlations
!
          if (MF == 33 .and. flagMTint(ichan)) then
            i1 = MTintindex(ichan)
            do MT2 = MT, nummt
              do k = 2, NL33(MF, MT)
                if (MT2 /= MT33(ichan, k)) cycle
                call hrwrite(real(XMF1(ichan, k)), real(XLFS1(ichan, k)), MAT1(ichan, k), MT33(ichan, k), &
 &                NC33(ichan, k), NI33(ichan, k), MAT, MF, MT, NS)
                call hrwrite(0., 0., LS(ichan, k), LB(ichan, k), NT33(ichan, k), NE33(ichan, k), MAT, MF, MT, NS)
!
! Write covariance data
!
                N = NT33(ichan, k)
                do i = 1, N
                  x(i) = b33(i1, k, i)
                enddo
                call xwrite(N, x, MAT, MF, MT, NS)
              enddo
            enddo
          endif
        enddo
      enddo
    endif
    write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  enddo
  if (mfexist(MF)) write(2, fmt = FEND) blank2, MAT, 0, 0, NS
  close (unit = 2)
  return
end subroutine write33
! Copyright A.J. Koning 2021
