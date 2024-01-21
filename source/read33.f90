subroutine read33(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF33 from existing ENDF-6 data library
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
!   dbl            ! double precision kind
! All global variables
!   numencovtot    ! number of energies for covariances
!   nummt          ! number of MT numbers
! Variables for input of specific ENDF data
!   adoptfile      ! name of library for MF information (per MT)
! Variables for covariances in ENDF format
!   flagMTint      ! flag for channel with inter - MT covariance data
!   MTindex        ! index for MT number
!   MTindexiso     ! index for isomer of MT number
!   MTintindex     ! index for MT number with inter - MT covariance data
!   Nchancov       ! number of channels with covariance data
! Variables for MF31_40
!   b33MTread      ! covariance matrix element
!   b33read        ! covariance matrix element
!   b8read         ! covariance matrix element
!   LBread         ! flag for meaning of numbers
!   LB8read        ! flag for meaning of numbers
!   LSread         ! symmetry flag
!   MAT1read       ! second MAT - number
!   MT33read       ! second MT - number
!   MTLread        ! lumped reaction identifier
!   NC33read       ! number of NC - type sub - subsections
!   NE33read       ! number of energies in energy array
!   NE8read        ! number of entries for LB = 8 section
!   NI33read       ! number of NI - type sub - subsections
!   NL33read       ! number of subsections
!   NT33read       ! total number of entries
!   NT8read        ! total number of entries
!   XLFS1read      ! second state discrete level number
!   XMF1read       ! second MF - number
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: i1        ! value
  integer           :: ichan     ! counter for channels
  integer           :: istat     ! error code
  integer           :: ii        ! counter
  integer           :: j         ! counter
  integer           :: k         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: MT2       ! MT number
  integer           :: nlin      ! number of lines
  real(dbl)         :: bdb(6)    ! covariance matrix element
!
! ****************** Read covariance parameters from MF33 **************
!
  MF = 33
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
  read(string, '(33x, i11, 11x, i11)', iostat = istat) MTLread(MT), NL33read(MF, MT)
  if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
  do ichan = 1, Nchancov
    if (MT /= MTindex(ichan)) cycle
    if (MTindexiso(ichan) /=  -1) cycle
!
! Intra-MT correlations
!
    read(3, '(2e11.6, 4i11)', iostat = istat) XMF1read(ichan, 1), XLFS1read(ichan, 1), &
 &    MAT1read(ichan, 1), MT33read(ichan, 1), NC33read(ichan, 1), NI33read(ichan, 1)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    read(3, '(22x, 4i11)', iostat = istat) LSread(ichan, 1), LBread(ichan, 1), NT33read(ichan, 1), NE33read(ichan, 1)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    if (NE33read(ichan, 1) > numencovtot) then
      write(*, '(" TEFAL-error: MF33 Increase numencovtot")')
      stop
    endif
    nlin = 1 + (NT33read(ichan, 1) - 1) / 6
    do i = 1, nlin
      ii = 6 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (bdb(j), j = 1, 6)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      do j = ii + 1, ii + 6
        if (abs(bdb(j - ii)) >= 1.e-38) b33MTread(ichan, j) = real(bdb(j - ii))
      enddo
    enddo
!
! Sub-subsection LB=8 for variances
!
    read(3, '(33x, 3i11)', iostat = istat) LB8read, NT8read(ichan), NE8read(ichan)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    nlin = 1 + (NT8read(ichan) - 1) / 6
    do i = 1, nlin
      ii = 6 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (bdb(j), j = 1, 6)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      do j = ii + 1, ii + 6
        if (abs(bdb(j - ii)) >= 1.e-38) b8read(ichan, j) = real(bdb(j - ii))
      enddo
    enddo
!
! Inter-MT correlations
!
    if (flagMTint(ichan)) then
      i1 = MTintindex(ichan)
      do MT2 = MT, nummt
        do k = 2, NL33read(MF, MT)
          if (MT2 /= MT33read(ichan, k)) cycle
          read(3, '(2e11.6, 4i11)', iostat = istat) XMF1read(ichan, k), XLFS1read(ichan, k), MAT1read(ichan, k), &
 &          MT33read(ichan, k), NC33read(ichan, k), NI33read(ichan, k)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          read(3, '(22x, 4i11)', iostat = istat) LSread(ichan, k), LBread(ichan, k), NT33read(ichan, k), NE33read(ichan, k)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
          nlin = 1 + (NT33read(ichan, k) - 1) / 6
          do i = 1, nlin
            ii = 6 * (i - 1)
            read(3, '(6e11.6)', iostat = istat) (bdb(j), j = 1, 6)
            if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
            do j = ii + 1, ii + 6
              if (abs(bdb(j - ii)) >= 1.e-38) b33read(i1, k, j) = real(bdb(j - ii))
            enddo
          enddo
        enddo
      enddo
    endif
  enddo
  close (unit = 3)
  return
end subroutine read33
! Copyright A.J. Koning 2021
