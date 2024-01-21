subroutine read12(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF12 from existing ENDF-6 data library
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
!   adopt          ! logical for existence of MF information (per MT)
!   adoptfile      ! name of library for MF information (per MT)
! Variables for partial cross sections in ENDF format
!   idchannel      ! identifier for channel
! Variables for initialization of ENDF format
!   idnum          ! number of different exclusive cross sections
!   MTid           ! channel identifier for MT - number
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NK             ! number of subsections
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF12_15
!   E12            ! incident energy (in ENDF - 6 format)
!   Eg             ! gamma energy
!   Egk            ! gamma energy
!   Esk            ! starting level (in ENDF - 6 format)
!   INTERg         ! interpolation scheme
!   LFg            ! photo energy distribution law
!   LG12           ! type setters
!   LO12           ! type setters
!   LPg            ! primary photon flag
!   NBTg           ! separation value for interpolation scheme
!   NPg            ! number of incident energies
!   NRg            ! number of interpolation ranges
!   xsgtotyield    ! total discrete photon multiplicity
!   xsgyield       ! gamma - ray multiplicity (in ENDF - 6 format)
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=5)  :: MTstr     ! string for MT number
  character(len=80) :: string    ! line with parameter value
  integer           :: i         ! counter
  integer           :: id        ! counter for deuterons
  integer           :: idc       ! help variable
  integer           :: ii        ! counter
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: k         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: nlin      ! number of lines
!
! ********* Read data from particular MT-numbers from MF12 files *******
!
  MF = 12
  if ( .not. adopt(MF, MT)) return
  open (unit = 3, file = adoptfile(MF, MT), status = 'old')
!
! Read until MT is found
!
! Total gamma yield
!
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
  read(string(23:33), '(i11)') LO12(MT)
  read(string(34:44), '(i11)') LG12(MT)
  read(string(45:55), '(i11)') NK(MF, MT)
  if (LO12(MT) == 1) then
    read(3, '(44x, 2i11)', iostat = istat) NR(MF, MT), NP(MF, MT)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    nlin = 1 + (NR(MF, MT) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6i11)', iostat = istat) (NBT(MF, MT, j), INTER(MF, MT, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
    do id = 0, idnum
      if (idchannel(id) == MTid(MT)) then
        idc = id
        exit
      endif
    enddo
    nlin = 1 + (NP(MF, MT) - 1) / 3
    do i = 1, nlin
      ii = 3 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (E12(idc, j), xsgtotyield(idc, j), j = ii+1, ii+3)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    enddo
!
! Gamma yield per discrete level
!
    if (NK(MF, MT) > 1) then
      do k = 1, NK(MF, MT)
        read(3, '(2e11.6, 4i11)', iostat = istat) Egk(idc, k), Esk(idc, k), LPg(MT, k), LFg(MT, k), NRg(MT, k), NPg(MT, k)
        nlin = 1 + (NRg(MT, k) - 1) / 3
        do i = 1, nlin
          ii = 3 * (i - 1)
          read(3, '(6i11)', iostat = istat) (NBTg(MT, k, j), INTERg(MT, k, j), j = ii + 1, ii + 3)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        enddo
        nlin = 1 + (NPg(MT, k) - 1) / 3
        do i = 1, nlin
          ii = 3 * (i - 1)
          read(3, '(6e11.6)', iostat = istat) (Eg(idc, k, j), xsgyield(idc, k, j), j = ii + 1, ii + 3)
          if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
        enddo
      enddo
    endif
  endif
  close (unit = 3)
  return
end subroutine read12
! Copyright A.J. Koning 2021
