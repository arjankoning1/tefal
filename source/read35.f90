subroutine read35(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read MF35 from existing ENDF-6 data library
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
!   dbl          ! double precision kind
! Variables for input of specific ENDF data
!   adoptfile    ! name of library for MF information (per MT)
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF31_40
!   b35          ! covariance matrix element
!   E35b         ! start energy of block
!   E35e         ! end energy of block
!   LB35         ! flag for meaning of numbers
!   LS35         ! symmetry flag
!   NE35         ! number of energies in energy array
!   NT35         ! total number of entries
!   numencov35   ! number of incident energies for FNS covariances
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
  integer           :: istat     ! error code
  integer           :: j         ! counter
  integer           :: k         ! counter
  integer           :: MF        ! MF-number
  integer           :: MT        ! MT-number
  integer           :: nlin      ! number of lines
  real(dbl)         :: bdb(6)    ! covariance matrix element
!
! ****************** Read covariance parameters from MF35 **************
!
  MF = 35
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
  read(string, '(44x, i11)', iostat = istat) NK(MF, MT)
  do k = 1, NK(MF, MT)
    read(3, '(2e11.6, 4i11)', iostat = istat) E35b(k), E35e(k), LS35, LB35, NT35(k), NE35(k)
    if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
    if (NT35(k) > numencov35) then
      write(*, '(" TEFAL-error: MF35 Increase numencov35", " to at least ", i5)') NT35(k)
      stop
    endif
    nlin = 1 + (NT35(k) - 1) / 6
    do i = 1, nlin
      ii = 6 * (i - 1)
      read(3, '(6e11.6)', iostat = istat) (bdb(j), j = 1, 6)
      if (istat /= 0) call read_error(adoptfile(MF, MT), istat)
      do j = ii + 1, ii + 6
        if (abs(bdb(j - ii)) >= 1.e-38) then
          b35(k, j) = real(bdb(j - ii))
        else
          b35(k, j) = 0.
        endif
      enddo
    enddo
  enddo
  close (unit = 3)
  return
end subroutine read35
! Copyright A.J. Koning 2021
