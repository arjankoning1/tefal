subroutine talystotal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read total cross sections from TALYS
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
! Variables for TALYS info
!   eninc           ! incident energy
!   flagtalysdet    ! flag for detailed ENDF - 6 information from TALYS
!   k0              ! index of incident particle
!   numinc          ! number of incident energies
! Constants
!   pi          ! pi
! Variables for ENDF limits, switches and tolerances
!   Eswitch     ! energy where ENDF - 6 representation is switched (in MeV)
! Variables for total cross sections in ENDF format
!   eninc6      ! incident energy for total cross section
!   numcut6     ! number of energies before high - energy format for total cross sections
!   numinc6     ! number of incident energies for total cross sections
!   Rprime      ! potential scattering radius
!   xselas      ! total elastic cross section
!   xselas6     ! total elastic cross section
!   xsnonel     ! nonelastic cross section
!   xsnonel6    ! nonelastic cross section
!   xstot6      ! total cross section
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: line
  character(len=132) :: key
  character(len=12) :: xsfile         ! file with cross sections
  integer :: istat            ! logical for file access
  integer :: N                ! neutron number of residual nucleus
  integer :: keyix          !
  integer :: nin              ! counter for incident energy
!
! ********************** Read total cross sections *********************
!
! A. Total cross sections for quantities in other files.
!
  xsfile='all.tot'
  open (unit = 1, file = xsfile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(xsfile, istat)
  do
    read(1,'(a)', iostat = istat) line
    if (istat == -1) exit
    key='entries'
    keyix=index(line,trim(key))
    if (keyix > 0) then
      read(line(keyix+len_trim(key)+2:80),*, iostat = istat) N
      if (istat /= 0) call read_error(xsfile, istat)
      read(1,'(/)')
      if (k0 == 1) then
        do nin = 1, N
          read(1, '(15x, 2es15.6)', iostat = istat) xsnonel(nin), xselas(nin)
          if (istat == -1) exit
          if (istat > 0) call read_error(xsfile, istat)
        enddo
      else
        do nin = 1, N
          read(1, '(15x, 2es15.6)', iostat = istat) xsnonel(nin)
          if (istat == -1) exit
          if (istat > 0) call read_error(xsfile, istat)
        enddo
      endif
      exit
    endif
  enddo
  close (unit = 1)
!
! A. Total cross sections for MF3/MT1-3.
!    For storage of the total and total non-elastic cross sections on an ENDF-6 file, we want a finer grid than that used for
!    the TALYS calculations which produces all the partial cross sections.
!    Variable names with a 6 (for ENDF-6 format) in their name are used for this.
!
  if (k0 == 0) then
    numinc6 = numinc
    do nin = 1, numinc
      eninc6(nin) = eninc(nin)
      xsnonel6(nin) = xsnonel(nin)
    enddo
  else
    numinc6 = 0
    xsfile='endf.tot'
    open (unit = 1, file = xsfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(xsfile, istat)
    do
      read(1,'(a)', iostat = istat) line
      if (istat == -1) exit
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) numinc6
        if (istat /= 0) call read_error(xsfile, istat)
        read(1,'(/)')
        if (k0 == 1) then
          do nin = 1, numinc6
            read(1, '(4es15.6)', iostat = istat) eninc6(nin), xsnonel6(nin), xselas6(nin), xstot6(nin)
            if (istat /= 0) call read_error(xsfile, istat)
          enddo
        else
          do nin = 1, numinc6
            read(1, '(2es15.6)', iostat = istat) eninc6(nin), xsnonel6(nin)
            if (istat /= 0) call read_error(xsfile, istat)
          enddo
        endif
        exit
      endif
    enddo
    close (unit = 1)
  endif
  numcut6 = numinc6
  do nin = 1, numinc6
    if (eninc6(nin) > Eswitch) then
      numcut6 = nin - 1
      exit
    endif
  enddo
  numcut6 = max(numcut6, 1)
!
! ************ For neutrons: read potential scattering radius **********
!
! We use Rprime in case we construct the most simple resonance file.
!
  if (k0 == 1 .and. flagtalysdet) then
    Rprime = 10. * sqrt(0.001 * max(real(xselas(1)), 0.) / (4.*pi))
    open (unit = 1, file = 'spr.opt', status = 'old', iostat=istat)
    if (istat == 0) then
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='Rprime [fm]'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Rprime
          if (istat /= 0) call read_error('spr.opt', istat)
          exit
        endif
      enddo
      close (unit = 1)
    endif
  endif
  return
end subroutine talystotal
! Copyright A.J. Koning 2021
