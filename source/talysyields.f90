subroutine talysyields
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read yields from TALYS
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
! Constants
!   parsym    ! symbol of particle
! Variables for yields in ENDF format
!   xsprod    ! particle production cross section
!   yieldp    ! particle production yield
!
! *** Declaration of local data
!
  implicit none
  logical          :: lexist       ! logical to determine existence
  character(len=9) :: prodfile     ! file with total particle production cross sections
  character(len=132) :: line       ! 
  character(len=132) :: key
  integer          :: N            ! neutron number of residual nucleus
  integer          :: nin          ! counter for incident energy
  integer          :: type         ! particle type
  integer          :: istat        ! 
  integer          :: keyix
!
! ****************** Read total production yields **********************
!
! 1. Production cross sections and particle yields
!
  do type = 0, 6
    prodfile = ' prod.tot'
    write(prodfile(1:1), '(a1)') parsym(type)
    inquire (file = prodfile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = prodfile, status = 'old')
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) N
          if (istat /= 0) call read_error(prodfile, istat)
          read(1,'(/)')
          do nin = 1, N
            read(1, '(15x, 2es15.6)', iostat=istat) xsprod(type, nin), yieldp(type, nin)
            if (istat == -1) exit
            if (istat > 0) call read_error(prodfile, istat)
          enddo
          exit
        endif
      enddo
      close (unit = 1)
    endif
  enddo
  return
end subroutine talysyields
! Copyright A.J. Koning 2021
