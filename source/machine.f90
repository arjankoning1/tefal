subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Machine dependent statements
!
! Author    : Arjan Koning
!
! 2023-10-27: Original code
! 2026-08-25: Runtime definition of TEFAL directory
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! Variables for path names
!   month      ! month
!   path       ! directory containing files to be read
!   year       ! year
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist          ! logical to determine existence
  character(len=3)   :: monthC(12)      ! month
  character(len=1024):: code_dir        ! code directory
  character(len=1024):: tefal_dir       ! code directory runtime defined
  integer            :: envstat
  integer            :: i               ! counter
  integer            :: n
  integer            :: values(8)       ! date and time values
!
! ********************* Set routine for date ***************************
!
!
  monthC = (/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
  call date_and_time(VALUES = values)
  month = monthC(values(2))
  year = mod(values(1), 100)
!
! ************************ Set directories *****************************
!
! The preferred option is to set an environment variable TEFAL_DIR,
! e.g. put in your ~/.profile or ~/.zshrc file.
!
! export TEFAL_DIR=/path/to/tefal
!
! If TEFAL_DIR is not set, get_environment_variable will simply return
! an empty string.
!
  call get_environment_variable('TEFAL_DIR', tefal_dir, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    code_dir = trim(tefal_dir)
  else
!
! If for some reason the above does not work, the code directory can be
! changed here manually.
!
    code_dir = '/path/to/tefal/'
  endif
  i = len_trim(code_dir)
  if (i > 0) then
    if (code_dir(i:i) /= '/') code_dir = trim(code_dir)//'/'
  endif
  path = code_dir
!
! Test to check accessibility of TEFAL data files
!
  inquire (file = trim(path)//'misc/endf_n.txt', exist = lexist)
  if (.not. lexist) then
    write(*, '(a)') 'TEFAL error: misc database not found.'
    write(*, '(2a)') 'Expected file: ', trim(path)//'misc/endf_n.txt'
    write(*, '(a)') 'Set the TEFAL_DIR environment variable:'
    write(*, '(a)') '  export TEFAL_DIR=/path/to/tefal'
    write(*, '(a)') 'Alternatively, edit code_dir in source/machine.f90'
    write(*, '(a)') 'and rebuild TEFAL.'
    error stop 77
  endif
  return
end subroutine machine
! Copyright A.J. Koning 2026
