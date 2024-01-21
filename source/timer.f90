subroutine timer
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of execution time
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
!   sgl          ! single precision kind
! Variables for input of ENDF library type
!   endffile     ! name of ENDF file
!   flagclean    ! flag to clean up double points
! Variables for initialization of ENDF format
!   iclean       ! number of cleaned points
!
! *** Declaration of local data
!
  implicit none
  integer   :: hour                 ! number of hours
  integer   :: hundred              ! number of 1/100th of seconds
  integer   :: minute               ! number of minutes
  integer   :: second               ! number of seconds
  real(sgl) :: etime                ! time function
  real(sgl) :: tarray(2)            ! help variable
  real(sgl) :: time                 ! time
!
! ****** Get elapsed time in seconds from beginning of execution *******
!
! etime    : time function
!
! The returned time should be "charge time" (e.g., cp+pp+sys). This
! could be machine dependent.
!
  if (flagclean) then
    write(8, '(/" Number of deleted points:", i6)') iclean
    close (unit = 8)
    write(*, '(/" Number of deleted points:", i6, " (see file tefal.clean)")') iclean
  endif
  time = etime(tarray)
  hour = int(time / 3600.)
  minute = int((time - hour * 3600) / 60.)
  second = int(time - hour * 3600 - minute * 60)
  hundred = int(100 * (time - int(time)))
  write(*, '(/" Execution time:", i3, " hours ", i2, " minutes ", i2, ".", i2.2, " seconds ")') hour, minute, second, hundred
  write(*, '(/" TEFAL has created the following ENDF-6 file: ", a)') trim(endffile)
  return
end subroutine timer
! Copyright A.J. Koning 2021
