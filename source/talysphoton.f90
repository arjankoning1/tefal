subroutine talysphoton
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read gamma decay scheme from TALYS
!
! Author    : Arjan Koning
!
! 2023-04-26: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
  use A1_error_handling_mod
!
! All global variables
!   numlevels      ! maximum number of discrete levels
! Constants
!   parsym         ! symbol of particle
! Variables for photon production in ENDF format
!   branchlevel    ! level to which branching takes place
!   branchratio    ! branch ratio
!   Egamdis        ! energy of gamma ray
!   Nbranch        ! number of branches for level
!   Ngamdis        ! number of gamma ray lines per level
!   nlev           ! number of excited levels for nucleus
!   yieldg         ! total discrete gamma yield per level
!   yieldratio     ! yield ratio for level
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist          ! logical to determine existence
  character(len=7)  :: decayfile       ! decay data file
  character(len=10) :: discretefile    ! file with discrete level data
  character(len=132) :: line           !
  character(len=132) :: key           !
  integer           :: type            ! particle type
  integer           :: i               ! counter
  integer           :: keyix
  integer           :: j               ! counter
  integer           :: NNL             ! number of levels
  integer           :: istat           !
!
! ********************* Read gamma decay schemes ***********************
!
  do type = 0, 6
    decayfile = 'decay.'//parsym(type)
    inquire (file = decayfile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = decayfile, status = 'old')
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) NNL
          if (istat /= 0) call read_error(decayfile, istat)
          read(1,'(/)')
          do i = 1, min(NNL, numlevels)
            read(1, '(15x, i6, 9x, es15.6)') Ngamdis(type, i), yieldg(type, i)
            do j = 1, Ngamdis(type, i)
              read(1, '(2es15.6)') Egamdis(type, i, j), yieldratio(type, i, j)
            enddo
          enddo
          exit
        endif
      enddo
      close (unit = 1)
    endif
  enddo
!
! ************** Read gamma ray transition probabilities ***************
!
  do type = 0, 6
    discretefile = 'discrete.'//parsym(type)
    inquire (file = discretefile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = discretefile, status = 'old')
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nlev(type)
          if (istat /= 0) call read_error(discretefile, istat)
          read(1,'(/)')
          do i = 0, min(nlev(type), numlevels)
            read(1, '(35x, i2)') Nbranch(type, i)
            do j = 1, Nbranch(type, i)
              read(1, '(37x, i3, f10.4)') branchlevel(type, i, j), branchratio(type, i, j)
            enddo
          enddo
          exit
        endif
      enddo
      close (unit = 1)
    endif
  enddo
  return
end subroutine talysphoton
! Copyright A.J. Koning 2023
