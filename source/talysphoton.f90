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
  character(len=6)  :: resstring       ! ZA string
  character(len=30) :: decayfile       ! decay data file
  character(len=30) :: discretefile    ! file with discrete level data
  character(len=132) :: line           !
  character(len=132) :: key           !
  integer           :: type            ! particle type
  integer           :: i               ! counter
  integer           :: keyix
  integer           :: j               ! counter
  integer           :: istat           !
  integer           :: Z
  integer           :: A
!
! ********************* Read gamma decay schemes ***********************
!
  do type = 0, 6
    Z = Zinit - parZ(type)
    A = Ainit - parZ(type) - parN(type)
    write(resstring(1:3),'(i3.3)') Z
    write(resstring(4:6),'(i3.3)') A
    decayfile = 'gamma'//resstring//'.tot'
    inquire (file = decayfile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = decayfile, status = 'old')
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(1,'(/)')
          do
            read(1, '(i6, 9x, 15x, i6, 9x, es15.6)', iostat = istat) i, Ngamdis(type, i), yieldg(type, i)
            if (istat == -1) exit
            do j = 1, Ngamdis(type, i)
              read(1, '(90x, 2es15.6)') Egamdis(type, i, j), yieldratio(type, i, j)
            enddo
          enddo
          if (i == numlevels) exit
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
    Z = Zinit - parZ(type)
    A = Ainit - parZ(type) - parN(type)
    write(resstring(1:3),'(i3.3)') Z
    write(resstring(4:6),'(i3.3)') A
    discretefile = 'levels'//resstring//'.tot'
    inquire (file = discretefile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = discretefile, status = 'old')
      do
        read(1,'(a)', iostat = istat) line
        if (istat == -1) exit
        key='number of excited levels'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nlev(type)
        if (istat /= 0) call read_error(discretefile, istat)
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
!
! Skip ground state
!
          read(1,'(//)')
          do i= 1, nlev(type)
            read(1, '(t65, i6)') Nbranch(type, i)
            do j = 1, Nbranch(type, i)
              read(1, '(t80, i6, 2x, f15.4)') branchlevel(type, i, j), branchratio(type, i, j)
            enddo
          enddo
          if (i == numlevels) exit
          exit
        endif
      enddo
      close (unit = 1)
    endif
  enddo
  return
end subroutine talysphoton
! Copyright A.J. Koning 2023
