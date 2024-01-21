subroutine talysfission
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read fission data from TALYS
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
! Variables for input of ENDF structure
!   flagmulti      ! flag to include multi - chance fission
! Variables for partial cross sections in ENDF format
!   flagfission    ! flag for fission
!   idchannel      ! identifier for channel
!   reacstring     ! string with reaction information
!   xsexcl         ! exclusive cross section
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist      ! logical to determine existence
  character(len=11) :: fisfile     ! fission file
  character(len=12) :: fiscfile    ! fission file
  character(len=132) :: line       ! 
  character(len=132) :: source
  character(len=132) :: title
  character(len=132) :: key
  integer           :: idc         ! help variable
  integer           :: in          ! counter for neutrons
  integer           :: istat       ! logical for file existence
  integer           :: keyix
  integer           :: nen         ! energy counter
  integer           :: nin         ! counter for incident energy
!
! ******************** Read fission cross sections ********************
!
  idchannel(-1) = -99
  flagfission = .false.
  fisfile = 'fission.tot'
  inquire (file = fisfile, exist = lexist)
  if (lexist) then
    open (unit = 1, file = fisfile, status = 'old')
    do
      read(1,'(a)', iostat = istat) line
      if  (istat == -1) exit
      key='title'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) title
      key='source'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) source
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nen
        if (istat /= 0) call read_error(fisfile, istat)
        read(1,'(/)')
        do nin = 1, nen
          read(1, '(15x, es15.6)') xsexcl(-1, nin)
        enddo
        flagfission = .true.
        exit
      endif
    enddo
    close (unit = 1)
    reacstring(-1, 0) = ' '//trim(source)//' '//title(1:len_trim(title))
  endif
!
! ***************** Multichance fission cross sections *****************
!
  if (flagfission .and. flagmulti) then
    do in = 0, 3
      idc = -2 - in
      idchannel(idc) = -98 + in
      fiscfile = 'xs000000.fis'
      write(fiscfile(3:3), '(i1)') in
      inquire (file = fiscfile, exist = lexist)
      if (lexist) then
        open (unit = 1, file = fiscfile, status = 'old')
        do
          read(1,'(a)', iostat = istat) line
          if  (istat == -1) exit
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nen
            if (istat /= 0) call read_error(fiscfile, istat)
            read(1,'(/)')
            do nin = 1, nen
              read(1, '(15x, es15.6)', iostat=istat) xsexcl(idc, nin)
              if  (istat == -1) exit
            enddo
            exit
          endif
        enddo
        close (unit = 1)
      endif
    enddo
  endif
  return
end subroutine talysfission
! Copyright A.J. Koning 2021
