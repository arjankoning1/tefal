subroutine countchannels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Count number of exclusive reaction channels from TALYS
!
! Author    : Arjan Koning
!
! 2026-09-03: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! All global variables
!   nummt          ! number of MT numbers
! Variables for input of ENDF structure
!   flagmtextra    ! flag to include extra MT numbers up to MT200
! Variables for initialization of ENDF format
!   idnum          ! number of different exclusive cross sections
!   MTnum          ! channel identifier for MT-number
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist         ! logical to determine existence
  character(len=12) :: xsfile         ! file with cross sections
  integer           :: ia             ! counter for alpha particles
  integer           :: iaend          ! end of particle summation
  integer           :: id             ! counter for deuterons
  integer           :: id0            ! channel identifier
  integer           :: idend          ! end of particle summation
  integer           :: ih             ! counter for helions
  integer           :: ihend          ! end of particle summation
  integer           :: imt            ! MT counter
  integer           :: in             ! counter for neutrons
  integer           :: inend          ! end of particle summation
  integer           :: ip             ! counter for protons
  integer           :: ipend          ! end of particle summation
  integer           :: it             ! counter for tritons
  integer           :: itend          ! end of particle summation
  integer           :: maxparticle    ! maximum number of ejectiles
  integer           :: npart          ! number of ejectiles
!
! ******************** Count exclusive channels ***********************
!
! Use exactly the same channel limits as in talyschannels.
!
  if (flagmtextra) then
    inend = 8
    ipend = 3
    maxparticle = 8
  else
    inend = 4
    ipend = 2
    maxparticle = 4
  endif
  idend = 1
  itend = 1
  ihend = 1
  iaend = 3
!
! idnum is the highest channel index. Since the first channel has
! index 0, start at -1.
!
  idnum = -1
  do ih = 0, ihend
    do it = 0, itend
      do id = 0, idend
        do ia = 0, iaend
          do ip = 0, ipend
            do in = 0, inend
              npart = in + ip + id + it + ih + ia
              if (npart > maxparticle) cycle
              id0 = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
              do imt = 1, nummt
                if (MTnum(imt) == id0) then
                  xsfile = 'xs      .tot'
                  write(xsfile(3:3), '(i1)') in
                  write(xsfile(4:4), '(i1)') ip
                  write(xsfile(5:5), '(i1)') id
                  write(xsfile(6:6), '(i1)') it
                  write(xsfile(7:7), '(i1)') ih
                  write(xsfile(8:8), '(i1)') ia
                  inquire (file = xsfile, exist = lexist)
                  if (lexist) idnum = idnum + 1
                  exit
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine countchannels
! Copyright A.J. Koning 2026
