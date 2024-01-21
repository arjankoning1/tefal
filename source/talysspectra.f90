subroutine talysspectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read spectra and pre-equilibrium ratios from TALYS
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
!   sgl            ! single precision kind
! All global variables
!   numen2         ! number of emission energies
! Constants
!   parsym         ! symbol of particle
!   pi             ! pi
! Variables for reaction initialization
!   Especindex     ! enegy index for spectra
!   Nenspec        ! number of incident energies for spectra
! Variables for TALYS info
!   eninc          ! incident energy
!   flagblock      ! flag to block spectra, angle and damage files
! Variables for input of ENDF library type
!   flagbreakup    ! breakup flag
! Variables for input of ENDF structure
!   flagtabddx     ! flag to give explicit DDX in MF6
! Variables for spectra in ENDF format
!   buratio        ! break - up ratio
!   ddxemis        ! double - differential emission spectra
!   Eocum          ! emission energies for total production spectra
!   Eoddx          ! emission energies for double - differential spectra
!   ncumddx        ! number of emission energies for DDX spectra
!   ncumout        ! number of emission energies for total production spectra
!   Nddx           ! number of angles
!   preeqratio     ! pre - equilibrium ratio
!   rmuddx         ! cosine of angle
!   xsemis         ! total production emission spectra
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist      ! logical to determine existence
  character(len=17) :: specfile    ! file with composite particle spectra
  character(len=23) :: ddxfile     ! file for DDX
  character(len=132) :: line       !
  character(len=132) :: key        !
  integer           :: iang        ! running variable for angle
  integer           :: iddx        ! counter
  integer           :: istat       ! logical for file existence
  integer           :: nen         ! energy counter
  integer           :: nen2        ! energy counter
  integer           :: Nfile       ! help variable
  integer           :: type        ! particle type
  integer           :: keyix          !
  real(sgl)         :: ang         ! angle
  real(sgl)         :: Ein         ! incident energy
  real(sgl)         :: Efile       ! incident energy
  real(sgl)         :: Afile       ! angle
!
! ******************* Read spectra and pre-equilibrium ratios **********
!
! 1. Spectra and pre-equilibrium ratios
!
  do type = 0, 6
    do nen = 1, Nenspec
      Ein = eninc(Especindex(nen))
      if (flagblock) then
        specfile = ' spec.tot'
      else
        specfile = ' spec0000.000.tot'
        write(specfile(6:13), '(f8.3)') Ein
        write(specfile(6:9), '(i4.4)') int(Ein)
      endif
      write(specfile(1:1), '(a1)') parsym(type)
      inquire (file=specfile, exist=lexist)
      if (lexist) then
        open (unit=1, file=specfile, status='old')
        do
          read(1,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(specfile, istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(specfile, istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nfile
            if (istat /= 0) call read_error(specfile, istat)
            read(1,'(/)')
            if (flagblock) then
              if (abs(Ein-Efile) <= 1.e-4) then
                ncumout(type,nen) = min(Nfile,numen2-1)
                if (flagbreakup) then
                  do nen2 = 1, ncumout(type,nen)
                    read(1, '(2es15.6,60x,2es15.6)') Eocum(type,nen,nen2) , xsemis(type,nen,nen2), &
 &                    preeqratio(type,nen,nen2), buratio(type,nen,nen2)
                  enddo
                else
                  do nen2 = 1, ncumout(type,nen)
                    read(1, '(2es15.6,60x,es15.6)') Eocum(type,nen,nen2), xsemis(type,nen,nen2), preeqratio(type,nen,nen2)
                  enddo
                endif
                exit
              else
                do nen2 = 1, Nfile
                  read(1, '()')
                enddo
                cycle
              endif
            else
              ncumout(type,nen) = min(Nfile,numen2-1)
              if (flagbreakup) then
                do nen2 = 1, ncumout(type,nen)
                  read(1, '(2es15.6,60x,2es15.6)') Eocum(type,nen,nen2), xsemis(type,nen,nen2), &
 &                  preeqratio(type,nen,nen2), buratio(type,nen,nen2)
                enddo
              else
                do nen2 = 1, ncumout(type,nen)
                  read(1, '(2es15.6,60x,es15.6)') Eocum(type,nen,nen2), xsemis(type,nen,nen2), preeqratio(type,nen,nen2)
                enddo
              endif
              exit
            endif
          endif
        enddo
        close (unit=1)
      endif
    enddo
  enddo
!
! 2. For incident deuterons: Double-differential nucleon spectra
!
  if (flagtabddx) then
    Nddx = 0
    do type = 1, 2
      do nen = 1, Nenspec
        Ein = eninc(Especindex(nen))
        iddx=0
        do iang = 0, 180
          ang = real(iang)
          if (flagblock) then
            ddxfile = ' ddx.deg'
          else
            ddxfile = ' ddxE0000.000A000.0.deg'
            write(ddxfile(6:13), '(f8.3)') Ein
            write(ddxfile(6:9), '(i4.4)') int(Ein)
            write(ddxfile(15:19), '(f5.1)') ang
            write(ddxfile(15:17), '(i3.3)') int(ang)
          endif
          write(ddxfile(1:1), '(a1)') parsym(type)
          inquire (file=ddxfile, exist=lexist)
          if (lexist) then
            open (unit=1, file=ddxfile, status='old')
            do
              read(1,'(a)', iostat = istat) line
              if (istat == -1) exit
              if (istat > 0) call read_error(ddxfile, istat)
              key='E-incident [MeV]'
              keyix=index(line,trim(key))
              if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
              if (istat /= 0) call read_error(ddxfile, istat)
              key='Angle [deg]'
              keyix=index(line,trim(key))
              if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Afile
              if (istat /= 0) call read_error(ddxfile, istat)
              key='entries'
              keyix=index(line,trim(key))
              if (keyix > 0) then
                read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nfile
                if (istat /= 0) call read_error(ddxfile, istat)
                read(1,'(/)')
                if (flagblock) then
                  if (abs(Ein-Efile) <= 1.e-4 .and. abs(ang-Afile) <= 1.e-4) then
                    iddx = iddx + 1
                    ncumddx(type,nen,iddx) = min(Nfile,numen2)
                    rmuddx(iddx) = cos(ang * pi / 180.)
                    do nen2 = 1, ncumddx(type,nen,iddx)
                      read(1, '(2es15.6)') Eoddx(type,nen,iddx,nen2), ddxemis(type,nen,iddx,nen2)
                    enddo
                    exit
                  else
                    do nen2 = 1, Nfile
                      read(1,'()')
                    enddo
                    cycle
                  endif
                else
                  iddx = iddx + 1
                  ncumddx(type,nen,iddx) = min(Nfile,numen2)
                  rmuddx(iddx) = cos(ang * pi / 180.)
                  do nen2 = 1, ncumddx(type,nen,iddx)
                    read(1, '(2es15.6)') Eoddx(type,nen,iddx,nen2), ddxemis(type,nen,iddx,nen2)
                  enddo
                  exit
                endif
              endif
            enddo
            close (unit = 1)
          endif
        enddo
        Nddx = iddx
      enddo
    enddo
  endif
  return
end subroutine talysspectra
! Copyright A.J. Koning 2021
