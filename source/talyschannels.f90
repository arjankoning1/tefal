subroutine talyschannels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read exclusive channel data from TALYS
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
!   numenrec       ! number of incident energies for recoils
!   numlevels      ! maximum number of discrete levels
!   nummt          ! number of MT numbers
! Variables for TALYS info
!   eninc          ! incident energy
!   flagblock      ! flag to block spectra, angel and gamma files
! Variables for input of ENDF library type
!   flaggpf        ! flag for general purpose library
! Variables for input of ENDF structure
!   flagmtextra    ! flag to include extra MT numbers up to MT200
!   flagrecoil     ! flag to include recoil information
! Variables for reaction initialization
!   Egamindex      ! enegy index for gamma cross sections
!   Especindex     ! enegy index for spectra
!   Nengam         ! number of incident energies for gamna c.s.
!   Nenspec        ! number of incident energies for spectra
!   nlevmax        ! number of included discrete levels
! Variables for initialization of ENDF format
!   idnum          ! number of different exclusive cross sections
!   MTnum          ! channel identifier for MT - number (ENDF format)
! Variables for partial cross sections in ENDF format
!   Egammadis      ! gamma energy
!   Eout           ! emission energy
!   Erec           ! recoil energy
!   Estartdis      ! starting level
!   Ethexcl        ! threshold energy
!   Ethexcliso     ! threshold energy for isomer
!   idchannel      ! identifier for channel
!   isoexist       ! flag for existence of isomer
!   Nisomer        ! number of isomers
!   nout           ! number of emission energies
!   noutrec        ! number of recoil energies
!   Qexcl          ! Q - value
!   Qexcliso       ! Q - value for isomer
!   reacstring     ! string with reaction information
!   recexcl        ! exclusive recoils
!   specexcl       ! exclusive spectra
!   xsexcl         ! exclusive cross section
!   xsexcliso      ! exclusive cross section for isomer
!   xsgamdis       ! exclusive discrete gamma - ray cross section
!   xsgamexcl      ! exclusive gamma cross section
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist         ! logical to determine existence
  character(len=12) :: isofile        ! file with isomeric cross section
  character(len=12) :: gamfile        ! giant resonance parameter file
  character(len=12) :: xsfile         ! file with cross sections
  character(len=21) :: recfile        ! file with recoil spectra
  character(len=21) :: spfile         ! file with spectra
  character(len=132) :: title
  character(len=132) :: source
  character(len=132) :: line        ! 
  character(len=132) :: key        ! 
  integer           :: i1             ! value
  integer           :: i2             ! value
  integer           :: ia             ! counter for alpha particles
  integer           :: iaend          ! end of particle summation
  integer           :: id             ! counter for deuterons
  integer           :: id0            ! help variable
  integer           :: idc            ! help variable
  integer           :: idend          ! end of particle summation
  integer           :: igam           ! counter for gammas
  integer           :: ih             ! counter for helions
  integer           :: ihend          ! end of particle summation
  integer           :: imt            ! MT counter
  integer           :: in             ! counter for neutrons
  integer           :: inend          ! end of particle summation
  integer           :: ip             ! counter for protons
  integer           :: ipend          ! end of particle summation
  integer           :: istat          ! error code
  integer           :: it             ! counter for tritons
  integer           :: itend          ! end of particle summation
  integer           :: keyix          ! 
  integer           :: maxparticle    ! maximum number of ejectiles
  integer           :: N              ! neutron number of residual nucleus
  integer           :: nen            ! energy counter
  integer           :: nen2           ! energy counter
  integer           :: nex            ! discrete level
  integer           :: nex0           ! discrete level
  integer           :: Nfile          ! help variable
  integer           :: Ng             ! number of discrete gamma-rays
  integer           :: nin            ! counter for incident energy
  integer           :: npart          ! number of ejectiles
  integer           :: type           ! particle type
  real(sgl)         :: Efinal         ! final level
  real(sgl)         :: Ein            ! incident energy
  real(sgl)         :: eps
  real(sgl)         :: Erat
  real(sgl)         :: Efile          ! incident energy
  real(sgl)         :: Egam           ! gamma energy
  real(sgl)         :: Est            ! starting level
  real(sgl)         :: xsga           ! exclusive discrete gamma-ray cross section
!
! ****************** Read exclusive cross sections *********************
!
! We loop over all exclusive reaction channels that are possible within the ENDF-6 format.
!
! For the high-energy EAF format, extra MT-numbers with more particles are allowed.
! For the other cases, we stick to the official ENDF-6 format.
!
  eps = 1.e-6
  isoexist = .false.
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
  idc = -1
  idchannel(idc) = -99
  do ih = 0, ihend
    do it = 0, itend
      do id = 0, idend
        do ia = 0, iaend
          do ip = 0, ipend
            do in = 0, inend
!
! Various checks for existence of channel
!
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
!
! Read cross sections
!
                  inquire (file = xsfile, exist = lexist)
                  if (.not. lexist) cycle
                  open (unit = 1, file = xsfile, status = 'old', iostat = istat)
                  if (istat /= 0) call read_error(xsfile, istat)
                  idc = idc + 1
                  idchannel(idc) = id0
                  do
                    read(1,'(a)',iostat = istat) line
                    if (istat == -1) exit
                    key='title'
                    keyix=index(line,trim(key))
                    if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) title
                    if (istat /= 0) call read_error(xsfile, istat)
                    key='source'
                    keyix=index(line,trim(key))
                    if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) source
                    if (istat /= 0) call read_error(xsfile, istat)
                    key='Q-value [MeV]'
                    keyix=index(line,trim(key))
                    if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Qexcl(idc)
                    if (istat /= 0) call read_error(xsfile, istat)
                    key='E-threshold [MeV]'
                    keyix=index(line,trim(key))
                    if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ethexcl(idc)
                    if (istat /= 0) call read_error(xsfile, istat)
                    key='entries'
                    keyix=index(line,trim(key))
                    if (keyix > 0) then
                      read(line(keyix+len_trim(key)+2:80),*, iostat = istat) N
                      if (istat /= 0) call read_error(xsfile, istat)
                      read(1,'(/)')
                      do nin = 1, N
                        read(1, '(15x, 2es15.6)', iostat = istat) xsexcl(idc, nin), xsgamexcl(idc, nin)
                        if (istat == -1) exit
                        if (istat > 0) call read_error(xsfile, istat)
                      enddo
                      exit
                    endif
                  enddo
                  close (unit = 1)
                  reacstring(idc, 0) = ' '//trim(source)//' '//title(1:len_trim(title))
!
! Isomers
!
                  isofile(1:9) = xsfile(1:9)
                  isofile(10:10) = 'L'
                  Nisomer(idc) = 0
                  do nex0 = 0, numlevin
                    write(isofile(11:12), '(i2.2)') nex0
                    inquire (file = isofile, exist = lexist)
                    if (lexist) then
                      open (unit = 1, file = isofile, status = 'old', iostat = istat)
                      if (istat /= 0) call read_error(isofile, istat)
                      nex = min (nex0, nlevmax)
                      isoexist(idc, nex) = .true.
                      do
                        read(1,'(a)',iostat = istat) line
                        if (istat == -1) exit
                        key='title'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) title
                        key='source'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) source
                        key='Q-value [MeV]'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Qexcliso(idc,nex)
                        if (istat /= 0) call read_error(xsfile, istat)
                        key='E-threshold [MeV]'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ethexcliso(idc,nex)
                        if (istat /= 0) call read_error(xsfile, istat)
                        key='entries'
                        keyix=index(line,trim(key))
                        if (keyix > 0) then
                          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) N
                          if (istat /= 0) call read_error(xsfile, istat)
                          read(1,'(/)')
                          do nin = 1, N
                            read(1, '(15x, es15.6)', iostat = istat) xsexcliso(idc, nex, nin)
                            if (istat == -1) exit
                            if (istat > 0) call read_error(xsfile, istat)
                          enddo
                          exit
                        endif
                      enddo
                      close (unit = 1)
                      Nisomer(idc) = Nisomer(idc) + 1
                      reacstring(idc, nex) = ' '//trim(source)//' '//title(1:len_trim(title))
                    endif
                  enddo
!
! Exclusive discrete gamma-rays (general purpose files only)
!
                  if (flaggpf) then
                    gamfile = xsfile
                    gamfile(10:12) = 'gam'
                    inquire (file = gamfile, exist = lexist)
                    if (lexist) then
                      open (unit = 1, file = gamfile, status = 'old', iostat = istat)
                      if (istat /= 0) call read_error(gamfile, istat)
                      do
                        read(1,'(a)',iostat = istat) line
                        if (istat == -1) exit
                        key='E-incident [MeV]'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
                        if (istat /= 0) call read_error(gamfile, istat)
                        key='entries'
                        keyix=index(line,trim(key))
                        if (keyix > 0) then
                          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ng
                          if (istat /= 0) call read_error(gamfile, istat)
                          read(1,'(/)')
                          do nen = 1, Nengam
                            Ein = eninc(Egamindex(nen))
                            Erat = Efile/Ein
                            if (Erat >= 1. - eps .and. Erat <= 1. + eps) then
                              do igam = 1, Ng
                                read(1, '(2(i6, 9x), 4e15.6)', iostat = istat) i1, i2, xsga, Est, Efinal, Egam
                                if (istat /= 0) call read_error(gamfile, istat)
                                if (i1 <= numlevels .and. i2 <= numlevels) then
                                  xsgamdis(idc, nen, i1, i2) = xsga
                                  Estartdis(idc, i1, i2) = Est
                                  Egammadis(idc, i1, i2) = Egam
                                endif
                              enddo
                            endif
                          enddo
                        endif
                      enddo
                      close (unit = 1)
                    endif
!
! Exclusive spectra (general purpose files only)
!
                    spfile = 'sp      E0000.000.tot'
                    spfile(3:8) = xsfile(3:8)
                    do nen = 1, Nenspec
                      Ein = eninc(Especindex(nen))
                      nout(idc, nen) = 0
                      Eout(idc, nen, 0) = 0.
                      do type = 0, 6
                        specexcl(idc, nen, type, 0) = 0.
                      enddo
                      if (flagblock) then
                        spfile = 'sp      .tot'
                      else
                        spfile = 'sp      E0000.000.tot'
                        write(spfile(10:17), '(f8.3)') Ein
                        write(spfile(10:13), '(i4.4)') int(Ein)
                      endif
                      spfile(3:8) = xsfile(3:8)
                      inquire (file = spfile, exist = lexist)
                      if (lexist) then
                        open (unit=1, file=spfile, status='old')
                        do
                          read(1,'(a)', iostat = istat) line
                          if (istat == -1) exit
                          if (istat > 0) call read_error(spfile, istat)
                          key='E-incident [MeV]'
                          keyix=index(line,trim(key))
                          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
                          if (istat /= 0) call read_error(spfile, istat)
                          key='entries'
                          keyix=index(line,trim(key))
                          if (keyix > 0) then
                            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nfile
                            if (istat /= 0) call read_error(spfile, istat)
                            read(1,'(/)')
                            if (flagblock) then
                              Erat = Efile/Ein
                              if (Erat >= 1. - eps .and. Erat <= 1. + eps) then
                                nout(idc, nen) = min(Nfile, numen2-1)
                                do nen2 = 1, nout(idc,nen)
                                  read(1, '(8es15.6)', iostat = istat) Eout(idc,nen,nen2), (specexcl(idc,nen,type,nen2),type=0,6)
                                  if (istat > 0) call read_error(spfile, istat)
                                enddo
                                exit
                              else
                                do nen2 = 1, Nfile
                                  read(1,'()', iostat = istat)
                                  if (istat > 0) call read_error(spfile, istat)
                                enddo
                                cycle
                              endif
                            else
                              nout(idc, nen) = min(Nfile, numen2-1)
                              do nen2=1,nout(idc,nen)
                                read(1,'(8es15.6)', iostat = istat) Eout(idc,nen,nen2), (specexcl(idc,nen,type,nen2),type=0,6)
                                if (istat > 0) call read_error(spfile, istat)
                              enddo
                              exit
                            endif
                          endif
                        enddo
                        close (unit=1)
                      endif
                    enddo
!
! Exclusive recoils (general purpose files only)
!
                    if (flagrecoil) then
                      do nen = 1, Nenspec
                        Ein = eninc(Especindex(nen))
                        noutrec(idc,nen) = 0
                        if (flagblock) then
                          recfile = 'sp      .rec'
                        else
                          recfile = 'sp      E0000.000.rec'
                          write(recfile(10:17), '(f8.3)') Ein
                          write(recfile(10:13), '(i4.4)') int(Ein)
                        endif
                        recfile(3:8) = xsfile(3:8)
                        inquire (file=recfile, exist=lexist)
                        if (lexist) then
                          open (unit=1, file=recfile, status='old')
                          do
                            read(1,'(a)', iostat = istat) line
                            if (istat == -1) exit
                            if (istat > 0) call read_error(recfile, istat)
                            key='E-incident [MeV]'
                            keyix=index(line,trim(key))
                            if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
                            if (istat /= 0) call read_error(recfile, istat)
                            key='entries'
                            keyix=index(line,trim(key))
                            if (keyix > 0) then
                              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nfile
                              if (istat /= 0) call read_error(recfile, istat)
                              read(1,'(/)')
                              if (flagblock) then
                                Erat = Efile/Ein
                                if (Erat >= 1. - eps .and. Erat <= 1. + eps) then
                                  noutrec(idc,nen) = min(Nfile,numenrec)
                                  do nen2 = 1,noutrec(idc,nen)
                                    read(1, '(2es15.6)', iostat = istat) Erec(idc,nen,nen2), recexcl(idc,nen,nen2)
                                    if (istat > 0) call read_error(recfile, istat)
                                  enddo
                                  exit
                                else
                                  do nen2 = 1, Nfile
                                    read(1, '()', iostat = istat)
                                    if (istat > 0) call read_error(recfile, istat)
                                  enddo
                                  cycle
                                endif
                              else
                                noutrec(idc,nen) = min(Nfile,numenrec)
                                do nen2 = 1, noutrec(idc,nen)
                                  read(1, '(2es15.6)', iostat = istat) Erec(idc,nen,nen2), recexcl(idc,nen,nen2)
                                  if (istat > 0) call read_error(recfile, istat)
                                enddo
                                exit
                              endif
                            endif
                          enddo
                          close (unit=1)
                        endif
                      enddo
                    endif
                  endif
                  exit
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
!
! Count the channels
!
  idnum = idc
  return
end subroutine talyschannels
! Copyright A.J. Koning 2021
