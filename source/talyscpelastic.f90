subroutine talyscpelastic
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read charged-particle elastic scattering data from TALYS
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
!   sgl        ! single precision kind
! Variables for TALYS info
!   eninc      ! incident energy
!   flagblock    ! flag to block spectra, angle and gamma files
!   k0         ! index of incident particle
!   Ltarget    ! excited level of target
!   numinc     ! number of incident energies
! Constants
!   parsym     ! symbol of particle
! Variables for angular distributions in ENDF format
!   cpang      ! differential cross section
!   elasni     ! nuclear+interference term
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist     ! logical to determine existence
  character(len=17) :: angfile    ! name of file with angular distributions
  character(len=132) :: line    ! 
  character(len=132) :: key
  integer           :: iang       ! running variable for angle
  integer           :: istat      ! logical for file access
  integer           :: keyix
  integer           :: nang1      ! number of angles
  integer           :: nin        ! counter for incident energy
  real(sgl)         :: Ein        ! incident energy
  real(sgl)         :: Efile      ! incident energy
!
! ******** Read charged-particle elastic angular distributions *********
!
  do nin = 1, numinc
    Ein = eninc(nin)
    if (flagblock) then
      angfile = '  ang.L00'
      write(angfile(8:9), '(i2.2)') Ltarget
    else
      angfile = '          ang.L00'
      write(angfile(3:10), '(f8.3)') Ein
      write(angfile(3:6), '(i4.4)') int(Ein)
      write(angfile(16:17), '(i2.2)') Ltarget
    endif
    write(angfile(1:2), '(2a1)') parsym(k0), parsym(k0)
    inquire (file=angfile, exist=lexist)
    if (lexist) then
      open (unit = 1, file = angfile, status = 'old')
      do
        read(1,'(a)',iostat = istat) line
        if (istat == -1) exit
        key='E-incident [MeV]'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
        if (istat /= 0) call read_error(angfile, istat)
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nang1
          if (istat /= 0) call read_error(angfile, istat)
          read(1,'(/)')
          if (flagblock) then
            if (abs(Ein-Efile) <= 1.e-4) then
              do iang = 0, nang1 - 1
                read(1,'(15x, es15.6, 45x, es15.6)') cpang(nin, iang), elasni(nin,iang)
              enddo
              exit
            else
              do iang = 1, nang1
                read(1, '()')
              enddo
              cycle
            endif
          else
            do iang = 0, nang1 - 1
              read(1,'(15x, es15.6, 45x, es15.6)') cpang(nin, iang), elasni(nin, iang)
            enddo
            exit
          endif
        endif
      enddo
      close (unit=1)
    endif
  enddo
  return
end subroutine talyscpelastic
