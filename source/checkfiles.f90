subroutine checkfiles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in external data files
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-03-29: Turned non-existing MF1, MF2 and MF5 into warnings
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! All global variables
!   nummf         ! number of MF numbers
!   nummt         ! number of MT numbers
! Variables for input of ENDF MF1
!   endftext      ! file with MF1 information
! Variables for input of ENDF library type
!   flageaf       ! flag for EAF - formatted activation library
!   flaggpf       ! flag for general purpose library
! Variables for input of ENDF structure
!   flagres       ! flag to include resonance parameters
! Variables for input of specific ENDF data
!   adopt         ! logical for existence of MF information (per MT)
!   adoptfile     ! name of library for MF information (per MT)
!   background    ! file with background cross sections
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist       ! logical to determine existence
  integer            :: i            ! counter
  integer            :: iend         ! counter
  integer            :: imf          ! MF counter
  integer            :: imt          ! MT counter
  character(len=132) :: filename     ! filename
!
! ******************* Check for data files to be adopted ***************
!
  do imf = 1, nummf
    do imt = 1, nummt
!
! ***** Check of existence of adopt files for various input options ****
!
      if (imf == 2 .and. imt == 151 .and. flaggpf) then
        iend = 2
      else
        iend = 1
      endif
      do i = 1, iend
        if (i == 1) then
          filename = adoptfile(imf, imt)
        else
          if (background(1:1) == ' ') then
            filename = adoptfile(imf, imt)
          else
            filename = background
          endif
        endif
        if (adopt(imf, imt)) then
          inquire (file = filename, exist = lexist)
          if (lexist) exit
          if (imf == 1 .and. (imt == 452 .or. imt == 455 .or. imt == 456 .or. imt == 458)) then
            write(*, '(" TEFAL-warning: ",a, " not present")') trim(filename)
            write(*,'("   MF ",i3," MT ",i3," omitted")')  imf,imt
            adopt(imf, imt) = .false.
            cycle
          endif
          if (imf == 2 .and. imt == 151) then
            write(*, '(" TEFAL-warning: ",a, " not present")') trim(filename)
            write(*,'("   MF ",i3," MT ",i3," potential scattering radius used")')  imf, imt
            adopt(imf, imt) = .false.
            cycle
          endif
          if (imf == 5 .and. (imt == 18 .or. imt == 455)) then
            write(*, '(" TEFAL-warning: ",a, " not present")') trim(filename)
            write(*,'("   MF ",i3," MT ",i3," omitted")')  imf,imt
            adopt(imf, imt) = .false.
            cycle
          endif
          write(*,'(" TEFAL-error: ",a, " not present")') trim(filename)
          stop
        endif
      enddo
!
! ***************** Check for existence of MF/MT file ******************
!
      if ( .not. adopt(imf, imt)) cycle
    enddo
  enddo
!
! Check or set default for ENDF-6 text files
!
  if (endftext(1:1) == ' ') then
    if (flaggpf) then
      if (k0 == 0) endftext = trim(path)//'misc/endf_g.txt'
      if (k0 == 1) endftext = trim(path)//'misc/endf_n.txt'
      if (k0 > 1) endftext = trim(path)//'misc/endf_p.txt'
    else
      if ( .not. flaggpf .and. .not. flageaf) then
        if (k0 == 0) endftext = trim(path)//'misc/endf_act_g.txt'
        if (k0 == 1) endftext = trim(path)//'misc/endf_act_n.txt'
        if (k0 > 1) endftext = trim(path)//'misc/endf_act_p.txt'
      endif
    endif
  else
    inquire (file = endftext, exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TEFAL-error: Non-existent MF1 file: ", a132)') endftext
      stop
    endif
  endif
  return
end subroutine checkfiles
! Copyright A.J. Koning 2021
