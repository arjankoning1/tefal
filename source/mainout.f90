subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2025-10-12: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
  write(*,'(/"    TEFAL-2.14 (Version: October 12, 2025)"/)')
  write(*, '(10x, " Creating ENDF-6 files with TALYS")')
  write(*, '(/" Copyright (C) 2025  A.J. Koning     ")')
!
! ***************** Write input file and default parameters ************
!
! inputout: subroutine to write the input parameters
!
  call inputout
  return
end subroutine mainout
! Copyright A.J. Koning 2025
