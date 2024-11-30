subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2024-11-21: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
  write(*,'(/"    TEFAL-2.1 (Version: November 30, 2024)"/)')
  write(*, '(10x, " Creating ENDF-6 files with TALYS")')
  write(*, '(/" Copyright (C) 2024  A.J. Koning     ")')
!
! ***************** Write input file and default parameters ************
!
! inputout: subroutine to write the input parameters
!
  call inputout
  return
end subroutine mainout
! Copyright A.J. Koning 2024
