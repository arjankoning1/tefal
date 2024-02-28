subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  write(*,'(/"    TEFAL-2.01 (Version: February 25, 2024)"/)')
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
! Copyright A.J. Koning 2023
