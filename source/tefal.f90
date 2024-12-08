program tefal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main program
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2024-12-08: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
!
!   |-------------------------------------------------------|
!   |                 Arjan Koning                          |
!   |                                                       |
!   | Email: A.Koning@@iaea.org                             |
!   |-------------------------------------------------------|
!
! MIT License
!
! Copyright (c) 2024 Arjan Koning
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! Input : - user input
!   - various files created by TALYS
!
! ********** Input, initialization and creation of ENDF-6 file *********
!
! machine       : subroutine for machine dependent statements
! tefalconstants: subroutine for constants and basic properties of particles
! tefalinput    : subroutine for input
! tefalinitial  : subroutine for initialization
! talysfiles    : subroutine to read data files from TALYS
! process       : subroutine to process data from TALYS
! tefalmake     : subroutine for construction of ENDF-6 data
!
  call machine
  call tefalconstants
  call tefalinput
  call tefalinitial
  call talysfiles
  call process
  call tefalmake
end program tefal
! Copyright A.J. Koning 2024
