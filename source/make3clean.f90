subroutine make3clean
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Clean up MF3
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tefal_mod
!
! All global variables
!   nummt      ! number of MT numbers
! Variables for info from TALYS
!   k0         ! index of incident particle
! Variables for reaction initialization
!   EHres      ! upper energy in resonance range
! Variables for initialization of ENDF format
!   iclean     ! number of cleaned points
!   mtexist    ! flag for existence of MT - number
! Variables for ENDF format
!   NBT        ! separation value for interpolation scheme
!   NP         ! number of incident energies
! Variables for MF3
!   E3         ! incident energy for MF3 (in ENDF - 6 format)
!   xs         ! cross section
!
! *** Declaration of local data
!
  implicit none
  logical :: flagdel              ! flag to delete a line
  integer :: i                    ! counter
  integer :: j                    ! counter
  integer :: MF                   ! MF-number
  integer :: MT                   ! MT-number
  integer :: MT2                  ! MT number
  integer :: MTb                  ! MT number
  integer :: N                    ! neutron number of residual nucleus
!
! Many actions may be required to construct MF3: linearization,
! adoption from other libraries, taking care of the overlap with the resonance range etc. (i.e. setting cross sections to zero).
! To make sure we get a "clean" MF3 we do some final checking and take actions if necessary.
! Of course, we report such actions in the output file.
!
! ****************** Remove redundant points ***************************
!
  open (unit = 8, file = 'tefal.clean', status = 'replace')
  write(8, '(/, 20x, " Deleted points with zero or equal values"/)')
  MF = 3
  do MT = 1, nummt
    if ( .not. mtexist(MF, MT)) cycle
    N = NP(MF, MT)
!
! Delete decreasing energies
!
    do
      flagdel = .false.
      do i = 2, N
        if (E3(MT, i) < E3(MT, i - 1)) then
          do j = i, N - 1
            E3(MT, j) = E3(MT, j + 1)
            xs(MT, j) = xs(MT, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted point: MF=", i3, " MT=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &          MF, MT, E3(MT, i+1), xs(MT, i + 1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if ( .not. flagdel) exit
      N = N - 1
    enddo
!
! Delete multiple ( > 2) points with equal incident energy
!
    do
      flagdel = .false.
      do i = 1, N - 2
        if (E3(MT, i) >= E3(MT, i + 1) .and. E3(MT, i) == E3(MT, i + 2)) then
          do j = i + 1, N - 1
            E3(MT, j) = E3(MT, j + 1)
            xs(MT, j) = xs(MT, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted point: MF=", i3, " MT=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &          MF, MT, E3(MT, i+1), xs(MT, i + 1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if ( .not. flagdel) exit
      N = N - 1
    enddo
!
! Delete multiple ( > 2) points with equal cross section or double points with equal energies and cross sections
!
    do
      flagdel = .false.
      do i = 1, N - 2
        if ((xs(MT, i) == xs(MT, i + 1) .and. xs(MT, i) == xs(MT, i + 2)) .or. &
          (E3(MT, i) == E3(MT, i + 1) .and. xs(MT, i) == xs(MT, i + 1))) then
          do j = i + 1, N - 1
            E3(MT, j) = E3(MT, j + 1)
            xs(MT, j) = xs(MT, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted point: MF=", i3, " MT=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &          MF, MT, E3(MT, i+1), xs(MT, i + 1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if ( .not. flagdel) exit
      N = N - 1
    enddo
!
! Delete second energy point if it has a zero cross section, for negative Q-value reactions. do this only above the RRR.
!
    flagdel = .false.
    if (E3(MT, 1) > EHres .and. xs(MT, 1) == 0..and.xs(MT, 2) == 0. .and. MT /= 5) then
      do j = 2, N - 1
        E3(MT, j) = E3(MT, j + 1)
        xs(MT, j) = xs(MT, j + 1)
      enddo
      flagdel = .true.
      write(8, '(" Deleted point: MF=", i3, " MT=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &      MF, MT, E3(MT, 2), xs(MT, 2)
      iclean = iclean + 1
    endif
    if (flagdel) N = N - 1
!
! Set new number of points
!
    NP(MF, MT) = N
    NBT(MF, MT, 1) = N
!
! Check for, and remove, negative cross sections
!
    if (k0 <= 1) then
      do i = 1, N
        if (xs(MT, i) < 0.) then
          if ((MT /= 1 .and. MT /= 2 .and. MT /= 18 .and. MT /= 102) .or. E3(MT, i) > 1.e6) then
            if (i > 1) then
              xs(MT, i) = xs(MT, i - 1)
            else
              xs(MT, i) = 0.
            endif
          endif
        endif
      enddo
    endif
  enddo
!
! If (n,p) (or (n,alpha) etc.) is explicitly set to zero, so are all the partials
!
  do MT = 103, 107
    if ( .not. mtexist(MF, MT)) then
      MTb = 600 + (MT - 103) * 50
      do MT2 = MTb, MTb + 40
        mtexist(MF, MT2) = .false.
      enddo
      mtexist(MF, MTb + 49) = .false.
    endif
  enddo
  return
end subroutine make3clean
! Copyright A.J. Koning 2021
