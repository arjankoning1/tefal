subroutine make9_10clean
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Clean up MF9 and MF10
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
! Variables for initialization of ENDF format
!   iclean     ! number of cleaned points
!   mtexist    ! flag for existence of MT - number
! Variables for MF8_10
!   E10        ! incident energy (in ENDF - 6 format)
!   E10ZA      ! incident energy (in ENDF - 6 format)
!   NBT10      ! separation value for interpolation scheme
!   NBTZA      ! separation value for interpolation scheme
!   NP10       ! number of incident energies
!   NPZA       ! number of incident energies
!   NSt        ! number of final states
!   xsiso      ! cross section for isomer (in ENDF - 6 format)
!   xsrpZA     ! cross section for residual production (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  logical :: flagdel              ! flag to delete a line
  integer :: i                    ! counter
  integer :: iso                  ! counter for isomer
  integer :: iza                  ! counter for Z,A combinations
  integer :: j                    ! counter
  integer :: MF                   ! MF-number
  integer :: MT                   ! MT-number
  integer :: N                    ! neutron number of residual nucleus
!
! To make sure we get a "clean" MF9-10 we do some final checking and take actions if necessary.
! Of course, we report such actions in the output file
!
! ****************** Remove redundant points ***************************
!
  do MF = 9, 10
    do MT = 1, nummt
      if ( .not. mtexist(MF, MT)) cycle
      if (MF == 10 .and. MT == 5) then
        do iza = 1, NSt(MT)
          N = NPZA(iza)
!
! Delete decreasing energies
!
          do
            flagdel = .false.
            do i = 2, N
              if (E10ZA(iza, i) < E10ZA(iza, i - 1)) then
                do j = i, N - 1
                  E10ZA(iza, j) = E10ZA(iza, j + 1)
                  xsrpZA(iza, j) = xsrpZA(iza, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iso=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iza, E10ZA(iza, i + 1), xsrpZA(iza, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Delete multiple ( > 2) points with equal incident energy
!
          do
            flagdel = .false.
            do i = 1, N - 2
              if (E10ZA(iza, i) >= E10ZA(iza, i + 1) .and. E10ZA(iza, i) == E10ZA(iza, i + 2)) then
                do j = i + 1, N - 1
                  E10ZA(iza, j) = E10ZA(iza, j + 1)
                  xsrpZA(iza, j) = xsrpZA(iza, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iso=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iza, E10ZA(iza, i + 1), xsrpZA(iza, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Delete multiple ( > 2) points with equal cross section or double points with equal energies and cross sections
!
          do
            flagdel = .false.
            do i = 1, N - 2
              if ((xsrpZA(iza, i) == xsrpZA(iza, i + 1) .and. xsrpZA(iza, i) == xsrpZA(iza, i + 2)) .or. &
                (E10ZA(iza, i) == E10ZA(iza, i + 1) .and. xsrpZA(iza, i) == xsrpZA(iza, i + 1))) then
                do j = i + 1, N - 1
                  E10ZA(iza, j) = E10ZA(iza, j + 1)
                  xsrpZA(iza, j) = xsrpZA(iza, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iza=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iza, E10ZA(iza, i + 1), xsrpZA(iza, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Set new number of points
!
          NPZA(iza) = N
          NBTZA(iza, 1) = N
        enddo
      else
        do iso = 1, NSt(MT)
          N = NP10(MT, iso)
!
! Delete decreasing energies
!
          do
            flagdel = .false.
            do i = 2, N
              if (E10(MT, iso, i) < E10(MT, iso, i - 1)) then
                do j = i, N - 1
                  E10(MT, iso, j) = E10(MT, iso, j + 1)
                  xsiso(MT, iso, j) = xsiso(MT, iso, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iso=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iso, E10(MT, iso, i + 1), xsiso(MT, iso, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Delete multiple ( > 2) points with equal incident energy
!
          do
            flagdel = .false.
            do i = 1, N - 2
              if (E10(MT, iso, i) >= E10(MT, iso, i + 1) .and. E10(MT, iso, i) == E10(MT, iso, i + 2)) then
                do j = i + 1, N - 1
                  E10(MT, iso, j) = E10(MT, iso, j + 1)
                  xsiso(MT, iso, j) = xsiso(MT, iso, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iso=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iso, E10(MT, iso, i + 1), xsiso(MT, iso, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Delete multiple ( > 2) points with equal cross section or double points with equal energies and cross sections
!
          do
            flagdel = .false.
            do i = 1, N - 2
              if ((xsiso(MT, iso, i) == xsiso(MT, iso, i + 1) .and. xsiso(MT, iso, i) == xsiso(MT, iso, i + 2)) .or. &
 &              (E10(MT, iso, i) == E10(MT, iso, i + 1) .and. xsiso(MT, iso, i) == xsiso(MT, iso, i + 1))) then
                do j = i + 1, N - 1
                  E10(MT, iso, j) = E10(MT, iso, j + 1)
                  xsiso(MT, iso, j) = xsiso(MT, iso, j + 1)
                enddo
                flagdel = .true.
                write(8, '(" Deleted point: MF=", i3, " MT=", i3, " iso=", i3, " E=", es12.5, " eV xs=", es12.5, " b")') &
 &                MF, MT, iso, E10(MT, iso, i + 1), xsiso(MT, iso, i + 1)
                iclean = iclean + 1
                exit
              endif
            enddo
            if (.not. flagdel) exit
            N = N - 1
          enddo
!
! Set new number of points
!
          NP10(MT, iso) = N
          NBT10(MT, iso, 1) = N
        enddo
      endif
    enddo
  enddo
  return
end subroutine make9_10clean
! Copyright A.J. Koning 2021
