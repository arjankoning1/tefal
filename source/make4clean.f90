subroutine make4clean(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Clean up MF6
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
! Variables for initialization of ENDF variables
!   iclean    ! number of cleaned points
! Variables for ENDF format
!   NBT       ! separation value for interpolation scheme
! Variables for MF4
!   E4        ! incident energy for MF4 (in ENDF - 6 format)
!   f4        ! angular distribution
!   leg       ! Legendre coefficients (in ENDF - 6 format)
!   NBT4      ! separation value for interpolation scheme
!   NE        ! number of incident energies
!   NEh       ! number of incident energies (MF4 only)
!   NL        ! Legendre order or number of cosines
!   NP4       ! number of incident energies
!   x4        ! cosine of the angle
!
! *** Declaration of local data
!
  implicit none
  logical :: flagdel              ! flag to delete a line
  integer :: i                    ! counter
  integer :: iE                   ! energy counter
  integer :: j                    ! counter
  integer :: L                    ! counter for Legendre coefficients
  integer :: MF                   ! MF-number
  integer :: MT                   ! MT-number
  integer :: N                    ! neutron number of residual nucleus
!
! To make sure we get a "clean" MF4 we do some final checking and take actions if necessary.
! Of course, we report such actions in the output file
!
! ****************** Remove redundant points ***************************
!
  MF = 4
  if (MT == 18) return
  N = NE
!
! Delete multiple ( > 2) points with trivial Legendre coefficients
!
  do
    flagdel = .false.
Loop1:  do i = 1, N - 2
      if (NL(MF, MT, i) == NL(MF, MT, i + 1) .and. NL(MF, MT, i) == NL(MF, MT, i + 2)) then
        do L = 1, NL(MF, MT, i)
          if (leg(i, L) /= leg(i + 1, L)) cycle Loop1
          if (leg(i, L) /= leg(i + 2, L)) cycle Loop1
        enddo
        do j = i + 1, N - 1
          E4(j) = E4(j + 1)
          NL(MF, MT, j) = NL(MF, MT, j + 1)
          do L = 1, NL(MF, MT, j)
            leg(j, L) = leg(j + 1, L)
          enddo
        enddo
        flagdel = .true.
        write(8, '(" Deleted Legendre coefficients: MF=", i3, " MT=", i3, " E=", es12.5, " eV", es12.5)') MF, MT, E4(i+1)
        iclean = iclean + 1
        exit
      endif
    enddo Loop1
    if (.not. flagdel) exit
    N = N - 1
  enddo
!
! Set new number of points
!
  NE = N
  NBT(MF, MT, 1) = N
!
! Delete multiple ( > 2) points with trivial angular distributions
!
  do iE = 1, NEh
    N = NP4(iE)
    do
      flagdel = .false.
      do i = 1, N - 2
        if (f4(iE, i) == f4(iE, i + 1) .and. f4(iE, i) == f4(iE, i + 2)) then
          do j = i + 1, N - 1
            x4(iE, j) = x4(iE, j + 1)
            f4(iE, j) = f4(iE, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted angular distribution points: MF=", i3," MT=", i3, " E=", es12.5, " eV", es12.5)') &
 &          MF, MT, x4(iE, i+1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if (.not.flagdel) exit
      N = N - 1
    enddo
!
! Set new number of points
!
    NP4(iE) = N
    NBT4(iE, 1) = N
  enddo
  return
end subroutine make4clean
! Copyright A.J. Koning 2021
