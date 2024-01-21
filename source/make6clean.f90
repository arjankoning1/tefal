subroutine make6clean(MT)
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
! Variables for input of ENDF structure
!   flagrecoil    ! flag to include recoil information
! Variables for initialization of ENDF format
!   iclean        ! number of cleaned points
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF6
!   b6            ! energy - angle values
!   b6gam         ! energy - angle values for photons
!   b6rec         ! energy - angle values for recoils
!   E6            ! incident energy (in ENDF - 6 format) for distribution
!   Ey            ! incident energy for yields (in ENDF - 6 format)
!   flagrec       ! flag to state that b6 element concerns recoil grid
!   kpart         ! section number for particles
!   LAW           ! flag for distribution function
!   NA            ! number of angular parameters
!   NBT6ea        ! separation value for interpolation scheme
!   NBT6y         ! separation value for interpolation scheme
!   ND            ! number of discrete energies
!   NE6ea         ! number of incident energies for distribution
!   NEP           ! number of secondary energy points
!   NP6y          ! number of incident energies for yields
!   NW            ! number of words
!   Y             ! product yield (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  integer :: dtype                ! data type
  logical :: flagdel              ! flag to delete a line
  integer :: i                    ! counter
  integer :: j                    ! counter
  integer :: k                    ! counter
  integer :: kk                   ! counter
  integer :: l                    ! counter
  integer :: MF                   ! MF-number
  integer :: MT                   ! MT-number
  integer :: N                    ! neutron number of residual nucleus
  integer :: Ncont                ! number of continuum energies
  integer :: Ncont1               ! number of continuum energies
!
! Many actions may be required to construct MF6.
! To make sure we get a "clean" MF6 we do some final checking and
! take actions if necessary. Of course, we report such actions in
! the output file
!
! ****************** Remove redundant points ***************************
!
  MF = 6
  do k = 1, NK(MF, MT)
    dtype = 1
    if (MT == 5 .or. MT == 18) then
      if (k > kpart .and. k < NK(MF, MT) .and. flagrec(k)) dtype = 2
      if (k == NK(MF, MT)) dtype = 3
    endif
    N = NP6y(k)
!
! Delete decreasing energies
!
    do
      flagdel = .false.
      do i = 2, N
        if (Ey(k, i) < Ey(k, i - 1)) then
          do j = i, N - 1
            Ey(k, j) = Ey(k, j + 1)
            Y(k, j) = Y(k, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted decreasing energy: MF=", i3, " MT=", i3, " k=", i3, &
 &          " E=", es12.5, " eV  Y=", es12.5)') MF, MT, k, Ey(k, i+1), Y(k, i + 1)
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
        if (Ey(k, i) >= Ey(k, i + 1) .and. Ey(k, i) == Ey(k, i + 2)) then
          do j = i + 1, N - 1
            Ey(k, j) = Ey(k, j + 1)
            Y(k, j) = Y(k, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted equal incident energy: MF=", i3, " MT=", i3, " k=", i3, &
 &          " E=", es12.5, " eV  Y=", es12.5)') MF, MT, k, Ey(k, i+1), Y(k, i + 1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if (.not. flagdel) exit
      N = N - 1
    enddo
!
! Delete multiple ( > 2) points with equal yield
!
    do
      flagdel = .false.
      do i = 1, N - 2
        if (Y(k, i) == Y(k, i + 1) .and. Y(k, i) == Y(k, i + 2)) then
          do j = i + 1, N - 1
            Ey(k, j) = Ey(k, j + 1)
            Y(k, j) = Y(k, j + 1)
          enddo
          flagdel = .true.
          write(8, '(" Deleted equal yield: MF=", i3, " MT=", i3, " k=", i3, &
 &          " E=", es12.5, " eV  Y=", es12.5)') MF, MT, k, Ey(k, i+1), Y(k, i + 1)
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
    NP6y(k) = N
    NBT6y(k, 1) = N
!
! Check for, and remove, negative yields
!
    do i = 1, N
      if (Y(k, i) < 0.) then
        if (i > 1) then
          Y(k, i) = Y(k, i - 1)
        else
          Y(k, i) = 0.
        endif
      endif
    enddo
!
! Delete multiple ( > 2) sets with trivial energy-angle distributions
!
    if (k > kpart .and. k < NK(MF, MT) .and. .not. flagrecoil) cycle
    if (LAW(k) > 1) cycle
    N = NE6ea(k)
    do
      flagdel = .false.
Loop1: do i = 1, N - 2
        if (NW(k, i) == NW(k, i + 1) .and. NW(k, i) == NW(k, i + 2)) then
          if (dtype == 1) then
            do kk = 1, NW(k, i)
              if (b6(k, i, kk) /= b6(k, i + 1, kk)) cycle Loop1
              if (b6(k, i, kk) /= b6(k, i + 2, kk)) cycle Loop1
            enddo
          endif
          if (dtype == 2) then
            do kk = 1, NW(k, i)
              if (b6rec(k, i, kk) /= b6rec(k, i + 1, kk)) cycle Loop1
              if (b6rec(k, i, kk) /= b6rec(k, i + 2, kk)) cycle Loop1
            enddo
          endif
          if (dtype == 3) then
            do kk = 1, NW(k, i)
              if (b6gam(i, kk) /= b6gam(i + 1, kk)) cycle Loop1
              if (b6gam(i, kk) /= b6gam(i + 2, kk)) cycle Loop1
            enddo
          endif
          do j = i + 1, N - 1
            E6(k, j) = E6(k, j + 1)
            NEP(k, j) = NEP(k, j + 1)
            NA(k, j) = NA(k, j + 1)
            ND(k, j) = ND(k, j + 1)
            NW(k, j) = NW(k, j + 1)
            if (dtype == 1) then
              do kk = 1, NW(k, j)
                b6(k, j, kk) = b6(k, j + 1, kk)
              enddo
            endif
            if (dtype == 2) then
              do kk = 1, NW(k, j)
                b6rec(k, j, kk) = b6rec(k, j + 1, kk)
              enddo
            endif
            if (dtype == 3) then
              do kk = 1, NW(k, j)
                b6gam(j, kk) = b6gam(j + 1, kk)
              enddo
            endif
          enddo
          flagdel = .true.
          write(8, '(" Deleted energy distribution: MF=", i3, " MT=", &
 &          i3, " k=", i3, " E=", es12.5, " eV")') MF, MT, k, E6(k, i+1)
          iclean = iclean + 1
          exit
        endif
      enddo Loop1
      if (.not. flagdel) exit
      N = N - 1
    enddo
    do
      flagdel = .false.
      do i = 1, N - 1
        if (E6(k, i) == E6(k, i + 1)) then
          do j = i, N - 1
            E6(k, j) = E6(k, j + 1)
            NEP(k, j) = NEP(k, j + 1)
            NA(k, j) = NA(k, j + 1)
            ND(k, j) = ND(k, j + 1)
            NW(k, j) = NW(k, j + 1)
            if (dtype == 1) then
              do kk = 1, NW(k, j)
                b6(k, j, kk) = b6(k, j + 1, kk)
              enddo
            endif
            if (dtype == 2) then
              do kk = 1, NW(k, j)
                b6rec(k, j, kk) = b6rec(k, j + 1, kk)
              enddo
            endif
            if (dtype == 3) then
              do kk = 1, NW(k, j)
                b6gam(j, kk) = b6gam(j + 1, kk)
              enddo
            endif
          enddo
          flagdel = .true.
          write(8, '(" Deleted energy distribution: MF=", i3, " MT=", &
 &          i3, " k=", i3, " E=", es12.5, " eV")') MF, MT, k, E6(k, i+1)
          iclean = iclean + 1
          exit
        endif
      enddo
      if (.not. flagdel) exit
      N = N - 1
    enddo
    NE6ea(k) = N
    NBT6ea(k, 1) = N
!
! Use reasonable values at threshold energy
!
    if (LAW(k) == 1) then
      do i = 1, N - 1
        Ncont = NEP(k, i) - ND(k, i)
        if (Ncont == 2) then
          do l = i + 1, N - 1
            Ncont1 = NEP(k, l) - ND(k, l)
            if (Ncont1 > 2) then
              kk = l
              j = 2 * ND(k, i) + NA(k, i) + 3
              if (dtype == 1) then
                if (b6(k, kk, j) > 0.) then
                  b6(k, i, j) = b6(k, kk, j)
                  if (b6(k, i, 2 * ND(k, i) + 2) == 1.) b6(k, i, 2 * ND(k, i) + 2) = 1. / b6(k, i, j)
                endif
              endif
              if (dtype == 2) then
                if (b6rec(k, kk, j) > 0.) then
                  b6rec(k, i, j) = b6rec(k, kk, j)
                  b6rec(k, i, 2 * ND(k, i) + 2) = 1. / b6rec(k, i, j)
                endif
              endif
              if (dtype == 3) then
                if (b6gam(kk, j) > 0.) then
                  b6gam(i, j) = b6gam(kk, j)
                  b6gam(i, 2 * ND(k, i) + 2) = 1. / b6gam(i, j)
                endif
              endif
              if (i > 1) then
                if (NW(k, i) == NW(k, i - 1) .and. NEP(k, i) == NEP(k, i - 1)) then
                  do j = 1, NW(k, i)
                    if (dtype == 1) b6(k, i, j) = b6(k, i - 1, j)
                    if (dtype == 2) b6rec(k, i, j) = b6rec(k, i - 1, j)
                    if (dtype == 3) b6gam(i, j) = b6gam(i - 1, j)
                  enddo
                endif
              endif
              exit
            endif
          enddo
        endif
      enddo
!
! Set first MT5 energy to trivial value
!
      if (MT == 5 .and. dtype == 1 .and. E6(k, 1) == EmineV .and. NW(k, 1) == 6) then
        b6(k, 1, 4) = 1.e-6
        b6(k, 1, 2) = 1.e+6
      endif
    endif
  enddo
  do k = 1, NK(MF, MT)
    if (Ey(k, 2) == Ey(k, 3) .and. Y(k, 2) == 0.) then
      Y(k, 2) = Y(k, 3)
      Y(k, 1) = Y(k, 3)
    endif
  enddo
  return
end subroutine make6clean
! Copyright A.J. Koning 2021
