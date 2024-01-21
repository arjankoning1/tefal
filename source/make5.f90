subroutine make5
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF5
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
!   nummt       ! number of MT numbers
! Variables for input of specific ENDF data
!   adopt       ! logical for existence of MF information (per MT)
! Variables for info from TALYS
!   Atarget     ! mass number of nucleus
!   Ztarget     ! charge number of nucleus
! Variables for initialization of ENDF format
!   MAT         ! MAT number
!   mfexist     ! flag for existence of MF - number
!   mtexist     ! flag for existence of MT - number
!   blank2      ! blank string
!   FEND        ! ENDF - 6 format
! Variables for ENDF format
!   NK          ! number of subsections
! Variables for MF1
!   EMAX        ! upper limit of energy range for evaluation
! Variables for MF3
!   E3          ! incident energy for MF3 (in ENDF - 6 format)
! Variables for MF5
!   E5          ! incident energy for MF5 (in ENDF - 6 format)
!   E5p         ! incident energy for which tabulated distribution is given
!   EFH         ! constant used in the energy - dependent fission neutron spectrum
!   EFL         ! constant used in the energy - dependent fission neutron spectrum
!   gE5         ! energy - spectrum values
!   INTER5      ! interpolation scheme
!   INTER5e     ! interpolation scheme
!   INTER5e2    ! interpolation scheme
!   LF          ! flag for energy distribution law
!   NBT5        ! separation value for interpolation scheme
!   NBT5e       ! separation value for interpolation scheme
!   NBT5e2      ! separation value for interpolation scheme
!   NE5e        ! number of incident energies for distribution
!   NF          ! number of secondary energy points
!   NP5         ! number of incident energies
!   NR5         ! number of interpolation ranges
!   NR5e        ! number of interpolation ranges
!   NR5e2       ! number of interpolation ranges
!   pE          ! fractional part of cross section
!   TM5         ! Effective (7) or maximum (12) temperature FNS parameter, T(E)
!   U           ! constant for upper energy limit
!
! *** Declaration of local data
!
  implicit none
  integer :: i                 ! counter
  integer :: k                 ! counter
  integer :: MF                ! MF-number
  integer :: MT                ! MT-number
  integer :: Nnew              ! help variable
!
! ***************************** Make MF5 *******************************
!
! read5   : subroutine to read MF5 from existing ENDF-6 data library
! write5  : subroutine to write MF5
!
  MF = 5
  open (unit = 2, file = 'MF5', status = 'replace')
  do MT = 1, nummt
    if ( .not. mtexist(3, MT) .and. .not. (MT == 455 .and. mtexist(1, MT))) cycle
    NK(5, MT) = 0
    if (adopt(MF, MT)) call read5(MT)
!
! Test and adjust energy boundaries
!
    do k = 1, NK(5, MT)
!
! Incident photons: shift energies by Q-value
!
      if (k0 == 0) then
        do i = 1, NE5e(k)
          E5(k, i) = E5(k, i) - QM(4)
        enddo
        do i = 1, NP5(k)
          E5p(k, i) = E5p(k, i) - QM(4)
        enddo
        E5(k, 1) = E3(18, 1)
        E5p(k, 1) = E3(18, 1)
      endif
!
! A. Remove highest energies
!
      if (E5p(k, NP5(k)) > EMAX) then
        E5p(k, NP5(k)) = EMAX
        Nnew = NE5e(k)
        do i = 1, NE5e(k)
          if (E5(k, i) > EMAX) Nnew = Nnew - 1
        enddo
        NE5e(k) = Nnew
        NBT5e(k, NR5e(k)) = Nnew
!
! B. Add high energies
!
      else
        if (E5p(k, NP5(k)) < EMAX) then
          E5p(k, NP5(k)) = EMAX
          if (E5(k, NE5e(k)) < EMAX) then
            NE5e(k) = NE5e(k) + 1
            NBT5e(k, NR5e(k)) = NBT5e(k, NR5e(k)) + 1
            E5(k, NE5e(k)) = EMAX
            NR5e2(k, NE5e(k)) = NR5e2(k, NE5e(k) - 1)
            do i = 1, NR5e2(k, NE5e(k))
              NBT5e2(k, NE5e(k), i) = NBT5e2(k, NE5e(k) - 1, i)
              INTER5e2(k, NE5e(k), i) = INTER5e2(k, NE5e(k) - 1, i)
            enddo
            NF(k, NE5e(k)) = NF(k, NE5e(k) - 1)
            do i = 1, 2 * NF(k, NE5e(k))
              gE5(k, NE5e(k), i) = gE5(k, NE5e(k) - 1, i)
            enddo
            TM5(k, 2 * NE5e(k) - 1) = 1.00001 * TM5(k, 2 * NE5e(k) - 3)
            TM5(k, 2 * NE5e(k)) = 0.
          endif
        endif
      endif
    enddo
!
! Default Maxwellian spectrum if no better alternative exists
!
    U(1) = -2.e7
    if (MT == 18 .and. NK(5, MT) == 0) then
      NK(5, MT) = 1
      LF(1) = 7
      NR5(1) = 1
      NP5(1) = 2
      NBT5(1, 1) = 2
      INTER5(1, 1) = 2
      E5p(1, 1) = E3(18, 1)
      E5p(1, 2) = EMAX
      pE(1, 1) = 1.
      pE(1, 2) = 1.
      NR5e(1) = 1
      NE5e(1) = 2
      EFL = 0.
      EFH = 0.
      NBT5e(1, 1) = 2
      INTER5e(1, 1) = 2
      TM5(1, 1) = E3(18, 1)
      TM5(1, 2) = 37000. * (Ztarget **2) / (real(Atarget))
      TM5(1, 3) = EMAX
      TM5(1, 4) = TM5(1, 2)
      mtexist(5, MT) = .true.
      mfexist(5) = .true.
    endif
!
! Write to MF5
!
    if (mtexist(5, MT)) call write5(MT)
  enddo
  write(2, fmt = FEND) blank2, MAT, 0, 0, 0
  close (unit = 2)
  return
end subroutine make5
! Copyright A.J. Koning 2021
