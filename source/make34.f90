subroutine make34
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make MF34
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
!   numenin     ! number of incident energies
! Variables for initialization of ENDF format
!   mfexist     ! flag for existence of MF - number
!   mtexist     ! flag for existence of MT - number
! Variables for covariances in ENDF format
!   Eleg        ! energy for Legendre covariance data
!   Nchanleg    ! number of energies with Legendre covariance dat
!   Nleg34      ! number of Legendre coefficients with covariances
!   Rleg        ! covariance matrix element for Legendre coeff.
! Variables for MF31_40
!   b34         ! covariance matrix element
!   LB34        ! flag for meaning of numbers
!   LS34        ! symmetry flag
!   LTT34       ! representation
!   LVT34       ! specification of transformation matrix
!   MAT34       ! material number for MF34
!   MT34        ! MT number for MF34
!   NE34        ! number of energies in energy array
!   NI34        ! number of NI - type sub - subsections
!   NL34        ! number of Legendre coefficients with covariances
!   NL341       ! number of Legendre coefficients with covariances
!   NMT34       ! total number of MT sections
!   NT34        ! total number of entries
!
! *** Declaration of local data
!
  implicit none
  integer :: covstep                    ! energy step for covariances
  integer :: i                          ! counter
  integer :: icov(numenin)              ! index for covariance
  integer :: iE                         ! energy counter
  integer :: iE0                        ! counter of energies
  integer :: iE1                        ! counter of energies
  integer :: j                          ! counter
  integer :: L1                         ! integration limits
  integer :: L2                         ! integration limits
  integer :: MT                         ! MT-number
!
! ***************************** Make MF34 ******************************
!
  MT = 2
  covstep = 2
  do L1 = 1, Nleg34
    do L2 = L1, Nleg34
      iE = 0
      do i = 1, Nchanleg
        if (mod(i, covstep) == 1 .and. i /= 1 .and. i /= Nchanleg) cycle
        iE = iE + 1
        icov(iE) = i
        b34(L1, L2, iE) = Eleg(i) * 1.e6
      enddo
      NE34(L1, L2) = iE
      do iE0 = 1, NE34(L1, L2) - 1
        i = icov(iE0)
        do iE1 = iE0, NE34(L1, L2) - 1
          j = icov(iE1)
          iE = iE + 1
          b34(L1, L2, iE) = Rleg(i, L1, j, L2)
        enddo
      enddo
      NI34(L1, L2) = 1
      NT34(L1, L2) = (NE34(L1, L2) * (NE34(L1, L2) + 1)) / 2
    enddo
  enddo
!
! 2. ENDF-6 parameters
!
  LVT34 = 0
  LTT34 = 1
  NMT34 = 1
  MAT34 = 0
  MT34 = 2
  NL34 = Nleg34
  NL341 = Nleg34
  LS34 = 1
  LB34 = 5
  mtexist(34, MT) = .true.
  mfexist(34) = .true.
  return
end subroutine make34
! Copyright A.J. Koning 2021
