subroutine write12(MT)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write MF12
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
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numenin        ! number of incident energies
!   numint         ! number of interpolation sections
!  Variables for partial cross sections in ENDF format
!   idchannel      ! identifier for channel
! Variables for initialization of ENDF format
!   AWR            ! standard mass parameter
!   blank2         ! blank string
!   idnum          ! number of different exclusive cross sections
!   MAT            ! MAT number
!   MTid           ! channel identifier for MT - number
!   SEND           ! ENDF - 6 format
!   ZA             ! standard charge parameter
! Variables for ENDF format
!   INTER          ! interpolation scheme
!   NBT            ! separation value for interpolation scheme
!   NK          ! number of subsections
!   NP             ! number of incident energies
!   NR             ! number of interpolation ranges
! Variables for MF12_15
!   E12            ! incident energy (in ENDF - 6 format)
!   Eg             ! gamma energy
!   Egk            ! gamma energy
!   ES12           ! energy of level
!   Esk            ! starting level (in ENDF - 6 format)
!   ESNS           ! energy of mother level
!   INTERg         ! interpolation scheme
!   LFg            ! photo energy distribution law
!   LG12           ! type setters
!   LO12           ! type setters
!   LP12           ! origin of photons
!   LPg            ! primary photon flag
!   NBTg           ! separation value for interpolation scheme
!   NPg            ! number of incident energies
!   NS12           ! number of levels below the present one
!   NT12           ! number of transitions for which data are given
!   TP12           ! probability of direct transition
!   xsgtotyield    ! total discrete photon multiplicity
!   xsgyield       ! gamma - ray multiplicity (in ENDF - 6 format)
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: id                        ! counter for deuterons
  integer   :: idc                       ! help variable
  integer   :: igam                      ! counter for gammas
  integer   :: ii                        ! counter
  integer   :: MT                        ! MT-number
  integer   :: MF                        ! MF-number
  integer   :: N                         ! neutron number of residual nucleus
  integer   :: NS                        ! line number
  integer   :: Nval(2*numint)            ! value
  real(sgl) :: x(2*numenin)              ! help variable
!
! ***************************** Write MF12 *****************************
!
! hrwrite  : subroutine to write header with real values
!
  MF = 12
  NS = 0
!
! LO=1: Multiplicities
!
  if (LO12(MT) == 1) then
    call hrwrite(ZA, AWR, LO12(MT), LG12(MT), NK(MF, MT), 0, MAT, MF, MT, NS)
    if (NK(MF, MT) > 1) then
      call hrwrite(0., 0., 0, 0, NR(MF, MT), NP(MF, MT), MAT, MF, MT, NS)
!
! 1. Total discrete gamma-ray multiplicities
!
! Write interpolation ranges
!
! kwrite     : subroutine to write integer value block
!
      do i = 1, NR(MF, MT)
        ii = 2 * i - 1
        Nval(ii) = NBT(MF, MT, i)
        Nval(ii + 1) = INTER(MF, MT, i)
      enddo
      N = 2 * NR(MF, MT)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write multiplicities
!
! xwrite     : subroutine to write real value block
!
      do id = 0, idnum
        if (idchannel(id) == MTid(MT)) then
          idc = id
          exit
        endif
      enddo
      do i = 1, NP(MF, MT)
        ii = 2 * i - 1
        x(ii) = E12(idc, i)
        x(ii + 1) = xsgtotyield(idc, i)
      enddo
      N = 2 * NP(MF, MT)
      call xwrite(N, x, MAT, MF, MT, NS)
    endif
!
! 2. Gamma-ray multiplicities per discrete transition
!
    do igam = 1, NK(MF, MT)
      call hrwrite(Egk(idc, igam), Esk(idc, igam), LPg(MT, igam), LFg(MT, igam), NR(MF, MT), NPg(MT, igam), MAT, MF, MT, NS)
!
! Write interpolation ranges
!
      do i = 1, NR(MF, MT)
        ii = 2 * i - 1
        Nval(ii) = NBTg(MT, igam, i)
        Nval(ii + 1) = INTERg(MT, igam, i)
      enddo
      N = 2 * NR(MF, MT)
      call kwrite(N, Nval, MAT, MF, MT, NS)
!
! Write cross sections
!
      do i = 1, NPg(MT, igam)
        ii = 2 * i - 1
        x(ii) = Eg(idc, igam, i)
        x(ii + 1) = xsgyield(idc, igam, i)
      enddo
      N = 2 * NPg(MT, igam)
      call xwrite(N, x, MAT, MF, MT, NS)
    enddo
  else
!
! LO=2: Transition probability arrays
!
    call hrwrite(ZA, AWR, LO12(MT), LG12(MT), NS12(MT), 0, MAT, MF, MT, NS)
    call hrwrite(ESNS(MT), 0., LP12(MT), 0, (LG12(MT) + 1) * NT12(MT), NT12(MT), MAT, MF, MT, NS)
    do i = 1, NT12(MT)
      ii = 2 * i - 1
      x(ii) = ES12(MT, i)
      x(ii + 1) = TP12(MT, i)
    enddo
    N = 2 * NT12(MT)
    call xwrite(N, x, MAT, MF, MT, NS)
  endif
  write(2, fmt = SEND) blank2, MAT, MF, 0, NS
  return
end subroutine write12
! Copyright A.J. Koning 2021
