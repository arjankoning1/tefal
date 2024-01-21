subroutine tefalconstants
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Constants and basic properties of particles
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
! Constants
!   fislim       ! mass above which nuclide fissions
!   light        ! mass number of lightest stable isotope
!   nuc          ! nuclide
!   parA         ! mass number of particle
!   parN         ! neutron number of particle
!   parmass      ! mass of particle in a.m.u.
!   parname      ! name of particle
!   parspin      ! spin of particle
!   parsym       ! symbol of particle
!   parZ         ! charge number of particle
!   pi           ! pi
!   xsepshigh    ! upper limit for cross sections in millibarns
!   xsepslow     ! lower limit for cross sections in millibarns
!
! *** Declaration of local data
!
  implicit none
!
! *********************** Fundamental constants ************************
!
  pi = 3.14159265358979323
!
! ****************** General properties of particles *******************
!
!   fission = -1
!   photon  = 0
!   neutron = 1
!   proton  = 2
!   deuteron= 3
!   triton  = 4
!   helium-3= 5
!   alpha   = 6
!
  parname = (/'fission ', 'photon  ', 'neutron ', 'proton  ', 'deuteron', 'triton  ', 'helium-3', 'alpha   '/)
  parsym = (/'f', 'g', 'n', 'p', 'd', 't', 'h', 'a'/)
  parZ = (/ 0, 0, 1, 1, 1, 2, 2 /)
  parN = (/ 0, 1, 0, 1, 2, 1, 2 /)
  parA = (/ 0, 1, 1, 2, 3, 3, 4 /)
  parmass =  (/ 0., 1.008664904, 1.007276470, 2.013553214, 3.016049268, 3.016029310, 4.002603250 /)
  parspin = (/ 0., 0.5, 0.5, 1., 0.5, 0.5, 0. /)
!
! ************************ Nuclear symbols *****************************
!
  nuc = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'B9', 'C0', 'C1', 'C2', 'C3', 'C4'/)
!
! Numbers for Po, At, Rn, Fr and nuclides heavier than Es are unknown.
!
  light = (/ 1,  3,  6,  9, 10, 12, 14, 16, 19, 20, &
     23, 24, 27, 28, 31, 32, 35, 36, 39, 40, 45, 46, 50, 50, 55, 54, 59, 58, 63, 64, &
     69, 70, 75, 74, 79, 78, 85, 84, 89, 90, 93, 92, 99, 96, 103, 102, 107, 106, 113, 112, &
    121, 120, 127, 124, 133, 130, 138, 136, 141, 142, 139, 144, 151, 152, 159, 156, 165, 162, 169, 168, &
    175, 174, 180, 180, 185, 184, 191, 190, 197, 196, 203, 204, 209, 206, 203, 211, 212, 223, 225, 227, &
    229, 234, 230, 235, 235, 240, 240, 240, 241, 240, 245, 248, 252, 255, 258, 261, 264, 267, 270, 273, &
    276, 279, 282, 285, 288, 291, 294, 297, 300, 303, 306, 309, 312, 315 /)
!
! *********************** Limit for cross sections *********************
!
! If the maximal cross section in an excitation cross section does not exceed xsepshigh, the MT-number is omitted altogether.
! Cross sections lower than xsepslow are assumed to have no physical meaning and are set to zero.
!
  xsepslow = 1.000001e-17
  xsepshigh = 1.e-6
  fislim = 215
end subroutine tefalconstants
! Copyright A.J. Koning 2021
