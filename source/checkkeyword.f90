subroutine checkkeyword
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in keywords
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
!   inline         ! input line
!   nlines         ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  integer, parameter  :: numkey=59         ! number of keywords
  integer             :: i                 ! counter
  integer             :: j                 ! counter
  character(len=132)  :: key               ! keyword
  character(len=132)  :: keyword(numkey)   ! keyword
  character(len=132)  :: word(40)          ! words on input line
!
! Although it is difficult to prevent the user from all possible input errors, we can check for the use of wrong keywords
! and for unphysical values for most of the input variables.
!
! *********************** Check for wrong keywords *********************
!
! TEFAL will stop if a keyword is incorrect
!
  data (keyword(i), i = 1, numkey) / ' ', 'addlow', 'adopt', 'author', 'background', 'breakup', &
    'capt6', 'clean', 'covariance', 'covdiscrete', 'covrp', 'cross', 'cuteps', 'diffweight', 'disc6', 'disclim', 'eaf', &
    'endfdetail', 'endffile', 'endftext', 'eswitch', 'eswitch4', 'exclude', 'fis10', 'gam13', 'gamdiscrete', &
    'gpf', 'high', 'identifier', 'include', 'intercor', 'lab', 'lssf', 'maxrp', 'mtall', 'mtextra', 'multichance', &
    'ngn', 'nmtmax', 'nocross', 'nomf', 'part6', 'partiala', 'partialcov', 'partiald', 'partialh', 'partialn', 'partialp', &
    'partialt', 'recoil', 'renorm', 'resonance', 'rp10', 'rp6', 'subfission', 'tabddx', 'urrcomp', 'urrenergy', 'urrmode'/
!
! A keyword can be de-activated by putting a # in front of it.
! All first words of the input lines are checked against the list of keywords.
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified.
!
Loop1:  do i = 1, nlines
    call getkeywords(inline(i), word)
    key = word(1)
    if (key(1:1) == '#') cycle
    do j = 1, numkey
      if (keyword(j) == key) cycle Loop1
    enddo
    write(*, '(/" TEFAL-error: Wrong keyword: ", a20)') key
    stop
  enddo Loop1
  return
end subroutine checkkeyword
! Copyright A.J. Koning 2021
