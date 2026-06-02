!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module ascii  --  selected values of ASCII character codes  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     tilde        !<decimal value of ASCII code for tilde (126)
!
!
module ascii
   implicit none
   integer ::  null       !<decimal value of ascii code for null (0)
   integer :: tab         !<decimal value of ASCII code for tab (9)
   integer :: linefeed    !<decimal value of ASCII code for linefeed (10)     
   integer :: formfeed    !<decimal value of ASCII code for formfeed (12)
   integer :: carriage!<decimal value of ASCII code for carriage return (13)
   integer :: escape!<decimal value of ASCII code for escape (27)
   integer :: space!<decimal value of ASCII code for blank space (32)
   integer :: exclamation!<decimal value of ASCII code for exclamation (33)
   integer :: quote!<decimal value of ASCII code for double quote (34)
   integer :: pound!<decimal value of ASCII code for pound sign (35)
   integer :: dollar!<decimal value of ASCII code for dollar sign (36)
   integer :: percent!<decimal value of ASCII code for percent sign (37)
   integer :: ampersand!<decimal value of ASCII code for ampersand (38)
   integer :: apostrophe!<decimal value of ASCII code for single quote (39)
   integer :: asterisk!<decimal value of ASCII code for asterisk (42)
   integer :: plus!<decimal value of ASCII code for plus sign (43)
   integer :: comma!<decimal value of ASCII code for comma (44)
   integer :: minus!<decimal value of ASCII code for minus sign (45)
   integer :: period!<decimal value of ASCII code for period (46)
   integer :: frontslash!<decimal value of ASCII codd for frontslash (47)
   integer :: colon!<decimal value of ASCII code for colon (58)
   integer :: semicolon!<decimal value of ASCII code for semicolon (59)
   integer :: equal!<decimal value of ASCII code for equal sign (61)
   integer :: question!<decimal value of ASCII code for question mark (63)
   integer :: atsign!<decimal value of ASCII code for at sign (64)
   integer :: backslash!<decimal value of ASCII code for backslash (92)
   integer :: caret!<decimal value of ASCII code for caret (94)
   integer :: underbar!<decimal value of ASCII code for underbar (95)
   integer :: vertical!<decimal value of ASCII code for vertical bar (124)
   integer :: tilde!<decimal value of ASCII code for tilde (126)
   parameter (null=0)
   parameter (tab=9)
   parameter (linefeed=10)
   parameter (formfeed=12)
   parameter (carriage=13)
   parameter (escape=27)
   parameter (space=32)
   parameter (exclamation=33)
   parameter (quote=34)
   parameter (pound=35)
   parameter (dollar=36)
   parameter (percent=37)
   parameter (ampersand=38)
   parameter (apostrophe=39)
   parameter (asterisk=42)
   parameter (plus=43)
   parameter (comma=44)
   parameter (minus=45)
   parameter (period=46)
   parameter (frontslash=47)
   parameter (colon=58)
   parameter (semicolon=59)
   parameter (equal=61)
   parameter (question=63)
   parameter (atsign=64)
   parameter (backslash=92)
   parameter (caret=94)
   parameter (underbar=95)
   parameter (vertical=124)
   parameter (tilde=126)
   save

contains

!> @brief 
!> converts integer to strings
!> @param[in] value integer: integer value to be converted
!> @param[in] n_char_min integer: minimum number of char
!> @param[out] str char(:): character array resulting from conversion
   function int_to_str(value,n_char_min) result(str)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: value
      INTEGER, INTENT(in), OPTIONAL :: n_char_min
      CHARACTER(:), ALLOCATABLE :: str
      INTEGER :: n_char
      CHARACTER(10) :: n_char_char

      if(value==0) then
         n_char=1
      else
         n_char=int(log10(real(value)))+1
      endif
      if (present(n_char_min)) then
         if(n_char<n_char_min) n_char=n_char_min
      endif
      !write(0,*) n_char
      allocate(character(n_char) :: str)
      write(n_char_char,'(i10.10)') n_char
      ! write(0,*) n_char_char
      write(str,'(i'//trim(n_char_char)//'.'//trim(n_char_char)//')')&
      &value
   end function int_to_str

end
