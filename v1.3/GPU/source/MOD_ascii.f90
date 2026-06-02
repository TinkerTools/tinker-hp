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
!     null         decimal value of ASCII code for null (0)
!     tab          decimal value of ASCII code for tab (9)
!     linefeed     decimal value of ASCII code for linefeed (10)
!     formfeed     decimal value of ASCII code for formfeed (12)
!     carriage     decimal value of ASCII code for carriage return (13)
!     escape       decimal value of ASCII code for escape (27)
!     space        decimal value of ASCII code for blank space (32)
!     exclamation  decimal value of ASCII code for exclamation (33)
!     quote        decimal value of ASCII code for double quote (34)
!     pound        decimal value of ASCII code for pound sign (35)
!     dollar       decimal value of ASCII code for dollar sign (36)
!     percent      decimal value of ASCII code for percent sign (37)
!     ampersand    decimal value of ASCII code for ampersand (38)
!     apostrophe   decimal value of ASCII code for single quote (39)
!     asterisk     decimal value of ASCII code for asterisk (42)
!     plus         decimal value of ASCII code for plus sign (43)
!     comma        decimal value of ASCII code for comma (44)
!     minus        decimal value of ASCII code for minus sign (45)
!     period       decimal value of ASCII code for period (46)
!     frontslash   decimal value of ASCII codd for frontslash (47)
!     colon        decimal value of ASCII code for colon (58)
!     semicolon    decimal value of ASCII code for semicolon (59)
!     equal        decimal value of ASCII code for equal sign (61)
!     question     decimal value of ASCII code for question mark (63)
!     atsign       decimal value of ASCII code for at sign (64)
!     backslash    decimal value of ASCII code for backslash (92)
!     caret        decimal value of ASCII code for caret (94)
!     underbar     decimal value of ASCII code for underbar (95)
!     vertical     decimal value of ASCII code for vertical bar (124)
!     tilde        decimal value of ASCII code for tilde (126)
!
!
module ascii
   implicit none
   integer null,tab
   integer linefeed,formfeed
   integer carriage,escape
   integer space,exclamation
   integer quote,pound
   integer dollar,percent
   integer ampersand
   integer apostrophe
   integer asterisk,plus
   integer comma,minus
   integer period,frontslash
   integer colon,semicolon
   integer equal,question
   integer atsign,backslash
   integer caret,underbar
   integer vertical,tilde
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

   function str_to_bool(input) result(res)
      character*(*), INTENT(IN) :: input
      character*1 :: res
      integer :: i,ntrue,nfalse
      integer, parameter :: bool_size = 10
      character*(bool_size) :: str,testval
      character*(*), parameter :: true_vals='yes       '&
         &//'y         '&
         &//'true      '&
         &//'t         '&
         &//'.true.    '
      character*(*), parameter :: false_vals='no        '&
         &//'n         '&
         &//'false     '&
         &//'f         '&
         &//'.false.   '

      str = trim(input)
      call lowcase(str)

      res="U"

      ntrue = len(true_vals)/bool_size
      do i=1,ntrue
         testval = true_vals((i-1)*bool_size+1:i*bool_size)
         if(str==testval) then
            res="T"
            return
         endif
      enddo

      nfalse = len(false_vals)/bool_size
      do i=1,nfalse
         testval = false_vals((i-1)*bool_size+1:i*bool_size)
         if(str==testval) then
            res="F"
            return
         endif
      enddo

   end function str_to_bool

   function str_to_logical(input) result(res)
      character*(*), INTENT(IN) :: input
      logical :: res
      character*1 :: bool

      bool = str_to_bool(input)
      if(bool=="T") then
         res=.true.
      else
         res=.false.
         if(bool=="U") then
            write(0,*) "Warning: input is not a valid logical value"
            write(0,*) "         input string: ",input
         endif
      endif

   end function str_to_logical

end
