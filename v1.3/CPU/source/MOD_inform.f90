!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  module inform  --  control values for I/O and program flow  ##
!     ##                                                              ##
!     ##################################################################
!
!
module inform
   implicit none
   integer :: digits !<decimal places output for energy and coordinates
   integer :: iprint !<steps between status printing (0=no printing)
   integer :: iwrite !<steps between coordinate dumps (0=no dumps)
   integer :: isend !<steps between socket communication (0=no sockets)
   logical :: silent !<logical flag to turn off all information printing
   logical :: verbose !<logical flag to turn on extra information printing
   logical :: dd_verbose !<logical flag to turn on extra information printing about dd
   logical :: verbosepmetime !<logical flag to turn on extra timings of PME printing
   logical :: verboseforcestime !<logical flag to turn on extra timings of forces
   logical :: debug !<logical flag to turn on full debug printing
   logical :: holdup !<logical flag to wait for carriage return on exit
   logical :: abort !<logical flag to stop execution at next chance

   ! All program list
   enum,bind(C)
      enumerator analyze_a
      enumerator bar_a
      enumerator dynamic_a
      enumerator minimize_a
      enumerator testgrad_a
      enumerator pimd_a
   end enum
   integer:: app_id=dynamic_a ! Only to be modifed inside a program

   ! Debug static counter
   enum,bind(C)
      enumerator tindPath,tindForce,tindEnergy,tindAtom
   end enum
   logical deb_Path,deb_Force,deb_Energy,deb_Atom


   ! Inform separated Subroutines
   interface
      module subroutine initDebugEnv
      end subroutine
   end interface

   interface
      module subroutine info_minmax_pva(opt)
         integer,optional::opt
      end subroutine
   end interface
   interface
      module subroutine info_dyn()
      end subroutine
   end interface

   interface minmaxone
      module subroutine minmaxonei( vector,sz,name )
         implicit none
         integer sz
         integer vector(*)
         character(*),optional,intent(in)::name
      end subroutine
      module subroutine minmaxonet( vector,sz,name )
         implicit none
         integer sz
         real*8 vector(*)
         character(*),optional,intent(in)::name
      end subroutine
   end interface

   interface normp
      module subroutine normt( array,n,val,p )
         integer n
         integer,optional::p
         real*8 array(*)
         real*8 val
      end subroutine
   end interface
end
