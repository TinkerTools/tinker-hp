!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ###############################################################
!     ##                                                           ##
!     ##  module restart  --  path-explicit DYN persistence codec  ##
!     ##                                                           ##
!     ###############################################################
!
#include "tinker_macro.h"
module restart
   use status  ,only: status_t
   implicit none

   interface
      module subroutine write_restart(filename,title,xbox,ybox,zbox,alpha,beta,&
         &gamma,x,y,z,v,a,aalt,status_value)
         character(*) ,intent(in) :: filename,title
         real(r_p)   ,intent(in) :: xbox,ybox,zbox
         real(t_p)   ,intent(in) :: alpha,beta,gamma
         real(r_p)   ,intent(in) :: x(:),y(:),z(:),v(:,:),a(:,:),aalt(:,:)
         type(status_t),intent(out) :: status_value
      end subroutine write_restart

      module subroutine read_restart(filename,expected_atoms,xbox,ybox,zbox,&
         &alpha,beta,gamma,x,y,z,v,a,aalt,aalt2,aalt2_present,status_value)
         character(*) ,intent(in) :: filename
         integer      ,intent(in) :: expected_atoms
         real(r_p)    ,intent(out) :: xbox,ybox,zbox
         real(t_p)    ,intent(out) :: alpha,beta,gamma
         real(r_p)    ,intent(inout) :: x(:),y(:),z(:),v(:,:),a(:,:),aalt(:,:),&
            &aalt2(:,:)
         logical      ,intent(out) :: aalt2_present
         type(status_t),intent(out) :: status_value
      end subroutine read_restart
   end interface
end module restart
