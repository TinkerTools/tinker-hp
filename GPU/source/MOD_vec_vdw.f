c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec_elec  --    Stuff for vectorization purpose      ##
c     ##                                                               ##
c     ###################################################################
c
c
#include "tinker_macro.h"
      module vec_vdw
      use sizes
      use couple
c     use domdec
      implicit none
!PGI$ ATTRIBUTES ALIGN:64::kvvec,kvvec1,kvvec2
      integer,dimension(maxvlst):: kvvec,kvvec1,kvvec2
!PGI$ ATTRIBUTES ALIGN:64:: kvlocvec,kvlocvec1,kvlocvec2
      integer,dimension(maxvlst):: kvlocvec,kvlocvec1,kvlocvec2
!PGI$ ATTRIBUTES ALIGN:64:: ktvec,ktvec2
      integer,dimension(maxvlst):: ktvec,ktvec2
!PGI$ ATTRIBUTES ALIGN:64::ivscalevec
      integer,dimension(maxvlst):: ivscalevec
!PGI$ ATTRIBUTES ALIGN:64::vscalevec
      real(t_p),dimension(maxvlst):: vscalevec
!PGI$ ATTRIBUTES ALIGN:64::radminvec
      real(t_p),dimension(maxvlst):: radminvec
!PGI$ ATTRIBUTES ALIGN:64::radmin4vec
      real(t_p),dimension(maxvlst):: radmin4vec
!PGI$ ATTRIBUTES ALIGN:64::epsilonvec
      real(t_p),dimension(maxvlst):: epsilonvec
!PGI$ ATTRIBUTES ALIGN:64::epsilon4vec
      real(t_p),dimension(maxvlst):: epsilon4vec
!PGI$ ATTRIBUTES ALIGN:64::rvvec,epsvec,rvvec2
      real(t_p),dimension(maxvlst):: rvvec,rvvec2,epsvec
!PGI$ ATTRIBUTES ALIGN:64::rv7vec,invrhovec,epsvec2
      real(t_p),dimension(maxvlst):: rv7vec,invrhovec,epsvec2
!PGI$ ATTRIBUTES ALIGN:64::rik2vec,rik6vec,invrikvec
      real(t_p),dimension(maxvlst):: rik2vec,rik6vec,invrikvec
!PGI$ ATTRIBUTES ALIGN:64::rik2vec1
      real(t_p),dimension(maxvlst):: rik2vec1
!PGI$ ATTRIBUTES ALIGN:64::rikvec,rik3vec,rik4vec
      real(t_p),dimension(maxvlst):: rikvec,rik3vec,rik4vec
!PGI$ ATTRIBUTES ALIGN:64::rik7vec,invtmpvec
      real(t_p),dimension(maxvlst):: rik7vec,invtmpvec
!PGI$ ATTRIBUTES ALIGN:64::tapervec,evec,rik5vec
      real(t_p),dimension(maxvlst):: tapervec,evec,rik5vec
!PGI$ ATTRIBUTES ALIGN:64:: devec,dtapervec,rv7orhovec
      real(t_p),dimension(maxvlst):: devec,dtapervec,rv7orhovec
!PGI$ ATTRIBUTES ALIGN:64:: dedxvec,dedyvec,dedzvec
      real(t_p),dimension(maxvlst):: dedxvec, dedyvec, dedzvec
!PGI$ ATTRIBUTES ALIGN:64:: redkvec,dtauvec
      real(t_p),dimension(maxvlst):: redkvec,dtauvec
!PGI$ ATTRIBUTES ALIGN:64:: tauvec,gtauvec,tau7vec
      real(t_p),dimension(maxvlst):: tauvec,gtauvec,tau7vec
!PGI$ ATTRIBUTES ALIGN:64:: redivec,p6vec
      real(t_p),dimension(maxvlst):: redivec,p6vec
!PGI$ ATTRIBUTES ALIGN:64:: devxvec,devyvec,devzvec
      real(t_p),dimension(maxvlst):: devxvec,devyvec,devzvec
!PGI$ ATTRIBUTES ALIGN:64:: devvxvec,devvyvec,devvzvec
      real(t_p),dimension(maxvlst):: devvxvec,devvyvec,devvzvec

      end
