c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###############################################################
c     ##                                                           ##
c     ##  neigh.i  --  pairwise neighbor list indices and storage  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lbuffer     width of the neighbor list buffer region
c     pbuffer     width of the preconditioner list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     pbuf2       square of half the preconditioner list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     vbufx       square of vdw cutoff plus twice the list buffer
c     mbufx       square of multipole cutoff plus twice the list buffer
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     dovlst      logical flag to rebuild vdw neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c
c
      integer, pointer :: nvlst(:)
      integer, pointer :: vlst(:,:)
      integer, pointer :: nelst(:)
      integer, pointer :: elst(:,:)
      integer, pointer :: ineignl(:)
      integer, pointer :: neigcell(:,:),numneigcell(:),repartcell(:)
      integer, pointer :: cell_len(:),indcell(:),bufbegcell(:)
      real*8, pointer :: xbegcell(:),ybegcell(:),zbegcell(:)
      real*8, pointer :: xendcell(:),yendcell(:),zendcell(:)

      integer ineigup
      integer ncell_tot
      real*8 lbuffer,pbuffer
      real*8 lbuf2,pbuf2
      real*8 vbuf2
      real*8 mbuf2
      real*8 torquebuf2
      logical dovlst,doclst
      logical domlst
      common /neigh/ lbuffer,pbuffer,lbuf2,pbuf2,vbuf2,mbuf2,
     &               torquebuf2,ineigup,
     &               nvlst,vlst,nelst,elst,
     &               dovlst,
     &               domlst,ineignl,ncell_tot,
     &               xbegcell,ybegcell,zbegcell,
     &               xendcell,yendcell,zendcell,neigcell,numneigcell,
     &               cell_len,indcell,bufbegcell,repartcell
