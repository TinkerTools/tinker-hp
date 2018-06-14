c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ##############################################################
c     ##                                                          ##
c     ##  atmlst.i  --  local geometry terms involving each atom  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     bndlist   list of the bond numbers involving each atom
c     anglist   list of the angle numbers centered on each atom
c
c
c      integer bndlist,anglist
c      common /atmlst/ bndlist(maxval,maxatm),
c     &                anglist(maxval*(maxval-1)/2,maxatm)
      integer, pointer :: bndlist(:,:),anglist(:,:)
      integer, pointer :: bndglob(:),angleglob(:),torsglob(:)
      integer, pointer :: bitorsglob(:),strbndglob(:)
      integer, pointer :: ureyglob(:),angangglob(:)
      integer, pointer :: opbendglob(:),opdistglob(:)
      integer, pointer :: impropglob(:),imptorglob(:)
      integer, pointer :: pitorsglob(:),strtorglob(:)
      integer, pointer :: tortorglob(:)
      integer, pointer :: vdwglob(:),poleglob(:),polerecglob(:)
      integer, pointer :: vdwglobnl(:),poleglobnl(:)
      integer, pointer :: molculeglob(:)
      integer, pointer :: npfixglob(:),ndfixglob(:),nafixglob(:)
      integer, pointer :: ntfixglob(:),ngfixglob(:),nchirglob(:)
      common /atmlst/ bndlist,anglist,bndglob,angleglob,torsglob,
     $                bitorsglob,strbndglob,ureyglob,angangglob,
     $                opbendglob,opdistglob,impropglob,imptorglob,
     $                pitorsglob,strtorglob,tortorglob,vdwglob,
     $                poleglob,polerecglob,vdwglobnl,poleglobnl,
     $                molculeglob,npfixglob,ndfixglob,nafixglob,
     $                ntfixglob,ngfixglob,nchirglob
