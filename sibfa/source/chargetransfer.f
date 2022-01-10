c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #########################################################################################
c     ##                                                                                     ##
c     ##  module chargetranfer  --  charge transfer caracteristics for the current structure ##
c     ##                                                                                     ##
c     #########################################################################################
c
c
c     nlp     total number of lone pairs of the current structure
c     naccept total number of acceptors of the current structure
c     nlplocnl localnl number of lone pairs in the system 
c     nlploc   local number of lone pairs in the system
c     nlpbloc   local+neighbors number of lone pairs in the system
c     nacceptlocnl localnl number of acceptors in the system 
c     nacceptloc local number of acceptors in the system 
c     nacceptbloc local+neighbors number of acceptors in the system 
c     lpatom   indexes of the carrier atom for each lone pair
c     lplocnl  global-localnl correspondance for lone pairs
c     acceptlocnl  global-localnl correspondance for acceptors
c     acceploc  global-local correspondance for acceptors
c     acceptor  atom pairs defining electron acceptors
c     izlp     index of the 'iz' atom used to rotate lone pairs 
c     ixlp     index of the 'ix' atom used to rotate lone pairs 
c     ilplst   atom-lone pair correspondance 
c     nilplst  number of lone pairs carried by each atom
c     lonepair  geometrical data of lone pairs in the local frame 
c     rlonepair geometrical data of lone pairs in the global frame
c     rvdwct1   ct-vdw radius for electron donnors of the system     
c     rvdwct2   ct-vdw radius for electron acceptors of the system     
c     tasct  coefficients to compute integrals for the current structure
c     tapct  coefficients to compute integrals for the current structure
c     mact  coefficients to compute integrals for the current structure
c     aect  ionization potential for the current structure
c     ahct  electronic affinity for the current structure
c     dincr_lpect  lone pair vdw increment for the current structure
c     hybrid_lp     hybridization coefficients for the lone pairs of the current structure
c     nbaccept number of acceptors before each atom
c     namelp    name of the corresponding atom for each lone pair
c
c     lpcharge  partial charges of the lone pairs
c     etact    eta parameter for ect MRW formula
c     Q        Q parameter for ect MRW formula
c     D1        D1 parameter for ect MRW formula
c     D2        D2 parameter for ect MRW formula
c     ialp : ionization potential of the lone pair hybrid : 12.78 ev for sp3-hybrid and 10.88 ev for a sp2-hybrid
c      
c          
      module chargetransfer
      implicit none
      integer nlp,naccept,nlplocnl,nacceptlocnl
      integer nlploc,nlpbloc,nacceptloc,nacceptbloc
      integer, allocatable :: lpatom(:),lplocnl(:),acceptlocnl(:)
      integer, allocatable :: acceptloc(:) 
      integer, allocatable :: acceptor(:,:),izlp(:),ixlp(:)
      integer, allocatable :: ilplst(:,:),nilplst(:)
      real*8, allocatable :: lonepair(:,:),lpcharge(:),rlonepair(:,:)
      real*8, allocatable :: rvdwct1(:),rvdwct2(:)
      real*8, allocatable :: tasct(:,:),tapct(:,:)
      real*8, allocatable :: mact(:,:),ahct(:),aelct(:)
      real*8, allocatable :: dincr_lpect(:)
      real*8, allocatable :: hybrid_lp(:,:)
      integer, allocatable :: nbaccept(:)
      character*8, allocatable :: namelp(:)
      real*8 etaect,Q,D1,D2,ialp
      parameter (etaect = 9.4146d0)
      parameter (Q = 1.0d0)
      parameter (D1 = 1.0d0)
      parameter (D2 = 1.0d0)
      save
      end
