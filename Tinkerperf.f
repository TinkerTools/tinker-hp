!----------------------------------------------------------------
!DIR$ ATTRIBUTES ALIGN:64::kglobvec1,kbisvec1
      integer kglobvec1(maxvlst),kbisvec1(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64:: kvvec1,kvlocvec1
      integer kvvec1(maxvlst),kvlocvec1(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rv7vec,rvvec2
      real*8 rv7vec(maxvlst),rvvec2(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rikvec,rik2vec,rik3vec
      real*8 rikvec(maxvlst),rik2vec(maxvlst),rik3vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rik4vec,rik5vec,rik6vec
      real*8 rik4vec(maxvlst),rik5vec(maxvlst),rik6vec(maxvlst)
!DIR$ ATTRIBUTES ALIGN:64::rik7vec,invrhovec,invtmpvec
      real*8 rik7vec(maxvlst),invrhovec(maxvlst),invtmpvec(maxvlst)
!----------------------------------------------------------------
nvloop8  = (int(nnvlst /  8) + 1) *  8
nvloop16 = (int(nnvlst / 16) + 1) * 16
nvloop8  = merge(nnvlst, nvloop8  , mod(nnvlst,  8) .eq. 0)
nvloop16 = merge(nnvlst, nvloop16 , mod(nnvlst, 16) .eq. 0)
!----------------------------------------------------------------
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
 do k = 1, nvloop16
    mask1(k)= (kvlocvec(k) /= 0).and.(kbisvec (k) <= nbloc)
&                               .and.(kvlocvec(k) <= nbloc)
 enddo
!----------------------------------------------------------------
         nnvlst1   = count(mask1)
         kglobvec1 = pack(kglobvec,mask1)
         kbisvec1  = pack(kbisvec ,mask1)
         kvvec1    = pack(kvvec   ,mask1)
         kvlocvec1 = pack(kvlocvec,mask1)
!----------------------------------------------------------------
         kk = 0
!DIR$ ASSUME (MOD(nvloop16,16).eq.0)
         do k = 1, nvloop16
            if (mask1(k)) then
               kk = kk + 1
               kglobvec1(kk) = kglobvec(k)
               kbisvec1 (kk) = kbisvec (k)
               kvvec1   (kk) = kvvec   (k)
               kvlocvec1(kk) = kvlocvec(k)
            endif
         enddo
         nnvlst1 = kk
!----------------------------------------------------------------
!DIR$ ASSUME (MOD(nvloop8,8).eq.0)
do k = 1, nvloop8
   rv7vec (k)   = rvvec2 (k)  ** 7
   rik3vec(k)   = rik2vec(k)  * rikvec(k)
   rik4vec(k)   = rik3vec(k)  * rikvec(k)
   rik5vec(k)   = rik4vec(k)  * rikvec(k)
   rik6vec(k)   = rik5vec(k)  * rikvec(k)
   rik7vec(k)   = rik6vec(k)  * rikvec(k)
   invrhovec(k) = (rik7vec(k) + ghal * rv7vec(k)) ** - one
   invtmpvec(k) = (rikvec (k) + dhal * rvvec2(k)) ** - one
enddo
