c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ###################################################################
c     ##                                                               ##
c     ##  module  vec_mpole  --    Stuff for vectorization purpose     ##
c     ##                                                               ##
c     ###################################################################

      module vec_mpole
      use sizes
      use couple
      use vec
      use vec_elec
      implicit none
      !DIR$ ATTRIBUTES ALIGN:64:: r2vec1
      real*8 r2vec1(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: invr2vec1
      real*8 invr2vec1(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: xposvec2,yposvec2,zposvec2
      real*8 xposvec2(maxelst),yposvec2(maxelst),zposvec2(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dikvec
      real*8 dikvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dikvecx, dikvecy, dikvecz
      real*8 dikvecx(maxelst),dikvecy(maxelst), dikvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dirvecx, dirvecy, dirvecz
      real*8 dirvecx(maxelst),dirvecy(maxelst), dirvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dkrvecx, dkrvecy, dkrvecz
      real*8 dkrvecx(maxelst),dkrvecy(maxelst), dkrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dikrvec,diqrkvec,dkqrivec
      real*8 dikrvec(maxelst),diqrkvec(maxelst),dkqrivec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: diqkvecx, diqkvecy, diqkvecz
      real*8 diqkvecx(maxelst), diqkvecy(maxelst), diqkvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: diqkrvecx, diqkrvecy, diqkrvecz
      real*8 diqkrvecx(maxelst), diqkrvecy(maxelst), diqkrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dkqivecx, dkqivecy, dkqivecz
      real*8 dkqivecx(maxelst), dkqivecy(maxelst), dkqivecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dkqirvecx, dkqirvecy, dkqirvecz
      real*8 dkqirvecx(maxelst), dkqirvecy(maxelst), dkqirvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: dqiqkvecx, dqiqkvecy, dqiqkvecz
      real*8 dqiqkvecx(maxelst), dqiqkvecy(maxelst), dqiqkvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qikvec,qrrikvec
      real*8 qikvec(maxelst),qrrikvec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qikrvecx, qikrvecy, qikrvecz
      real*8 qikrvecx(maxelst), qikrvecy(maxelst), qikrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qkirvecx, qkirvecy, qkirvecz
      real*8 qkirvecx(maxelst), qkirvecy(maxelst), qkirvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrirvecx, qrirvecy, qrirvecz
      real*8 qrirvecx(maxelst),qrirvecy(maxelst), qrirvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrkrvecx, qrkrvecy, qrkrvecz
      real*8 qrkrvecx(maxelst),qrkrvecy(maxelst), qrkrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qrrvecx, qrrvecy, qrrvecz
      real*8 qrrvecx(maxelst),qrrvecy(maxelst), qrrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qikrrvecx, qikrrvecy, qikrrvecz
      real*8 qikrrvecx(maxelst),qikrrvecy(maxelst), qikrrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: qkirrvecx, qkirrvecy, qkirrvecz
      real*8 qkirrvecx(maxelst),qkirrvecy(maxelst), qkirrvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: mscalevec
      real*8 mscalevec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: evec, devec
      real*8 evec(maxelst),devec(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: frcvecx,frcvecy,frcvecz
      real*8 frcvecx(maxelst),frcvecy(maxelst),frcvecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: ttmivecx,ttmivecy,ttmivecz
      real*8 ttmivecx(maxelst),ttmivecy(maxelst),ttmivecz(maxelst)
      !DIR$ ATTRIBUTES ALIGN:64:: ttmkvecx,ttmkvecy,ttmkvecz
      real*8 ttmkvecx(maxelst),ttmkvecy(maxelst),ttmkvecz(maxelst)

      end module
