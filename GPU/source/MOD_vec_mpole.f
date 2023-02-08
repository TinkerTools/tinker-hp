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
#include "tinker_macro.h"
      !PGI$ ATTRIBUTES ALIGN:64:: imscalevec
      integer imscalevec(maxelst)
      !PGI$ ATTRIBUTES ALIGN:64:: r2vec1
      real(t_p),dimension(maxelst)::r2vec1
      !PGI$ ATTRIBUTES ALIGN:64:: invr2vec1
      real(t_p),dimension(maxelst)::invr2vec1
      !PGI$ ATTRIBUTES ALIGN:64:: xposvec2,yposvec2,zposvec2
      real(t_p),dimension(maxelst)::xposvec2,yposvec2,zposvec2
      !PGI$ ATTRIBUTES ALIGN:64:: dikvec
      real(t_p),dimension(maxelst)::dikvec
      !PGI$ ATTRIBUTES ALIGN:64:: dikvecx, dikvecy, dikvecz
      real(t_p),dimension(maxelst)::dikvecx,dikvecy,dikvecz
      !PGI$ ATTRIBUTES ALIGN:64:: dirvecx, dirvecy, dirvecz
      real(t_p),dimension(maxelst)::dirvecx,dirvecy,dirvecz
      !PGI$ ATTRIBUTES ALIGN:64:: dkrvecx, dkrvecy,dkrvecz
      real(t_p),dimension(maxelst)::dkrvecx,dkrvecy,dkrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: dikrvec,diqrkvec,dkqrivec
      real(t_p),dimension(maxelst)::dikrvec,diqrkvec,dkqrivec
      !PGI$ ATTRIBUTES ALIGN:64:: diqkvecx, diqkvecy, diqkvecz
      real(t_p),dimension(maxelst)::diqkvecx, diqkvecy,diqkvecz
      !PGI$ ATTRIBUTES ALIGN:64:: diqkrvecx, diqkrvecy, diqkrvecz
      real(t_p),dimension(maxelst)::diqkrvecx, diqkrvecy,diqkrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: dkqivecx, dkqivecy, dkqivecz
      real(t_p),dimension(maxelst)::dkqivecx, dkqivecy,dkqivecz
      !PGI$ ATTRIBUTES ALIGN:64:: dkqirvecx, dkqirvecy, dkqirvecz
      real(t_p),dimension(maxelst)::dkqirvecx, dkqirvecy,dkqirvecz
      !PGI$ ATTRIBUTES ALIGN:64:: dqiqkvecx, dqiqkvecy, dqiqkvecz
      real(t_p),dimension(maxelst)::dqiqkvecx, dqiqkvecy,dqiqkvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qikvec,qrrikvec
      real(t_p),dimension(maxelst)::qikvec,qrrikvec
      !PGI$ ATTRIBUTES ALIGN:64:: qikrvecx, qikrvecy, qikrvecz
      real(t_p),dimension(maxelst)::qikrvecx, qikrvecy,qikrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qkirvecx, qkirvecy, qkirvecz
      real(t_p),dimension(maxelst)::qkirvecx, qkirvecy,qkirvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qrirvecx, qrirvecy, qrirvecz
      real(t_p),dimension(maxelst)::qrirvecx,qrirvecy,qrirvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qrkrvecx, qrkrvecy, qrkrvecz
      real(t_p),dimension(maxelst)::qrkrvecx,qrkrvecy,qrkrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qrrvecx, qrrvecy, qrrvecz
      real(t_p),dimension(maxelst)::qrrvecx,qrrvecy,qrrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qikrrvecx, qikrrvecy, qikrrvecz
      real(t_p),dimension(maxelst)::qikrrvecx,qikrrvecy,qikrrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: qkirrvecx, qkirrvecy, qkirrvecz
      real(t_p),dimension(maxelst)::qkirrvecx,qkirrvecy,qkirrvecz
      !PGI$ ATTRIBUTES ALIGN:64:: mscalevec
      real(t_p),dimension(maxelst)::mscalevec
      !PGI$ ATTRIBUTES ALIGN:64:: evec, devec
      real(t_p),dimension(maxelst)::evec,devec
      !PGI$ ATTRIBUTES ALIGN:64:: frcvecx,frcvecy,frcvecz
      real(t_p),dimension(maxelst)::frcvecx,frcvecy,frcvecz
      !PGI$ ATTRIBUTES ALIGN:64:: ttmivecx,ttmivecy,ttmivecz
      real(t_p),dimension(maxelst)::ttmivecx,ttmivecy,ttmivecz
      !PGI$ ATTRIBUTES ALIGN:64:: ttmkvecx,ttmkvecy,ttmkvecz
      real(t_p),dimension(maxelst)::ttmkvecx,ttmkvecy,ttmkvecz

      end module
