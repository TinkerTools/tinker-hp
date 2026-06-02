#ifndef DAMPING_H
#define DAMPING_H
#include "utils.h"

using namespace CodePrm;

template<int rorder> __device__
inline static void dampewald(const real& r, const real& r2, const real& ewald, const real& scale, restrict real* dmpe){
   const real ra      = ewald* r;
   const real exp2a   = f_exp(-ra*ra);
         real bni_1   = f_erfc(ra) / r;
              dmpe[0] = scale* bni_1;
   const real aesq2   = twor* ewald* ewald;
         real afac    = ewald>zeror ? f_inv(sqrtpi*ewald) : zeror;
   const real ir2     = f_inv(r2);
   const int  niter   = (rorder-1)>>1;
   #pragma unroll
   for( int i=0;i<niter;i++ ){
      afac      = aesq2 * afac;
      real bni  = (static_cast<real>(2*i+1)*bni_1+afac*exp2a) * ir2;
      dmpe[i+1] = scale * bni;
      bni_1     = bni;
   }
}


template<int rorder, int fea> __device__
inline static void dampthole(const real r, real& damp, const real pgamma, real* dmpik){
   const real p5 = 0.5; const real p65 = 0.65; const real p15 = 0.15;
   const real r9 = 9.0; const real r18 = 18.0; const real r35 = 35.0;
   dmpik[0] = oner; dmpik[1] = oner; dmpik[2] = oner; if (rorder>9) dmpik[3] = oner;

   //use alternate Thole model for AMOEBA+ direct polarization
   if (fea & DIRDAMP){
      if (damp!=zeror && pgamma!=zeror){
         damp = pgamma * f_pow(r/damp,1.5);
         if (damp<50.0){
            const real expdamp = f_exp(-damp);
                           dmpik[0] = oner - expdamp;
                           dmpik[1] = oner - expdamp*(oner +  p5*damp);
            if (rorder>=7) dmpik[2] = oner - expdamp*(oner + p65*damp + p15*damp*damp);
         }
      }
   }
   //use original AMOEBA Thole polarization damping factors
   else {
      if (damp!=zeror && pgamma!=zeror){
         const real r_damp=r/damp;
         damp = pgamma * r_damp*r_damp*r_damp;
         if (damp<50.0){
            const real expdamp = f_exp(-damp);
               dmpik[0] = oner - expdamp;
               dmpik[1] = oner - expdamp*(oner+damp);
            if (rorder>=7){
               const real damp2 = damp*damp;
               dmpik[2] = oner - expdamp*(oner +damp +0.6*damp2);
               if (rorder>=9) dmpik[3] = oner - expdamp*(oner +damp +(r18/r35)*damp2 +(r9/r35)*damp2*damp);
            }
         }
      }
   }
}
#endif
