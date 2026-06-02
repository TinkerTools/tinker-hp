#define IMAGE_CC

#include "utils.h"
#include "image.h"

real _xcell, _ycell, _zcell;
real _ixcell, _iycell, _izcell;
real eps_cell;
real beta_sin, beta_cos, beta_term, gamma_sin, gamma_cos, gamma_term;
real _box34;
int monoclinic=0,triclinic=0,octahedron=0;
BoxShape_e boxShape;
//extern const int tinkerdebug;

long long gcd(long long a, long long b)
{
   while (b != 0) {
      long long t = b;
      b = a % b;
      a = t;
   }
   return llabs(a);
}

EXTERN_C_BEG

void C_get_cell( real xcell_, real ycell_, real zcell_, real eps_cell_, real beta_sin_, real beta_cos_,
                 real beta_term_, real gamma_sin_, real gamma_cos_, real gamma_term_,
                 int mono_, int tric_, int octa_, real box34_ ){

   assert( mono_+tric_+octa_<2 );

   _xcell     = xcell_;
   _ycell     = ycell_;
   _zcell     = zcell_;
   _ixcell    = (real) CodePrm::onem/ (double)xcell_;
   _iycell    = (real) CodePrm::onem/ (double)ycell_;
   _izcell    = (real) CodePrm::onem/ (double)zcell_;
   eps_cell   = eps_cell_;
   beta_sin   = beta_sin_;
   beta_cos   = beta_cos_;
   beta_term  = beta_term_;
   gamma_sin  = gamma_sin_;
   gamma_cos  = gamma_cos_;
   gamma_term = gamma_term_;
   monoclinic = mono_;
   triclinic  = tric_;
   octahedron = octa_;
   _box34     = box34_;

                boxShape = BOX_ORTH;
   if ( mono_ ) boxShape = BOX_MONO;
   if ( tric_ ) boxShape = BOX_TRIC;
   if ( octa_ ) boxShape = BOX_OCTA;

   //if (tinkerdebug) {
   //   printf(" C_get_cell -beta %lf %lf -gamma %lf %lf -shape %d %d %d\n", beta_sin,beta_term,gamma_sin,gamma_term,monoclinic,triclinic,octahedron,_box34);
   //}

}

void best_rational_approx(const double x, const long long max_den,
                          long long *num, long long *den)
{
    long long h0 = 0, h1 = 1;
    long long k0 = 1, k1 = 0;

    double value = x;
    int sign = (x < 0) ? -1 : 1;
    value = fabs(value);

    while (1)
    {
        long long a = (long long)floor(value);

        long long h2 = a * h1 + h0;
        long long k2 = a * k1 + k0;

        if (k2 > max_den)
           break;

        h0 = h1; h1 = h2;
        k0 = k1; k1 = k2;

        double frac = value - a;
        if (frac < 1e-16)
            break;

        value = 1.0 / frac;
    }

    *num = sign * h1;
    *den = k1;

    long long g = gcd(*num, *den);
    *num /= g;
    *den /= g;
    //printf("  -best_rational_approx %24.15f -m_den %lu = ( %lu / %lu )\n",x,max_den,*num,*den);
}

EXTERN_C_END
