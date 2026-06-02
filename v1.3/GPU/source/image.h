#ifndef IMAGE_H
#define IMAGE_H

#ifndef ORTHOGONAL_BOX_SHAPE_ONLY
#define INC_ALL_BOX_SHAPE
#endif

#define TINKER_IMAGE_PARAMS                                  \
    const real xcell , const real ycell , const real zcell , \
    const real xcell2, const real ycell2, const real zcell2
#define TINKER_IMAGE1_PARAMS                                 \
    TINKER_IMAGE_PARAMS , const real beta_sin, const real beta_cos, const real beta_term, const real gamma_sin, const real gamma_cos, const real gamma_term, const int monoclinic, const int triclinic, const int octahedron, const real _box34
#define TINKER_PROC_BOX_PARAMS                               \
    const real p_xbeg, const real p_xend, const real p_ybeg, \
    const real p_yend, const real p_zbeg, const real p_zend
#define MIDPOINTIMAGE_PARAMS                                 \
    TINKER_IMAGE_PARAMS , TINKER_PROC_BOX_PARAMS
#define MIDPOINTIMAGE1_PARAMS                                \
    TINKER_IMAGE1_PARAMS , TINKER_PROC_BOX_PARAMS , const real eps_cell

#define TINKER_IMAGE_ARGS                                    \
    xcell,ycell,zcell,xcell2,ycell2,zcell2,beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term,monoclinic,triclinic,octahedron,_box34
#define TINKER_PROC_BOX_ARGS                                 \
    p_xbeg,p_xend,p_ybeg,p_yend,p_zbeg,p_zend
#define MIDPOINTIMAGE_ARGS                                   \
    TINKER_IMAGE_ARGS , TINKER_PROC_BOX_ARGS , eps_cell

#define Image(x,y,z) image_generic(x,y,z,TINKER_IMAGE_ARGS)
#define ImageOrth(x,y,z) image_orthogonal(x,y,z,TINKER_IMAGE_ARGS)
#define Midpointimage(xk,yk,zk,xr,yr,zr) midpointimage(xk,yk,zk,xr,yr,zr,MIDPOINTIMAGE_ARGS)

enum BoxShape_e { BOX_ORTH, BOX_MONO, BOX_TRIC, BOX_OCTA };

__device__
inline static void cart_to_frac( real& xr, real &yr, real &zr, TINKER_IMAGE1_PARAMS )
{
   zr =                      zr            / gamma_term;
   yr = (     yr           - zr*beta_term) / gamma_sin;
   xr =  xr - yr*gamma_cos - zr*beta_cos;
}

__device__
inline static void frac_to_cart( real& xr, real &yr, real &zr, TINKER_IMAGE1_PARAMS )
{
   xr = xr + yr*gamma_cos + zr*beta_cos;
   yr =      yr*gamma_sin + zr*beta_term;
   zr =                     zr*gamma_term;
}

__device__ 
inline static void image_orthogonal ( real& xr, real &yr, real &zr, TINKER_IMAGE1_PARAMS )
{
   if ( abs(xr)>xcell2 ) xr -= copysign(xcell,xr)*f_floor((abs(xr)-xcell2)/xcell+(real)1.0);
   if ( abs(yr)>ycell2 ) yr -= copysign(ycell,yr)*f_floor((abs(yr)-ycell2)/ycell+(real)1.0);
   if ( abs(zr)>zcell2 ) zr -= copysign(zcell,zr)*f_floor((abs(zr)-zcell2)/zcell+(real)1.0);
}

__device__ 
inline static void image_generic( real& xr, real &yr, real &zr, TINKER_IMAGE1_PARAMS )
{
#ifdef INC_ALL_BOX_SHAPE
   if ( monoclinic||triclinic ) cart_to_frac( xr,yr,zr,TINKER_IMAGE_ARGS );
#endif

   ImageOrth( xr,yr,zr );

#ifdef INC_ALL_BOX_SHAPE
   if ( monoclinic||triclinic ) {
      frac_to_cart( xr,yr,zr,TINKER_IMAGE_ARGS );
   }
   else if (octahedron) {
      if (abs(xr)+abs(yr)+abs(zr) > _box34) {
         xr -= copysign(xcell2,xr);
         yr -= copysign(ycell2,yr);
         zr -= copysign(zcell2,zr);
      }
   }
#endif
}

__device__
inline static int midpointimage ( real& restrict xk, real& restrict yk, real& restrict zk, real& restrict xr, real& restrict yr, real& restrict zr, MIDPOINTIMAGE1_PARAMS )
{
#ifdef INC_ALL_BOX_SHAPE
   if ( monoclinic||triclinic ){
      cart_to_frac( xr,yr,zr,TINKER_IMAGE_ARGS );
      cart_to_frac( xk,yk,zk,TINKER_IMAGE_ARGS );
   }
#endif

   image_orthogonal(xr,yr,zr,TINKER_IMAGE_ARGS);

#ifdef INC_ALL_BOX_SHAPE
   if (octahedron) {
      if (abs(xr)+abs(yr)+abs(zr) > _box34) {
         xr -= copysign(xcell2,xr);
         yr -= copysign(ycell2,yr);
         zr -= copysign(zcell2,zr);
      }
   }
#endif

   /* definition of middle point between i and k atoms */
   xk += 0.5*xr;
   yk += 0.5*yr;
   zk += 0.5*zr;

   image_orthogonal(xk,yk,zk,TINKER_IMAGE_ARGS);

#ifdef INC_ALL_BOX_SHAPE
   if ( monoclinic||triclinic ) frac_to_cart( xr,yr,zr,TINKER_IMAGE_ARGS );
   else if (octahedron) {
      if (abs(xk)+abs(yk)+abs(zk) > _box34) {
         xk -= copysign(xcell2,xr);
         yk -= copysign(ycell2,yr);
         zk -= copysign(zcell2,zr);
      }
   }
#endif

   /* Adjust mid point position if necessary */
   if ((xcell2-xk)<eps_cell) xk -= 4*eps_cell;
   if ((ycell2-yk)<eps_cell) yk -= 4*eps_cell;
   if ((zcell2-zk)<eps_cell) zk -= 4*eps_cell;

   /* Test if in the processor box */
   return ((zk<p_zbeg)||(zk>=p_zend)
         ||(yk<p_ybeg)||(yk>=p_yend)
         ||(xk<p_xbeg)||(xk>=p_xend))?  0: 1;
}

#ifndef IMAGE_CC
extern const real eps_cell;
extern const real _xcell,_ycell,_zcell;
extern const real _ixcell,_iycell,_izcell;
extern const real _box34;
extern const real beta_sin,beta_cos,beta_term,gamma_sin,gamma_cos,gamma_term;
extern const int monoclinic,triclinic,octahedron;
#endif
#endif
