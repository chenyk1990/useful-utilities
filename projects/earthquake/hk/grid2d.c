/*
grid2d.c:	2D grid-search
Description:	Search for minimum of function f(x,y) using grid-search

revision:	Lupei Zhu, 12-20-98, initial coding
		Lupei Zhu, 9-10-03, add covarance sig_xy

*/

/*#include <stdio.h>
#include <math.h>
#include <float.h>*/
#include "gmt.h"

float	grid2d(	float	*f,	/* IN: f(ix+iy*nx) at (ix,iy), ix faster */
		int	nx,	/* IN: dimension of x, nx > 2*/
		int	ny,	/* IN: dimension of y, ny > 2*/
		float	*x,	/* OUT: x of minimum */
		float	*sig_x,	/* OUT: D(f,x,2) ~ 2/sigma^2 of x */
		float	*y,	/* OUT: y of minimum */
		float	*sig_y,	/* OUT: D(f,y,2) ~ 2/sigma^2 of y */
		float 	*sig_xy, /* OUT: D(f,x,y) ~ 2/sigma_xy */
		int	*ixy,	/* OUT: location of the minimum */
		int	*flag	/* OUT: >0 if minimum on boundaries */
	)
{
     int	i, ix, iy, best_ix, best_iy;
     float	f0, eps;

     f0 = FLT_MAX;
     for(i=0,iy=0; iy<ny; iy++) {
        for(ix=0; ix<nx; ix++,i++) {
	   if( !GMT_is_fnan (f[i]) && f[i]<f0 ) {
		f0 = f[i];
		best_ix = ix;
		best_iy = iy;
	   }
        }
     }

     *ixy = i = best_ix + nx*best_iy;

     if (best_ix==0 || best_ix==nx-1 || best_iy==0 || best_iy==ny-1
	|| GMT_is_fnan (f[i-1]) || GMT_is_fnan (f[i+1])
	|| GMT_is_fnan (f[i-nx]) || GMT_is_fnan (f[i+nx])
	|| GMT_is_fnan (f[i-nx-1]) || GMT_is_fnan (f[i+nx-1])
	|| GMT_is_fnan (f[i-nx+1]) || GMT_is_fnan (f[i+nx+1]) ) {
	fprintf(stderr,"Warning: minium on boundaries\n");
        *flag = 1;
	*x = best_ix; *y= best_iy;
	*sig_x = *sig_y = 1; *sig_xy = 0.;
	return f0;
     }
     *flag = 0;
     *sig_x = f[i+1]+f[i-1]-2*f0; /*(f[i+nx+1]+f[i+nx-1]-2*f[i+nx]+f[i-nx+1]+f[i-nx-1]-2*f[i-nx]);*/
     if (*sig_x<0.) {fprintf(stderr,"Warning: p2f/px2=%e < 0\n",*sig_x);*sig_x=FLT_MIN;}
     *sig_y = f[i+nx]+f[i-nx]-2*f0; /*(f[i+nx+1]+f[i-nx+1]-2*f[i+1]+f[i+nx-1]+f[i-nx-1]-2*f[i-1]);*/
     if (*sig_y<0.) {fprintf(stderr,"Warning: p2f/py2=%e < 0\n",*sig_y);*sig_y=FLT_MIN;}
     /**sig_xy = 0.5*(f[i+nx+1]+2*f0+f[i-nx-1]-f[i+nx]-f[i+1]-f[i-1]-f[i-nx]);*/
     *sig_xy = 0.25*(f[i+nx+1]+f[i-nx-1]-f[i+nx-1]-f[i-nx+1]);

     eps = 4*(*sig_x)*(*sig_y)-(*sig_xy)*(*sig_xy);
     *x = best_ix - (2*(*sig_y)*(f[i+1]-f[i-1])-(*sig_xy)*(f[i+nx]-f[i-nx]))/eps;
     *y = best_iy - (2*(*sig_x)*(f[i+nx]-f[i-nx])-(*sig_xy)*(f[i+1]-f[i-1]))/eps;
     f0 -= 0.5*(*x-best_ix)*(*x-best_ix)*(*sig_x)
      	    +0.5*(*y-best_iy)*(*y-best_iy)*(*sig_y)
	    +(*x-best_ix)*(*y-best_iy)*(*sig_xy);
   
     return f0;

}

/* for an interior minimum:
     for(iy=1; iy<ny-1; iy++) {
        for(i=1+iy*nx,ix=1; ix<nx-1; ix++,i++) {
	  if(err[i]<err0 &&
	     err[i]<err[i-1] &&
	     err[i]<err[i+1] &&
	     err[i]<err[i-nx] &&
	     err[i]<err[i+nx] &&
	     err[i]<err[i-nx+1] &&
	     err[i]<err[i-nx-1] &&
	     err[i]<err[i+nx-1] &&
	     err[i]<err[i+nx+1]) {
		err0 = err[i];
		best_ix = ix;
		best_iy = iy;
	  }
        }
     }
*/
