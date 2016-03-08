/*******************************************************
 -		decon_sub.c
 -  subroutines for decovolving source function from 
 -  teleseismic records
 -
 -  revision history:
 -  06/22/97	Lupei Zhu	Initial revision
 -				water-level is now c*auto-
 -				correlation, consistent with
 -				wiener filter methods
 -  10/28/98	Lupei Zhu	correct amplitude change
 -				due to c by *(1+c)
********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Complex.h"
#include "sac.h"

/******************************************************************
* calculate deconvolution filter 1/(s[w]+water_level).
* the water-level is c*auto-correlation(s[t]). A low-pass gaussian filter
*	exp(-w^2/(4 sigma^2))
* is also applied.
******************************************************************/
void	water(
		   complex	*a,		/* IN/OUT, spectrum of s(t) */
		   int		nft,		/* number of pts of s(w), 2^N*/
		   float	dt,		/* sampling interval in sec */
		   float	c,		/* water-level in %*/
		   float	gauss		/* sigma in Hz */
		   )
{
  int      j;
  float    w, delw;
  double   *d2, agg, water, u;

  u = 1+c;
  d2 = (double *) malloc ((nft+1)*sizeof(double));
  delw = PI/(dt*nft);
  d2[0] = a[0].x*a[0].x;
  d2[nft] = a[0].y*a[0].y;
  water = d2[0]+d2[nft];
  for (j=1; j<nft; j++) {
    d2[j] = a[j].x*a[j].x+a[j].y*a[j].y;
    water += d2[j];
  }
  water = c*water/nft;		/* in terms of auto-correlation */
  a[0].x *= u/(d2[0]+water);
  for (w=delw,j=1; j<nft; j++,w+=delw) {
    agg = 0.5*w/gauss;
    a[j] = dmltp(u*exp(-agg*agg)/(d2[j]+water),conjg(a[j]));
  }
  agg = 0.5*w/gauss;
  a[0].y *= u*exp(-agg*agg)/(d2[j]+water);
  free(d2);
}


/* solving toeplitz equations using Levinson recursion*/
void	toeplitz(
		 const float	*r,		/* In: autocorrelation */
		 float		*g,		/* In: cross-correlation */
		 /* Out: deconvolution result */
		 int		n,		/* In: length of filter */
		 float		u		/* In: white noise level */
		 )
{
  int	i, j, k;
  float	dd, bb, v, *a, *aa, *x, *xx;

  i = n*sizeof(float);
  a = (float *) malloc(i);
  x = (float *) malloc(i);
  aa = (float *) malloc(i);
  xx = (float *) malloc(i);

  v=(1.+u)*r[0];
  aa[0]=1.;
  xx[0]=g[0]/v;
  for (j=1; j<n; j++) {

    for (dd=0.,bb=0.,i=0; i<j; i++) {
      k=j-i;
      dd+=r[k]*aa[i];
      bb+=r[k]*xx[i];
    }

    a[0]=1.;
    a[j]=-dd/v;
    for (i=1; i<j; i++)
      a[i]=aa[i]+a[j]*aa[j-i];
    v+=dd*a[j];

    x[j]=(g[j]-bb)/v;
    for (i=0; i<j; i++)
      x[i]=xx[i]+x[j]*a[j-i];

    for (i=0; i<=j; i++) {
      aa[i]=a[i];
      xx[i]=x[i];
    }

  }

  u=u+1;
  for (i=0; i<n; i++) g[i]=u*x[i];

  free(a);
  free(x);
  free(aa);
  free(xx);
}
