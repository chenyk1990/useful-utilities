/*********************************************************************
*	decon.c:
*		deconvolution of 3 components with estimated source
*		function using either freq. domain water-level
*		deconvolution or time domain Wiener filtering.
*
*	Author:  Lupei Zhu
*
*	Revision History
*		09/04/2000	Initial coding modified from water_src.c
*		11/08/2000	add cut option
*		11/10/2000	add filtering option
*		12/20/2002	combine with water_src.c
*		06/08/2003	add -T and modify -C to include other tmarks
*********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Complex.h"
#include "sac.h"

int main(int argc, char **argv) {
  int 		i, j, k, nn, nft, nft2, nshift;
  int		cut, hpass, norm, nCom, error, m, tref;
  char		outf[128];
  float         filterLen, c, gauss, t1, t2, hf1, hf2, dt, shift, taper;
  float		*src, *data, *acr, *fpt, *wndw, *trace;
  SACHEAD	hd;
  void		toeplitz(const float *, float *, int, float);
  void		water(complex *,int,float,float,float);
  
  cut = 0;
  hpass = 0;
  norm = 0;
  nCom = 3;
  c = 0.01;		/* water-level/white-noise-level in % */
  gauss = 5.;		/* Gaussian low-pass filter */
  filterLen = 0.;
  shift = 5.;
  error = 0;
  taper = 0.;
  src = NULL;
  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {

       case 'C':
	 cut = 1;
	 sscanf(&argv[i][2],"%d/%f/%f",&tref,&t1,&t2);
         break;
  
       case 'E':
         sscanf(&argv[i][2],"%f",&shift);
	 assert(shift>0.);
         break;
  
       case 'F':
	 sscanf(&argv[i][2],"%f",&filterLen);
	 assert(filterLen>0.);
         break;
  
       case 'G':
         sscanf(&argv[i][2],"%f",&gauss);
	 assert(gauss>0.);
         break;
  
       case 'H':
	 hpass = 1;
	 sscanf(&argv[i][2],"%f/%f",&hf1,&hf2);
	 break;

       case 'S':
	 strcpy(outf,&argv[i][2]);
	 break;
  
       case 'T':
         sscanf(&argv[i][2],"%f",&taper);
	 break;

       case 'W':
         sscanf(&argv[i][2],"%f",&c);
         break;

       default:
         error = TRUE;
         break;
  
       }
    } 
  }

  if (argc == 1 || error) {
    fprintf(stderr,"usage: %s -Ssrc [(-Flen|-Ggauss)] [-Cmark/t1/t2 -Eshift -Hf1/f2 \
-I -MmaxShift -N -Ttaper -Wc] data_filenames ...\n\
	-S: specify the name of the source-time or Green function\n\
	-F: use the Wiener filtering with length len (no)\n\
	-G: use water-level decon. with the Gaussian param. (5)\n\
	-C: window data [tmark+t1,tmark+t2] (off)\n\
		mark = -5(b), -3(o), -2(a), 0-9 (t0-t9)\n\
	-E: leave shift sec before zero time (5 sec)\n\
	-H: high-pass filter data (off)\n\
	-T: taper the trace with cosine window, taper=0-0.5 (0.)\n\
	-W: water-level or white noise level (0.01)\n",argv[0]);
    return -1;
  }

  if (cut) {
     src=read_sac2(outf, &hd, tref, t1, t2);
  } else {
     src=read_sac(outf, &hd);
     tref = -5;		/* begin of the trace */
     t1 = 0.;
     t2 = hd.npts*hd.delta;
  }
  nn = hd.npts;
  if (src==NULL) {
    fprintf(stderr,"fail to get data from %s\n",outf);
    return -1;
  }
  dt = hd.delta;
  if (taper>0.) {
     wndw = coswndw(nn,taper);
     for(j=0;j<nn;j++) src[j] *= wndw[j];
  }

  m = rint(filterLen/dt); if (m<1 || m>nn) m=nn;
  nshift = rint(shift/dt);
  nft = 1; while (nft < nn) nft *= 2; nft2=nft; nft *= 2;

  src = realloc(src, nft*sizeof(float));
  for(j=nn;j<nft;j++) src[j]=0.;
  fftr((complex *)src, nft2, dt);
  if (hpass) filter((complex *)src,nft2,hf1,hf2,dt,1);

  if (filterLen>0.) {	/* wiener filtering */
     if ( (acr = (float *) malloc(nft*sizeof(float))) == NULL) {
       fprintf(stderr,"fail to allocate memory for auto-corr\n");
       return -1;
     }
     memcpy(acr, src, nft*sizeof(float));
     cor((complex *)src, (complex *)acr, dt, nft2);
     acr += nft2;
  } else {	/* water-level deconvolution */
     water((complex *)src,nft2,dt,c,gauss);
  }
  
  k = 0;
  data = malloc(nft*sizeof(float));
  while ( ++k<argc ) {
    
    if (argv[k][0] == '-') continue;
    /* data input */
    fprintf(stderr,"%s\n",argv[k]);
    strcpy(outf,argv[k]);
    if ( (trace=read_sac2(outf, &hd, tref, t1, t2))==NULL ||
	  hd.npts != nn ) continue;
    if (taper>0.) {
	for(j=0;j<nn;j++) data[j] = trace[j]*wndw[j];
    } else {
    	for(j=0;j<nn;j++) data[j] = trace[j];
    }
    for(j=nn;j<nft;j++) data[j]=0.;
    free(trace);
    fftr((complex *)data, nft2, dt);
    if (nshift>0) shiftSpec((complex *)data, nft2, nshift);
    if (hpass) filter((complex *)data,nft2,hf1,hf2,dt,1);

    if (filterLen>0.) {
       cor((complex *)src, (complex *)data, dt, nft2);
       fpt = data+nft2;
       toeplitz(acr, fpt, m, c);
    } else {
       specMul((complex *)data, (complex *)src, nft2);
       fftr((complex *)data, nft2, -dt);
       fpt = data;
    }
  
    /* output results */
    strcat(outf,"d");
    hd.npts = m;
    hd.b = *((float *)&hd + 10 + tref) - shift;
    hd.e = hd.b + hd.npts*hd.delta;
    hd.user1 = gauss;
    write_sac(outf,hd,fpt);
  }

  return 0;
  
}
