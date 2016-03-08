/*********************************************************************
*	iter_decon.c:
*	  Time-domain deconvolution based on Kikuchi and Kanamori
*		(1981) iterative deconvolution algorithm
*
*	Author:  Lupei Zhu
*
*	Revision History
*		01/24/2002	Initial coding modified from based
*				on Ammon's saciterd.f
*		06/08/2003	add -T option and modify -C option
*				to include using other time marks
*		04/12/2007	add positive (-P) and duration (-D)
*				constraint options
*		04/15/2007	add -I option to integrate data
*		04/25/2007	add option to input file names and
*				arrival times from stdin
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Complex.h"
#include "sac.h"
/*#include "buttbp.h"*/

int main(int argc, char **argv) {
  int 		verbose,i,j,k,nn,shift,nft,nft2,cut,positive,intg,diff,n_files;
  int		error,maxiter,saveRes,size,filter,tref,ntau,keepTime;
  char		line[512],tttt[128],gggg[128];
  float		*src,*res,*spikes,*wndw;
  float		amp,amp0,t1,t2,srcAutoCor,dt,pre,thr,rms0,rms,taper,tau,arr,delay;
  complex	*src_spec,*res_spec,*temp;
  SACHEAD	hd;
  int		zp=0, order;
  float		f1, f2, sn[30], sd[30];
  
  error = 0;
  maxiter=50;
  saveRes=0;
  filter=0;
  cut=0;
  pre=0.;
  positive=0;
  thr=0.0001;
  tau=-1.;
  taper=0.;
  intg=0;
  diff=0;
  keepTime=0;
  delay=0.;
  verbose=0;
  n_files=0;
  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
  
       case 'C':
	 cut = 1;
	 sscanf(&argv[i][2],"%d/%f/%f",&tref,&t1,&t2);
         break;
  
       case 'D':
	 sscanf(&argv[i][2],"%f",&tau);
         break;

       case 'F':
	 j = sscanf(&argv[i][2], "%d/%f/%f/%d",&filter,&f1,&f2,&order);
	 if (j<2||filter<0||filter>3||f1<0.) error=TRUE;
	 if (j<3) f2 = -1.;
	 if (j<4) order = 4;
	 if (f2<0.) pre=-f2;
	 break;

       case 'I':
         intg = 1;
	 diff = 0;
	 break;

       case 'J':
         diff = 1;
	 intg = 0;
	 break;

       case 'K':
         keepTime = 1;
	 j = sscanf(&argv[i][2],"%f",&dt);
	 if (j==1) delay = dt;
	 break;

       case 'N':
         sscanf(&argv[i][2],"%d",&maxiter);
         break;
  
       case 'P':
         positive = 1;
	 break;

       case 'T':
         sscanf(&argv[i][2],"%f",&taper);
         break;
  
       case 'V':
	 verbose=1;
	 break;

       case 'W':
	 saveRes=1;
	 break;

       default:
         error = TRUE;
         break;
  
       }
    } else 
       n_files++;
  }

  if (argc == 1 || error) {
    fprintf(stderr,"usage: %s [-Cmark/t1/t2 -Dtau -Fn/f1[/f2[/order]] (-I|-J) -K[delay] -Nnitr -P -Ttaper -V -W src data_filenames ...]\n\
    The program does iterative time-domain deconvolution. Names of the source time function (or Greens function) and data traces\n\
    can be specified in the command line or input from the stdin in the form of:\n\
            sac_file_name [arrival_time].\n\
    other options are:\n\
	-Cmark/t1/t2: windowing data with [tmark+t1,tmark+t2] (off)\n\
		mark = -5(b), -3(o), -2(a), 0-9 (t0-t9)\n\
	-Dtau	add duration constraint (off)\n\
	-F	apply a filter to output (n=1), src (n=2), and data (n=3) (0)\n\
		f1, f2: parameter of the Gaussian low-pass filter and\n\
			the pre-signal length in sec if f2 is negative or\n\
		        the corner freqs. of the Butterworth filter if f2 is positive\n\
		order: the order of the butterworth filter (4, must be < 10).\n\
	-I	integrate the data before deconvolution. The output will be diffed.\n\
	-J	differentiate the source time function before deconvolution. The output will be diffed\n\
	-K	keep the same time as in the data (plus a delay, default is to set the onset time to 0.)\n\
	-Nnitr  number of iterations (50)\n\
	-P:	add positive constraint (off)\n\
	-Ttaper taper the trace with cosine window, taper=0-0.5 (0.)\n\
	-V:	verbose\n\
	-W:	output residual trace (off)\n\
    Example\n\
	iter_decon -SKUL.z -G3 -N100 -C-2/-10/80 -T0.1 KUL.[r,t]\n",
    argv[0]);
    return -1;
  }

  k = 0;
  src = NULL;
  while( (n_files && ++k<argc) || (n_files==0 && fgets(line,512,stdin))) {

    
    /* get sac file name */
    if (n_files) {      /* from command line */
       if (argv[k][0] == '-') continue;
       else strcpy(tttt,argv[k]);
       arr = 0.;
    } else {		/* from stdin */
       i=sscanf(line, "%s %f",tttt,&arr);
       if (i>1) tref = -1;
       else arr = 0.;
    }
    fprintf(stderr,"%s\n",tttt);

    if (src==NULL) { 	/* input denominator */
       if (cut) {
          src=read_sac2(tttt, &hd, tref, arr+t1, arr+t2);
       } else {
          src=read_sac(tttt, &hd);
          tref=-5;
          t1=0.;
          t2=hd.npts*hd.delta;
       }
       if (src==NULL) {
         fprintf(stderr,"fail to get data from %s\n",tttt);
         return -1;
       }
       nn = hd.npts;
       dt = hd.delta;
       if (diff) diffrt(src,nn,dt);
/*       if (filter&&f2>0.) ButtDesign(order,dt,f1,f2,sn,sd);*/
       if (taper>0.) {
          wndw = coswndw(nn,taper);
          for(j=0;j<nn;j++) src[j] *= wndw[j];
       }

       /* memory allocation */
       nft = 1; while (nft < nn) nft *= 2; nft2 = nft; nft *= 2;
       ntau = nft2;
       if (tau>0.) {ntau = ceil(tau/dt); if (ntau>nft2) ntau = nft2;}
       size = nft*sizeof(float);
       spikes =     (float *) malloc(size);
       src_spec = (complex *) malloc(size);
       res_spec = (complex *) malloc(size);
       temp = (complex *) malloc(size);

       /* compute spectrum of the denominator */
       for(j=0;j<nft2;j++) src_spec[j]=Zero;
       memcpy(src_spec, src, nn*sizeof(float));
/*       if (filter>1&&f2>0.) apply_((float *)src_spec,&nft,&zp,sn,sd,&order);*/
       fftr(src_spec, nft2, dt);
       if (filter>1&&f2<0.) fltGauss(src_spec, nft2, f1*dt);
       srcAutoCor = specPwr(src_spec,nft2)/dt;

    } else { 		/* data input */
       if ( (res=read_sac2(tttt, &hd, tref, arr+t1, arr+t2))==NULL ||
	    hd.npts != nn ) continue;
       if (taper) for(j=0;j<nn;j++) res[j] = res[j]*wndw[j];
       if (intg) cumsum(res,nn,dt);
       for(j=0;j<nft2;j++) res_spec[j]=Zero;
       memcpy(res_spec, res, nn*sizeof(float));
/*       if (filter>2&&f2>0.) apply_((float *)res_spec,&nft,&zp,sn,sd,&order);*/
       fftr(res_spec,nft2,dt);
       if (filter>2&&f2<0.) fltGauss(res_spec, nft2, f1*dt);
       rms0=specPwr(res_spec,nft2);

       /* iterative deconvolution */
       for(j=0;j<nft;j++) spikes[j]=0.;
       i = 0; amp = 1.; amp0 = 1.;
       while(i++<maxiter && fabs(amp/amp0)>thr) {
	   memcpy(temp, res_spec, size);
           cor(src_spec, temp, dt, nft2);
	   if (positive) shift = findMax((float *)temp+nft2, ntau, &amp);
	   else shift = findMaxAbs((float *)temp+nft2, ntau, &amp);
	   amp = amp/srcAutoCor;
	   if (i==1) amp0 = amp;
	   memcpy(temp, src_spec, size);
	   shiftSpec(temp, nft2, (float) shift);
	   specScale(temp, -amp, nft2);
	   specAdd(res_spec, temp, nft2);
	   rms = sqrt(specPwr(res_spec, nft2)/rms0);
	   spikes[shift] += amp/dt;
	   if (verbose)
	       fprintf(stderr,"iteration %d amp=%6.2f shift=%6.2f rms=%4.2f\n",
			i, amp, shift*dt, rms);
       }
    
       /* filter the spikes and output results */
       if (filter) {
          if (f2>0.) {
/*             apply_(spikes,&nft,&zp,sn,sd,&order);*/
          } else {
             fftr((complex *)spikes, nft2, dt);
	     shiftSpec((complex *)spikes, nft2, pre/dt);
	     fltGauss((complex *)spikes, nft2, f1*dt);
	     fftr((complex *)spikes, nft2, -dt);
          }
       }
       strcpy(gggg,tttt);
       strcat(gggg,"i");
       hd.npts=nn; if (tau>0.&&!filter) hd.npts=ntau;
       hd.b = -pre;
       if (keepTime) {
          if (tref==-1) hd.b += arr+delay;
          else hd.b += *((float *)&hd + 10 + tref)+delay;
       }
       hd.e = hd.b + hd.npts*hd.delta;
       hd.user1 = f1;
       if (intg||diff) diffrt(spikes,nft2,dt);
       write_sac(gggg, hd, spikes);
       if (saveRes) {
          fftr(res_spec,nft2,-dt);
          strcat(gggg,".res");
          write_sac(gggg,hd,(float *)res_spec);
       }
       free(res);
    }
  }

  return 0;
  
}
