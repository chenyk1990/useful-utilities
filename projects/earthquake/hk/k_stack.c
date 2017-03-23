/*************************************************************
*	Name:		k_stack.c:
*
*	Description:	H-k stacking of receiver functions to
*			estimate crustal thickness and Vp/Vs ratio.
*			For refererence, see Zhu and Kanamori, JGR 105, 2000.
*
*	Author:		Lupei Zhu, March, 1997 at Caltech
*
* modification history:
*	10/6/01 LZ	add ability to search k with a fixed tps
*			by inputing the H and k as xmin/xmax
*	04/24/04 LZ	add -M to add more multiples (basin reverb).
*	12/30/04 LZ	use tps as input for fixed tps search.
*	05/20/09 LZ	add nth-root stacking option.
**************************************************************/
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
#include "gmt.h"
#include "sac.h"
#include "Complex.h"

#define NS 1000		/* max. number of traces */

int main(int argc, char **argv) {
  int	i,j,k,n,num_of_traces,npts,nFiles=0,fg;
  char	outf[128], line[128], inf[128];
  GMT_LONG	error = FALSE, addmore=FALSE;
  SACHEAD	hd;
  float	*r,*ss, *var, *data, *fpt, *buff[NS];
  float	p, vp=6.3, vs, eta_p, eta_s, t0, dt, w1, w2, w3, temp, am_cor;
  float	tps, dtps, tppps, dtppps, tpsps, dtpsps, smtht, cxy;
  struct GRD_HEADER grd;

  argc = GMT_begin(argc, argv);
  GMT_grd_init(&grd, argc, argv, FALSE);

  w1 = 0.7; w2 = 0.2; w3 = 0.1;
  grd.x_inc = 0.5; grd.y_inc = 0.01;
  smtht = 1.;
  int nRoot = 1;
  int cxWidth = 0;

  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
       case 'G':
         sscanf(&argv[i][2],"%s",outf);
         break;
       case 'I':
         sscanf(&argv[i][2],"%lf/%lf",&grd.x_inc,&grd.y_inc);
         break;
       case 'M':
         addmore = TRUE;
         break;
       case 'N':
         sscanf(&argv[i][2],"%d",&nRoot);
	 break;
       case 'R':
         sscanf(&argv[i][2],"%lf/%lf/%lf/%lf",&grd.x_min,&grd.x_max,&grd.y_min,&grd.y_max);
         break;
       case 'T':
         sscanf(&argv[i][2], "%f",&smtht);
	 break;
       case 'U':
         sscanf(&argv[i][2], "%d",&cxWidth);
	 break;
       case 'V':
         sscanf(&argv[i][2], "%f",&vp);
         break;
       case 'W':
         sscanf(&argv[i][2],"%f/%f/%f",&w1,&w2,&w3);
         break;
       default:
         error = TRUE;
         break;
       }
    } else {
       nFiles++;
    }
  }

  if (argc == 1 || error) {
     fprintf(stderr,"H-k stacking of receiver functions. \n\
usage: %s -Goutputfile -Rz1/z2/k1/k2 [-Idz/dk -M -Nn -Tsmth -Un -Vvp -Ww1/w2/w3] [file_list] \n\
  file_list is a list of SAC files of receiver functions (from stdin if not provided). The ray parameter is stored in the user0 field in the SAC header.\n\
  -R: crustal thickness Vp/Vs ranges \n\
     if k1=k2, output stack(tps) in SAC format for fixed k, tps is for p=0.06 s/km\n\
     if z1=z2, output stack(k) in SAC format for fixed tps=z1\n\
     else, output stack(H,k) in GMT grid format\n\
  -G: output file name \n\
  -I: grid spacing \n\
  -M: add more mutiples (off) \n\
  -N: use n-th root stacking (1) \n\
  -T: smooth in the time domain (1 sample) \n\
  -U: add using cc coefficents of phase Ps and PpPs and specify the half window length (off) \n\
  -V: crustal Vp (6.3 km/s) \n\
  -W: weights for Ps, PpPs and PsPs+PpSs (0.7/0.2/0.1)\n",argv[0]);
     return -1;
  }

  grd.node_offset = 0;
  grd.nx = rint((grd.x_max-grd.x_min)/grd.x_inc)+1;
  grd.ny = rint((grd.y_max-grd.y_min)/grd.y_inc)+1;
  temp = 2*smtht;
  w1 /= temp; w2 /= temp; w3 /= temp;

  n = grd.nx*grd.ny;
  if (n<1 || (ss = (float *) malloc(n*sizeof(float))) == NULL ||
  	    (var = (float *) malloc(n*sizeof(float))) == NULL ) {
     fprintf(stderr, "unable to allocate memory for size %d\n",n);
     return -1;
  }

  /** input traces one by one and do stacking **/
  vp = 1./(vp*vp);
  num_of_traces = 0;
  i = 0;
  while ( (nFiles && ++i<argc) || (!nFiles && fgets(line,128,stdin)) )  {
      if (nFiles) {
	 if (argv[i][0] == '-') continue;
	 strcpy(inf, argv[i]);
      } else
         sscanf(line,"%s",inf);
      if ( (r=read_sac(inf,&hd)) == NULL) continue;
      k= rint(-hd.b/hd.delta);
      if (k < 0) {free(r);continue;}
      data = r + k;
      npts = hd.npts - k;
      dt = grd.x_inc/hd.delta;
      t0 = tps = grd.x_min/hd.delta;
      if ( num_of_traces == NS || (buff[num_of_traces] = fpt = (float *) malloc(n*sizeof(float))) == NULL ) {
         fprintf(stderr, "too many traces or unable to allocate memory for variance\n",n);
         return -1;
      }
      p = hd.user0*hd.user0;
      am_cor = 151.5478*p + 3.2896*sqrt(p) + 0.2618;
      eta_p = sqrt(vp-p);
      num_of_traces++;
      for(vs=grd.y_max,j=0; j<grd.ny; j++,vs-=grd.y_inc) {
         eta_s = sqrt(vp*vs*vs-p);
	 if (grd.nx==1) {
	    t0 = tps/(eta_s - eta_p);
	 } else {
            tps = t0*(eta_s - eta_p);
	 }
         dtps   = dt*(eta_s - eta_p);
         tppps = t0*(eta_s + eta_p);
         dtppps = dt*(eta_s + eta_p);
         tpsps = t0*2.*eta_s;
         dtpsps = dt*2.*eta_s;
         for (k=0; k<grd.nx; k++, fpt++) {
	     *fpt =w1*amp(tps-smtht,tps+smtht, data, npts)/am_cor
	           +w2*amp(tppps-smtht,tppps+smtht, data, npts)
	           -w3*amp(tpsps-smtht,tpsps+smtht, data, npts);
	     if (cxWidth) {
		*fpt *= acc(data, npts, tps, tppps, cxWidth);
	     }
	     if (addmore) {
	        *fpt+=-w1*amp(tps+tpsps-smtht,tps+tpsps+smtht, data, npts)
	           -w2*amp(tppps+tpsps-smtht,tppps+tpsps+smtht, data, npts)
	           +w3*amp(tpsps+tpsps-smtht,tpsps+tpsps+smtht, data, npts)
	           +w1*amp(tps+2*tpsps-smtht,tps+2*tpsps+smtht, data, npts)
	           +w2*amp(tppps+2*tpsps-smtht,tppps+2*tpsps+smtht, data, npts)
	           -w3*amp(tpsps+2*tpsps-smtht,tpsps+2*tpsps+smtht, data, npts);
	     }
	     if ( grd.nx>1 ) tps += dtps;
	     tppps += dtppps;
	     tpsps += dtpsps;
         }
      }
      free(r);
  }

  if (num_of_traces < 1) return -1;

  /* compute mean and standard deviation*/
  for (k=0; k<n; k++) {
     temp = 0.;
     for(i=0;i<num_of_traces;i++) temp += (buff[i][k]>0. ? 1 : -1) * expf(logf(fabs(buff[i][k]))/nRoot);
     ss[k] = (temp>0. ? 1 : -1) * expf(nRoot*logf(fabs(temp)/num_of_traces));
     temp = 0.;
     for(i=0;i<num_of_traces;i++) {
         temp += (buff[i][k]-ss[k])*(buff[i][k]-ss[k]);
     }
     var[k] = sqrt(temp/(num_of_traces-1));
  }

  /* output result as SAC file if we only did H or k search */
  if ( grd.ny == 1 || grd.nx == 1 ) {
     if (grd.ny == 1) {
        hd.npts = grd.nx;
	p=0.06*0.06;
	temp = sqrt(vp*grd.y_max*grd.y_max-p)-sqrt(vp-p);
        hd.b = grd.x_min*temp;
        hd.delta = grd.x_inc*temp;
     } else {
        hd.npts = grd.ny;
        hd.b = grd.y_min;
        hd.delta = grd.y_inc;
	for(j=0,k=n-1;j<k;j++,k--) {
	   temp = ss[j];
	   ss[j] = ss[k];
	   ss[k] = temp;
	   temp = var[j];
	   var[j] = var[k];
	   var[k] = temp;
	}
     }
     write_sac(outf,hd,ss);
     strcat(outf,".std");
     write_sac(outf, hd, var);
  } else {
     for(k=0;k<n;k++) ss[k]=-ss[k];
     GMT_write_grd (outf,&grd,ss,0.,0.,0.,0.,GMT_pad,FALSE);
     strcat(outf,".std");
     GMT_write_grd (outf,&grd,var,0.,0.,0.,0.,GMT_pad,FALSE);
     GMT_end(argc, argv);
  }

  return 0;

}
