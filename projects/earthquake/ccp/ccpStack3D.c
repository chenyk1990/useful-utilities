/**********************************************************************************************
*	Name:		ccpStack3D.c:
*
*	Description:	common-conversion-point stack receiver functions
*
*	Note:		ray parameter p is stored in SAC header user0
*			velocity model is in the form of
*			  thickness vp vp/vs
*			and can be terminated by setting thickness to 0
*
*	Author:		Lupei Zhu, March, 1999 at Caltech
*		4/30/99	initial coding
*		6/14/99	add amplitude correction wrt p
*		10/6/99 use Nth-root stacking (N=2)
*		02/27/2000 add -X standard deviation calculation
*		01/21/2001 add option -H to include top layer thick. variation
*		02/09/2002 add -F to do un-flattening
*		02/12/2002 add -E to use Fresnel zone size for smoothing
*		07/18/2002 change -H option for input vel of individual stn
*		12/31/2004 change -H option for input vel of individual stn
*		12/25/2009 add -T to use for S receiver functions
*		02/21/2010 change input velocity model format.
*		04/19/2012 user can specify dominant period for computing Fresnel zone size.
*			   use DoneThisRay to teminate a ray that bottoms inside the model.
*		08/29/2012 add -P to allow user the choose amplitude correction or not.
***********************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "sac.h"
#include "grd3d.h"
#include "Complex.h"

#define MAXS	500	/* max # of stations */
#define MAXL	800	/* max # of lays in the velocity model */
#define NRAYS	900	/* max # of rays in one bin */

int main(int argc, char **argv) {
  int		i,j,k,ix,iy,iz,n,layer,nlayers,nlayers0,nFiles,error,map,srf;
  int		iix,iiy,smthX,smthY,smthx,smthy,verberation,multiple,amc;
  int 		*nn, nRoot,std,stn_vel,nst,flat,verbose,nlyr[MAXS];
  char		line[128], fname[128], outf[128], hfile[128], stn[MAXS][16], *cpt;
  SACHEAD	hd;
  float		*r,*ss,*va[NRAYS],*rr;
  float		vp[MAXL],vs[MAXL],dpth[MAXL],fratio[MAXL];
  float		vp0[MAXL],vs0[MAXL],dpth0[MAXL],fratio0[MAXL];
  float		stnDp[MAXS][MAXL],stnVs[MAXS][MAXL],stnVp[MAXS][MAXL],stnFr[MAXS][MAXL];
  float		p,t,tJump,delZ,z,z0,h,x,delX,y,delY,am,temp,smtht,velo;
  float		con,deg,cosPhi0,cosAz,sinAz,cosDiff,sinDiff,am_cor;
  float		unflat_h,flat_min,delz,waveLen,tppps,delZppps,period;
  FILE		*fd;
  GRD3D		grd;
  void deverberation(float *, int, float, float);

  deg = 3.1415926/180.;
  con = 6371*deg;
  smthX=0;
  smthY=0;
  smtht=1.;
  flat = 0;
  period = 0.;
  nRoot=1;
  std = 0;
  stn_vel=0;
  nFiles=0;
  error=0;
  verbose=0;
  map = 0;
  verberation = 0;
  multiple = 0;
  srf = 0;
  amc = 0;

  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
  
       case 'A':
	 sscanf(&argv[i][2],"%f",&grd.az);
	 grd.az *= deg;
	 cosAz = con*cos(grd.az);
	 sinAz = con*sin(grd.az);
	 break;
       case 'C':
	 sscanf(&argv[i][2],"%f/%f",&grd.lam,&grd.phi);
	 cosPhi0 = cos(grd.phi*deg);
	 break;
       case 'D':	/*de-verberation*/
         verberation = 1;
	 break;
       case 'E':
	 sscanf(&argv[i][2],"%f",&period);
	 break;
       case 'F':
	 flat = 1;
	 break;
       case 'G':
         strcpy(outf,&argv[i][2]);
         break;
       case 'H':
         stn_vel = 1;
	 strcpy(hfile,&argv[i][2]);
	 break;
       case 'I':
	 sscanf(&argv[i][2],"%f/%f/%f",&grd.step[0],&grd.step[1],&grd.step[2]);
	 break;
       case 'M':
         map = 1;
	 break;
       case 'N':
	 sscanf(&argv[i][2],"%d",&nRoot);
	 break;
       case 'P':
	 amc = 1.;
	 break;
       case 'Q':
         verbose = 1;
	 break;
       case 'R':
         sscanf(&argv[i][2],"%f/%f/%f/%f/%f/%f",&grd.min[0],&grd.max[0],
						&grd.min[1],&grd.max[1],
						&grd.min[2],&grd.max[2]);
         break;
       case 'S':
	 j=sscanf(&argv[i][2],"%d/%d/%f",&smthX,&smthY,&temp);
	 if (j==3) smtht=temp;
	 break;
       case 'T':
         srf = 1;
	 break;
       case 'V':
 	 fd = fopen(&argv[i][2],"r");
         break;
       case 'X':
	 std = 1;
	 break;
       case 'Y':
	 multiple = 1;
	 break;
       default:
         error = 1;
         break;
       }
    } else {
       nFiles++;
    }
  }
  
  printf("Hello\n");

  if (argc == 1 || error || fd == NULL ) {
     fprintf(stderr,"usage: %s -Rx1/x2/y1/y2/z1/z2 -Idx/dy/dz -Vvelocity_model -Goutput_name -Clongitude/latitude -Aazimuth [-D -Es -F -Hfile -Nn -P -Q -SsmthX/smthY/[smtht] -T -X -Y RF_SAC_files] \n\
   SAC files should have stla,stlo,stel,baz,user0 (p),kstnm set in head.\n\
   If not sac file is given in the command line, it will read in from stdin.\n\
	-A:  profile azimuth \n\
	-C:  profile origin \n\
	-D:  de-reverberation (off) \n\
	-E:  use s as the dominant period to compute the Fresnel zone size for the stacking area (0.) \n\
	-F:  do earth un-flattening (off)\n\
	-G:  output file name \n\
	-H:  velocity models for individual stations. File format: \n\
		stationName n \n\
		thickness vp vp/vs\n\
		... \n\
		stationName n \n\
		... \n\
	     Note: use with -F has not be tested yet\n\
        -I:  pixel size, x is along the profile; y is ccw; z is down. \n\
	-M:  map RF traces to depth domain (off) \n\
	-N:  use n-th root stacking (linear) \n\
	-P:  do amplitude incident angle correction using an empirical (no).\n\
	-Q:  verbose \n\
	-R:  stacking volume \n\
	-S:  smooth lengths (in number of pixels/samples) (0/0/1) \n\
	-T:  use S receiver functions (0).\n\
	-V:  velocity model (thickness vp vp/vs) \n\
	-X:  compute standard deviation of stack (off, use number of rays) \n\
	-Y:  use PpPs phase instead of Ps\n",argv[0]);
     return -1;
  }

  /* input velocity model: layer_thickness, vp, vp/vs
     each layer is represented by 
	dpth: depth of the bottom
	vp:   1/vp^2
	vs:   1/vs^2
	fratio: un-flatten ratio, fratio*thickness = true thickness
 */
  layer = 0;
  h = 0.;
  unflat_h = 0.;
  flat_min = 0.; /* this is z_min after flattening */
  while ( fscanf(fd, "%f%f%f",&dpth[layer],&vp[layer],&vs[layer]) == 3 ) {
     printf("dpth[%d]=%g,vp[%d]=%g,vs[%d]=%g\n",layer,dpth[layer],layer,vp[layer],layer,vs[layer]);
     vp[layer] = 1./(vp[layer]*vp[layer]);
     vs[layer] = vp[layer]*(vs[layer]*vs[layer]);
     h += dpth[layer];
     fratio[layer] = 1.;
     if (flat) fratio[layer] = 1.-(unflat_h+0.5*dpth[layer])/6371.;
     unflat_h += dpth[layer]*fratio[layer];
     if (grd.min[2]>0. && flat_min==0. && unflat_h>grd.min[2]) {
	flat_min = unflat_h-(unflat_h-grd.min[2])/fratio[layer];
     }
     if (dpth[layer]<0.001*grd.step[2] || unflat_h>grd.max[2]) {
	layer++;
	break;
     }
     dpth[layer] = h;
     layer++;
  }
  
  
  
    printf("Hello2, h=%g\n",h);
  fclose(fd);
  fprintf(stderr,"model is terminated at depth %f %f %f\n",h,flat_min,unflat_h);
  dpth[layer-1] = 2*grd.max[2]; /* half space */
  nlayers = layer;

  /* input optional velocity model for individual stations */
  if (stn_vel) {
     /* first we need to save the master velocity model */
     nlayers0 = nlayers;
     for(layer=0;layer<nlayers;layer++) {
	dpth0[layer] = dpth[layer];
        vp0[layer] = vp[layer];
	vs0[layer] = vs[layer];
	fratio0[layer] = fratio[layer];
     }
     /* read in individual models */
     fd = fopen(hfile,"r");
     assert(fd != NULL);
     nst = 0;
     while (fscanf(fd, "%s %d",stn[nst],&nlyr[nst])==2) {
	for(h=0.,layer=0;layer<nlyr[nst];layer++) {
	   fscanf(fd,"%f%f%f",&z,&x,&y);
	   h += z;
	   stnDp[nst][layer] = h;
	   stnVp[nst][layer] = 1./x/x;
	   stnVs[nst][layer] = stnVp[nst][layer]*y*y;
	   stnFr[nst][layer] = 1.;
	}
        nst++;
     }
     fclose(fd);
  }
printf("Hello3\n");
  /* allocate memory and initializing */
  grd.n[0] = rint((grd.max[0]-grd.min[0])/grd.step[0])+1;
  grd.n[1] = rint((grd.max[1]-grd.min[1])/grd.step[1])+1;
  grd.n[2] = rint((grd.max[2]-grd.min[2])/grd.step[2])+1;
  n = grd.n[0]*grd.n[1]*grd.n[2];
  if ((ss = (float *) malloc(n*sizeof(float))) == NULL ||
      (va[0] = (float *) malloc(n*sizeof(float))) == NULL ||
      (rr = (float *) malloc(grd.n[2]*sizeof(float))) == NULL ||
      (nn =   (int *) malloc(n*sizeof(int))) == NULL ) {
     fprintf(stderr, "unable to allocate memory\n");
     return -1;
  }
  if (std) {
     for(i=1;i<NRAYS;i++) {
         if ((va[i] = (float *) malloc(n*sizeof(float))) == NULL) {
             fprintf(stderr,"unable to allocal memory for STD comp.\n");
             return -1;
         }
     }
  }
  for(k=0; k<n; k++) {
     ss[k] = 0.;
     va[0][k] = 0.;
     nn[k] = 0;
  }
printf("Hello4\n");
  /** input traces one by one and do stacking **/
  i=0;
  while ( (nFiles && ++i<argc) || (!nFiles && fgets(line,128,stdin)) ) {
      if (nFiles) {
	 if (argv[i][0] == '-') continue;
	 strcpy(fname, argv[i]);
      } else {
         sscanf(line,"%s",fname);
      }
      if ( (r=read_sac(fname,&hd)) == NULL) continue;
      if (stn_vel) {
	 layer = 0;
	 h = 0.;
	 cpt = strchr(hd.kstnm,' ');
	 if (cpt != NULL) {k = cpt-hd.kstnm; hd.kstnm[k] = '\0';} else k=strlen(hd.kstnm);
	 for (j=0;j<nst;j++) {
	    if ( k==strlen(stn[j]) && strncmp(hd.kstnm, stn[j], k)==0 ) break;
	 }
	 if (j<nst) {
	    for(;layer<nlyr[j];layer++) {
	       dpth[layer] = h = stnDp[j][layer];
	       vs[layer] = stnVs[j][layer];
	       vp[layer] = stnVp[j][layer];
	       fratio[layer] = stnFr[j][layer];
	    }
	 } else {
	    fprintf(stderr,"Warning: %s not found in the input velocity model\n",hd.kstnm);
	 }
         for(j=0;j<nlayers0;j++) {
	    if (dpth0[j]>h) {
	       dpth[layer] = dpth0[j];
               vp[layer] = vp0[j];
	       vs[layer] = vs0[j];
	       fratio[layer] = fratio0[j];
	       layer++;
	    }
         }
	 nlayers = layer;
	 if (verbose) {
	    for(layer=0;layer<nlayers;layer++) {
	       fprintf(stderr,"%s %d %f %f %f %f\n",hd.kstnm,layer,dpth[layer],1./sqrt(vs[layer]),1./sqrt(vp[layer]),fratio[layer]);
	    }
         }
      }
      p = hd.user0*hd.user0;
      am_cor = 1.; if (amc) am_cor = 151.5478*p + 3.2896*sqrt(p) + 0.2618;	/* am corr for 3.0/3.5*/
      /* t, x, y, z: values along the ray Ps, start from the surface */
      t = -hd.b/hd.delta;
      tppps = t;
      z0 = h = -0.001*hd.stel;
      if ( z0>flat_min ) {
	 if (verbose) fprintf(stderr,"%s discarded, station sits below z1\n",fname);
	 free(r);
	 continue;
      }
      cosDiff = cos(grd.az-hd.baz*deg);
      sinDiff = sin(grd.az-hd.baz*deg);
      x = (hd.stlo-grd.lam)*cosPhi0*sinAz+(hd.stla-grd.phi)*cosAz;
      y =-(hd.stlo-grd.lam)*cosPhi0*cosAz+(hd.stla-grd.phi)*sinAz;
      /* next point */
      z = flat_min;
      iz = 0;
      layer = -1; delZ=delX=delY=delZppps=1.;	/* the air layer above the surface */
      while(iz<grd.n[2]) {
	  if (z>h) { /* crossing a layer interface */
	      /* finish what's left between z0 and the intf right above z */
	      while (z>h) {
	         tJump = (h-z0)/delZ;
	         t += tJump;
	         x += tJump*delX;
	         y += tJump*delY;
		 tppps += (h-z0)/delZppps;
	         if (verbose) fprintf(stderr,"%3d %6.2f %6.2f %6.2f %6.2f\n",
	  			layer,t,h,x,y);
	         layer++;
	         /* del_XYZ are the changes to x,y,z for t incr by dt */
		 temp = vp[layer]-p;
		 if (temp<0.) goto DoneThisRay; 	/* the ray bottoms out, teminate it */
		 temp = sqrt(temp);
	         delZ     = hd.delta/(sqrt(vs[layer]-p)-temp);
		 delZppps = hd.delta/(sqrt(vs[layer]-p)+temp);
		 velo = vs[layer]; if (srf) velo = vp[layer];
	         temp = delZ/sqrt(velo/p-1.);
	         delX = temp*cosDiff;
	         delY = temp*sinDiff;
	         z0 = h;
	         h = dpth[layer];
		 waveLen = period/sqrt(velo);
		 //fprintf(stderr,"dominant wavelength is %f km\n",waveLen);
	      }
	      /* delz is for 1-grid true depth change  */
	      delz = grd.step[2]/fratio[layer];
          }
	  tJump = (z-z0)/delZ;
	  t += tJump;
	  x += tJump*delX;
	  y += tJump*delY;
	  tppps += (z-z0)/delZppps;
	  /* get the amplitude at t */
	  /*j = floor(t); if (j>=hd.npts-1) break; am = r[j]+(t-j)*(r[j+1]-r[j]);*/
	  if (multiple) am = amp(tppps-smtht,tppps+smtht,r,hd.npts)/(2*smtht);
	  else am = amp(t-smtht,t+smtht,r,hd.npts)/(2*smtht);
	  if (verberation) deverberation(r,hd.npts,tppps,am);
	  am = am/am_cor;
	  if (verbose) fprintf(stderr,"%3d %6.2f %6.2f %6.2f %6.2f %6.2f %f\n",
	  			layer,t,z,x,y,tppps,am);
	  if (map) rr[iz] = am;
	  ix = rint((x-grd.min[0])/grd.step[0]);
	  iy = rint((y-grd.min[1])/grd.step[1]);
	  smthx=smthX;
	  smthy=smthY;
	  temp = 0.5*sqrt(z*waveLen); /*half size of the Fresnel zone*/
	  smthx+=rint(temp/grd.step[0]);
	  smthy+=rint(temp/grd.step[1]);
	  for(iiy=iy-smthy;iiy<=iy+smthy;iiy++) {
	     if (iiy<0 || iiy>=grd.n[1] ) continue;
	     for(iix=ix-smthx;iix<=ix+smthx;iix++) {
	        if (iix<0 || iix>=grd.n[0] ) continue;
		k = iix+grd.n[0]*(iiy+grd.n[1]*iz);
	        ss[k] += ( am>0. ? 1 : -1) * expf(logf(fabs(am))/nRoot);
		if (std) {
		   if (nn[k] >= NRAYS) {
		      fprintf(stderr,"max. # of rays exceeds\n");
		      return -1;
		   }
		   va[nn[k]][k] = am;
		}
	        nn[k] ++;
	     }
	  }
	  z0 = z;
	  iz++;
	  z+=delz;
      }
      DoneThisRay:
      free(r);
      /* write RF in depth domain */
      if (map) {
         hd.delta = grd.step[2];
         hd.b = grd.min[2];
         hd.npts = grd.n[2];
         write_sac(strcat(fname,".map"), hd, rr);
      }
  }
printf("Hello5\n");
  /* output result */
  for (k=0; k<n; k++) {
     if (nn[k]<1) continue;
     ss[k] =  ( ss[k]>0. ? 1 : -1 ) * expf(nRoot*logf(fabs(ss[k]/nn[k])));
     if ( std ) {
        temp = 0.;
        for(i=0;i<nn[k];i++) {
	    temp += (va[i][k]-ss[k])*(va[i][k]-ss[k]);
        }
        if (nn[k]>2) va[0][k] = sqrt(temp/(nn[k]-1));
     } else {
        va[0][k] = nn[k];
     }
  }
  
  for(int j=0;j<n;j++)
  printf("ss[%d]=%g\n",j,ss[j]);
  
  /* write stacking */
  fd = fopen(outf,"wb");
  fwrite(&grd, sizeof(GRD3D), 1, fd);
  fwrite(ss, sizeof(float), n, fd);
  fclose(fd);
  /* write STD  */
  strcat(outf,".std");
  fd = fopen(outf,"wb");
  fwrite(&grd, sizeof(GRD3D), 1, fd);
  fwrite(va[0], sizeof(float), n, fd);
  fclose(fd);

  return 0;

}

/* using 4-point spreading to remove reverberation at time t in r[i] */
void deverberation(float *r, int n, float t, float am) {
  int i;
  float x;
  i = floor(t);
  if (i>0 && i<n-2) {
     x = t-i;
     am = 0.25*am;
     r[i-1] -= am*(1-x);
     r[i]   -= am*(2-x);
     r[i+1] -= am*(1+x);
     r[i+2] -= am*x;
  }
  return;
}
