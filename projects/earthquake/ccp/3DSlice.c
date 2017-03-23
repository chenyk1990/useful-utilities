/*************************************************************
*	Name:		3DSlice.c:
*
*	Description:	make vertical slice of a 3D grid
*
*	Note:		The binary 3D file (x fastest, z slowest)
*
*	Author:		Lupei Zhu, March, 1999 at USC
*		7/26/99	initial coding
*		12/3/00 add interpolation
**************************************************************/
/*#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>*/
#include "gmt.h"
#include "grd3d.h"

int main(int argc, char **argv) {
  int		i,j,k,ix,iy,n,error=0,ordr;
  char		outf[128];
  float		*a3,*ss,coef[4];
  float		con,deg,lam,phi,az,cosAz,sinAz,x0,y0,seg,x,y;
  FILE		*fd;
  struct GRD_HEADER grd;
  GRD3D		hd;
  void		intp2D(float, float, float *, int);

  GMT_grdio_init();
  GMT_grd_init(&grd, argc, argv, FALSE);

  deg = 3.1415926/180.;
  con = 6371*deg;

  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
  
       case 'R':
         sscanf(&argv[i][2],"%lf/%lf/%lf",&grd.x_min,&grd.x_max,&grd.x_inc);
         break;
  
       case 'C':
	 sscanf(&argv[i][2],"%f/%f",&lam,&phi);
	 break;

       case 'A':
	 sscanf(&argv[i][2],"%f",&az);
	 break;

       case 'G':
         sscanf(&argv[i][2],"%s",outf);
         break;
  
       default:
         error = TRUE;
         break;
  
       }
    } else {
       fd = fopen(argv[i],"rb");
    }
  }

  if (argc == 1 || error || fd == NULL ) {
     fprintf(stderr,"usage: %s 3D_file -Rx1/x2/dx -Goutput_name -Clam/phi -Aazimuth\n",argv[0]);
     return -1;
  }

  fread(&hd,sizeof(GRD3D),1,fd);
  if (hd.n[0]<2 || hd.n[1]<2) {
     fprintf(stderr,"the input is a 2D file\n");
     return -1;
  }
  x = hd.step[0]>hd.step[1] ? hd.step[0] : hd.step[1];
  ordr = log(2*x/grd.x_inc)/log(2.); if (ordr<1) ordr=1;
  fprintf(stderr," interpolation order = %d\n",ordr);
  /* allocate memory and initializing */
  grd.node_offset = 0;	/* node registration */
  grd.nx = rint((grd.x_max-grd.x_min)/grd.x_inc)+1;
  grd.ny = hd.n[2];
  n = grd.nx*grd.ny;
  if ((ss = (float *) malloc(n*sizeof(float))) == NULL ||
      (a3 = (float *) malloc(hd.n[0]*hd.n[1]*hd.n[2]*sizeof(float))) == NULL ) {
     fprintf(stderr, "unable to allocate memory\n");
     return -1;
  }
  for(k=0; k<n; k++) {
     ss[k] = 0.;
  }
  fread(a3,sizeof(float),hd.n[0]*hd.n[1]*hd.n[2],fd);
  fclose(fd);

  az = hd.az-az*deg;
  cosAz = cos(az);
  sinAz = sin(az);
  x0 = (lam-hd.lam)*cos(hd.phi*deg)*con;
  y0 = (phi-hd.phi)*con;
  az = x0*sin(hd.az) + y0*cos(hd.az) - hd.min[0];
  y0 = y0*sin(hd.az) - x0*cos(hd.az) - hd.min[1];
  x0 = az;
  for(seg=grd.x_min,i=0;i<grd.nx;i++,seg+=grd.x_inc) {
      x = (x0 + seg*cosAz)/hd.step[0];
      y = (y0 + seg*sinAz)/hd.step[1];
      ix = floor(x);	iy = floor(y);
      if (ix<0 || ix>hd.n[0]-1 || iy<0 || iy>hd.n[1]-1) {
	for(j=0; j<grd.ny; j++) ss[i+j*grd.nx]=GMT_f_NaN;
	continue;
      }
      x = x-ix;		y = y-iy;
      intp2D(x,y,coef,ordr);
      k = ix + iy*hd.n[0];
      for (j=0; j<grd.ny; j++,k+=hd.n[0]*hd.n[1]) {
	 ss[i+j*grd.nx] = coef[0]*a3[k+hd.n[0]]+coef[1]*a3[k+hd.n[0]+1]
				+coef[2]*a3[k]+coef[3]*a3[k+1];
      }
  }

  /* output result */
  grd.y_min = -hd.max[2];
  grd.y_max = -hd.min[2];
  grd.y_inc =  hd.step[2];
  GMT_write_grd (outf,&grd,ss,0.,0.,0.,0.,GMT_pad,FALSE);

  return 0;

}
