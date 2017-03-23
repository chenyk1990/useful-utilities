/*************************************************************
*	Name:		grdmin.c:
*
*	Description:	grid-search for minimum of a GMT grd file
*
*	Author:		Lupei Zhu, Dec., 1999 at USC
**************************************************************/
/*#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>*/
#include "gmt.h"
#include "inversion.h"

int main(int argc, char **argv) {
  int	i,j,n,fg;
  char	temp[128];
  GMT_LONG	error = FALSE, deviation = FALSE;
  float	x,vx,y,vy,vxy,p,dev,*a;
  struct GRD_HEADER grd;

  argc = GMT_begin(argc,argv);

  GMT_grd_init(&grd, argc, argv, FALSE);

  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
       case 'D':
         deviation = TRUE;
	 break;
       default:
         error = TRUE;
         break;
       }
    }
  }

  if (argc == 1 || error) {
     fprintf(stderr,"usage: %s [-D] grd_file_list\n\
	-D: output standard deviation at the min using grdfile.std\n\
	Output: best_x best_y P2S/P2x P2S/P2y P2S/PxPy min [deviation]\n",argv[0]);
     return -1;
  }

  i = 0;
  while ( ++i<argc )  {

      if (argv[i][0] == '-') continue;

      fprintf (stderr, "%s\n",argv[i]);
      if (GMT_read_grd_info (argv[i], &grd)) {
	  fprintf (stderr, "grdinfo: Error opening file %s\n", argv[i]);
	  continue;
      }
      n = grd.nx*grd.ny;

      a = (float *) malloc(sizeof(float)*n);
      if ( a==NULL ||
              GMT_read_grd(argv[i],&grd,a,0.,0.,0.,0.,GMT_pad,FALSE)
	      ) continue;

      /* do a grid search to find minimum */
      p = grid2d(a, grd.nx, grd.ny, &x, &vx, &y, &vy, &vxy, &j, &fg);
      vx  /= (grd.x_inc*grd.x_inc);
      vy  /= (grd.y_inc*grd.y_inc);
      vxy = -vxy/(grd.x_inc*grd.y_inc);

      dev = 1.;
      if ( deviation ) {
         strcpy(temp,argv[i]);
	 strcat(temp,".std");
	 if (GMT_read_grd_info (temp, &grd) ||
		GMT_read_grd(temp,&grd,a,0.,0.,0.,0.,GMT_pad,FALSE)) {
	    fprintf(stderr,"%s i/o error\n",temp);
	 } else dev = a[j];
      }

      printf("%s %8.3f %8.3f %8.3f %8.3f %8.3f %8.2e %8.2e\n",argv[i],
	      grd.x_min+x*grd.x_inc, grd.y_max-y*grd.y_inc,
	      vx,vy,vxy,p,dev);
      free(a);

  }

  GMT_end(argc,argv);

  return 0;

}
