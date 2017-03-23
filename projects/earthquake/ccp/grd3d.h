/* grd3d.h
   header for storing data (float) on a 3D grids, z pointing down, x is the fastest index
*/
#ifndef _GRD3D_
#define _GRD3D_

typedef struct {
	float lam,phi,az;	/* longitude, latitude of the origin and azimuth (in radians) of the x-axis */
	float min[3],max[3],step[3];	/* in km, in order of x, y and z */
	int   n[3];		/* number of grids in x, y, and z directions */
} GRD3D;

#endif
