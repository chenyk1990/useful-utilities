#ifndef __MY_INV__
  #define __MY_INV__

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*3D grid*/
typedef struct {
        int     n[3];           /* number of points */
        float   x0[3], step[3]; /* min and step */
        float   *err;
} GRID;

extern int	svdrs(float *, int, int, int, float *, int, float *);
extern float	iter(float *, float *, float *(*f)(float *), int, int, int, int, float, int);
extern float	marquardt(float *, float *, float*, float *(*f)(float *), int, int, int, int, float, int);
float		grid2d(float *,int,int,float *,float *,float *,float *,float *,int *,int *);
float		grid3d(float *,int *,float *, int *, int *, int *);

#endif
