/*--------------------------------------------------------
knt.h	structure declaration for knt.c
author	Lupei Zhu
----------------------------------------------------------*/

#ifndef __KNT__
  #define __KNT__

#include "Complex.h"

#define MAXL 500

/******************************************
* reflection/transmission coefficients
*    | pp ps |
*    | sp ss |
*******************************************/
typedef struct {
   complex pp;
   complex ps;
   complex sp;
   complex ss;
   complex sh;
} matrix;

/*******************************************
* use for both eigenvector matrix D
*    |  Mu Md |
*    |  Nu Nd |
* and refl/trans matrices
*    |  Ru Rd |
*    |  Tu Td |
*******************************************/
typedef struct {
   matrix ru;
   matrix rd;
   matrix tu;
   matrix td;
} Dmatrix;

typedef struct {
   float   beta, kapa, thik;	/* Vs, Vp/Vs ratio, thickness */
   float   dbet;	/* perturbed Vs for computing differentials */
   complex qa, qb;	/* vertical P and S slowness */
   complex dqa, dqb;	/* when the layer is perturbed */
   matrix  SRd, STu;	/* Rd and Tu between upper interf and the half space */
   matrix  URu, UTu;	/* Ru and Tu between lower interf and the free surf+ */
   Dmatrix D;		/* eigenvector matrix (p51, M->r,N->t) */
   Dmatrix Q;		/* refl/transm matrix for upper interface */
   Dmatrix dQ1, dQ2;	/* Q for upper/lower interf. when the layer perturbed */
} Layer;


/* declaration of matrix operations */
extern matrix plus(matrix a, matrix b);
extern matrix mltp(matrix a, matrix b);
extern matrix ngtv(matrix a);
extern matrix invs(matrix a);
extern matrix trns(matrix a);
extern matrix imlt(matrix a);

/* declaration of KNT functions */
extern matrix EEE(complex qa, complex qb, float wh);
extern void    ifmat(float p, int nlyrs, Layer *);
extern void   delifm(float p, int nlyrs, Layer *);
extern matrix rcvrfn(float w, int nlyrs, Layer *);
extern matrix rcvrtd(float w, int nlyrs, Layer *);
extern matrix delrcv(float w, Layer *);
extern Layer  *mdSetup(int, const float *, const float *, const float *,
			float db, float dh, int *nlyrs);
extern float  *partial(int, int, int, const float *, const float *, const float *,
			float, float, float, float, float, float);
extern void    respknt(int ps, int, int, const float *, const float *, const float *,
			float, float, complex *, complex *);
extern int    mdin(const char *, float *, float *, float *);

#endif
