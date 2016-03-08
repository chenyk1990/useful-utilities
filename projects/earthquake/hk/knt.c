/*************************************************************************
knt.c: compute response of layered model to teleseismic P

author:	Lupei Zhu

Revision History:
	Lupei Zhu  1991		initial coding from T. J. Ownens's fortran
				version and
				Kennett (1983), Seismic Wave Propagation in
				Stratified  Media, Cambridge U. Press.
	Lupei Zhu  2/24/1999	remove gloabal variables.
	Lupei Zhu  12/21/2009	add S wave incidence in respknt() and partial().

**************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "knt.h"

matrix	I         ={1.,0., 0.,0., 0.,0., 1.,0., 1.,0.};
matrix	ZeroMatrix={0.,0., 0.,0., 0.,0., 0.,0., 0.,0.};


/*****************************************************************************
Set up model from input layer thickness, Vs, and kapa.
Model index:
	layer 0 (air)
    ----------------------------------- free surface
    	layer 1
    ----------------------------------- 
		:
		:
    ----------------------------------- 
   	layer nlyrs, half space

******************************************************************************/
Layer	*mdSetup(int		m,	/* number of input layers */
		 const float	thik[],	/* layer thikness */
		 const float	beta[],	/* Vs */
		 const float	kapa[],	/* Vp/Vs */
		 float		db,	/* amount of Vs perturbation */
		 float		dh,	/* amount of depth perturabtion */
		 int		*nlyrs	/* number of layers in the model */
	 )
{
  int	j;
  Layer	*lyr, *lyr0;
  
  *nlyrs = m;
  if (dh>0.) *nlyrs = 2*m-1;
  if ( (lyr0=(Layer *)malloc((*nlyrs+1)*sizeof(Layer))) == NULL ) return lyr0;
  for(lyr=lyr0+1,j=0; j<m; j++,lyr++) {
     lyr->thik = thik[j] - dh;
     lyr->beta = beta[j];
     lyr->kapa = kapa[j];
     lyr->dbet = beta[j] + db;
     if (dh>0. && j<m-1) {
        lyr++;
        lyr->thik = dh;
        lyr->beta = beta[j];
        lyr->kapa = kapa[j];
        lyr->dbet = beta[j+1];
     }
  }

  return lyr0;

}


/* calculate eigenvector matrix D and qa, qb from layer parameter beta, kapa */
/* see page 46 (3.7) and 50 (3.36) and 51 (3.37 and 3.38)                    */
Dmatrix DDD(float p, float beta, float kapa, complex *qa, complex *qb) {
  Dmatrix D;
  float   rho;
  complex bet, pa, pb, t1, t2, tt;
  rho = 0.77 + 0.32*kapa*beta;
  bet = cmplx(beta, 0.5*beta/225.);
  tt = cmltp(bet, bet);
  t1 = cmplx(-p*p, 0.);
  t2 = cinvs(tt);
  *qb = Csqrt(cplus(t2, t1));
  if (qb->y>0.) *qb = cngtv(*qb);
  t2 = cmplx(kapa*beta, 0.5*kapa*beta/500.);
  t2 = cinvs(t2);
  t2 = cmltp(t2, t2); 
  *qa = Csqrt(cplus(t2, t1));
  if (qa->y>0.) *qa = cngtv(*qa);
  pa = cinvs(Csqrt(dmltp(2.*rho, *qa)));
  pb = cinvs(Csqrt(dmltp(2.*rho, *qb)));
  t1 = dmltp(2.*rho*p, tt);
  t2 = cplus(dmltp(p, t1), cmplx(-rho,0.));
  D.rd.pp = cmltp(IMAGE, cmltp(*qa, pa));	D.ru.pp = cngtv(D.rd.pp);
  D.rd.ps = dmltp(p, pb);			D.ru.ps =       D.rd.ps;
  D.rd.sp = dmltp(p, pa);			D.ru.sp =       D.rd.sp;
  D.rd.ss = cmltp(IMAGE, cmltp(*qb, pb));	D.ru.ss = cngtv(D.rd.ss);
  D.rd.sh = cmltp(cinvs(bet), pb);		D.ru.sh =       D.rd.sh;
  D.td.pp = cmltp(t2, pa);			D.tu.pp =       D.td.pp;
  D.td.ps = cmltp(t1, D.rd.ss);			D.tu.ps = cngtv(D.td.ps);
  D.td.sp = cmltp(t1, D.rd.pp);			D.tu.sp = cngtv(D.td.sp);
  D.td.ss = cmltp(t2, pb);			D.tu.ss =       D.td.ss;
  D.td.sh = cmltp(dmltp(rho,bet), D.rd.ss);	D.tu.sh = cngtv(D.td.sh);
  return(D);
}


/* calculate the refl/trans matrix Q of interface between two layers */
Dmatrix QQQ(Dmatrix D1, Dmatrix D2) {
  Dmatrix Q;
  matrix  T, L;
  T = plus(mltp(trns(D1.ru), D2.td), ngtv(mltp(trns(D1.tu), D2.rd)));
  L = invs(T);
  Q.td = imlt(L);	/* (5.19) on page 105 */
  Q.tu = trns(Q.td);	/* (5.24) on page 106 */
  T = plus(mltp(trns(D1.rd), D2.td), ngtv(mltp(trns(D1.td), D2.rd)));
  Q.rd = ngtv(mltp(T, L));	/* (5.19) on page 105 */
  T = plus(mltp(trns(D1.ru), D2.tu), ngtv(mltp(trns(D1.tu), D2.ru)));
  Q.ru = ngtv(mltp(L, T));	/* (5.24) on page 106 */
  return(Q);
}


/* calculate E matrix through the layer of thickness h at freq. w */
matrix EEE(complex qa, complex qb, float wh) {
  complex iwh;
  matrix  Ed;
  iwh = cmplx(0., -wh);
  Ed.pp = cphase(cmltp(iwh, qa));
  Ed.ss = cphase(cmltp(iwh, qb));
  Ed.ps = Zero;
  Ed.sp = Zero;
  Ed.sh = Ed.ss;
  return(Ed);
}


/* calculate reflection-transmission matrix of each interface in model *lyr */
void	ifmat(float p, int nlyrs, Layer *lyr)
{
  Layer	*halfSpace;
  halfSpace = lyr+nlyrs;
  lyr->D.ru=I; lyr->D.rd=lyr->D.td=lyr->D.tu=ZeroMatrix;	/* air layer */
  while ( lyr++ < halfSpace ) {
      lyr->D = DDD(p, lyr->beta, lyr->kapa, &(lyr->qa), &(lyr->qb));
      lyr->Q = QQQ((lyr-1)->D, lyr->D);
  }
}


/* calculate refl-trans matrix of each interf. when each layer is perturbed. */
void	delifm(float p, int nlyrs, Layer *lyr)
{
  Layer   *halfSpace;
  Dmatrix D;
  halfSpace = lyr+nlyrs;
  while ( (++lyr) < halfSpace ) {
      D = DDD(p, lyr->dbet, lyr->kapa, &(lyr->dqa), &(lyr->dqb));
      lyr->dQ1 = QQQ((lyr-1)->D, D);	/* upper interface */
      lyr->dQ2 = QQQ(D, (lyr+1)->D);	/* lower interface */
  }
}


/**************************************************************************
calculate displacement at free surface using bottom-up recursing.
The intermediate results are stored in STu and SRd in order to speed up the
computation of differential seismogram.
***************************************************************************/
matrix	rcvrfn(float w, int nlyrs, Layer *lyr0)
{
  Layer  *lyr;
  matrix L, Rd, Tu, Ed;

  lyr = lyr0 + nlyrs;	/* recursing begins from the half space */
  lyr->STu = I; lyr->SRd = ZeroMatrix;
  Rd = lyr->Q.rd;
  Tu = lyr->Q.tu;
  while ( (--lyr) > lyr0 ) {
      Ed = EEE(lyr->qa, lyr->qb, w*lyr->thik);
      Rd = mltp(Ed, mltp(Rd, Ed));
      Tu = mltp(Ed, Tu);
      lyr->SRd = Rd; lyr->STu = Tu;	/* save the two */
      L = mltp(lyr->Q.tu, invs(plus(I, ngtv(mltp(Rd, lyr->Q.ru)))));
      Rd = plus(lyr->Q.rd, mltp(lyr->Q.tu, mltp(Rd, trns(L))));
      Tu = mltp(L, Tu);
  }
  return(Tu);
}


/**************************************************************************
calculate displacement at free surface using up-to-bottom recursing.
The intermediate results are stored in UTu and URu in order to speed up the
computation of differential seismogram.
***************************************************************************/
matrix	rcvrtd(float w, int nlyrs, Layer *lyr)
{
  Layer  *halfSpace;
  matrix L, Ru, Tu, Ed;

  halfSpace = lyr+nlyrs;
  lyr->UTu = I; lyr->URu = ZeroMatrix;
  lyr++;	/* recursing begins from the first layer */
  Ru = lyr->Q.ru;
  Tu = lyr->Q.tu;
  while ( lyr < halfSpace ) {
      Ed = EEE(lyr->qa, lyr->qb, w*lyr->thik);
      Ru = mltp(Ed, mltp(Ru, Ed));
      Tu = mltp(Tu, Ed);
      lyr->URu = Ru; lyr->UTu = Tu;	/* save the two */
      lyr++;
      L = mltp(invs(plus(I, ngtv(mltp(lyr->Q.rd, Ru)))), lyr->Q.tu);
      Ru = plus(lyr->Q.ru, mltp(lyr->Q.td, mltp(Ru, L)));
      Tu = mltp(Tu, L);
  }
  return(Tu);
}


/* calculate the displacement at free surface when layer k is perturbed. */
matrix	delrcv(float w, Layer *k)
{
  matrix Tu, Rd, L, Ed;
  Layer *below, *above;

  below = k+1;
  above = k-1;

  Tu = below->STu;
  Rd = below->SRd;
  L  = mltp(k->dQ2.tu, invs(plus(I, ngtv(mltp(Rd, k->dQ2.ru)))));
  Tu = mltp(L, Tu);
  Rd = plus(k->dQ2.rd, mltp(k->dQ2.tu, mltp(Rd, trns(L))));

  Ed = EEE(k->dqa, k->dqb, w*k->thik);
  Rd = mltp(Ed, mltp(Rd, Ed));
  Tu = mltp(Ed, Tu);

  L  = mltp(k->dQ1.tu, invs(plus(I, ngtv(mltp(Rd, k->dQ1.ru)))));
  Tu = mltp(L, Tu);
  Rd = plus(k->dQ1.rd, mltp(k->dQ1.tu, mltp(Rd, trns(L))));

  Tu = mltp(above->UTu, mltp(invs(plus(I, ngtv(mltp(Rd, above->URu)))), Tu));
  return(Tu);
}


/************************************************************************
  calculate the receiver function and the differentials w.r.t Vs (if db>0)
  and interface depth (if dh>0).
  Return Value: array of size nft*n, the first nft elements is the receiver
  		function and the rest are the differentials w.r.t Vs (and
		depth of lower interface) of each layer.
		n=1 (db=0); m (db>0&&dh=0); 2m-1 (db>0&&dh>0)
*************************************************************************/
float *partial( int		ps,	/* 0=P; 1=S receiver function */
		int		nft,	/* number of pts for fft */
		int		m,	/* number of input layers */
		const float	thik[],	/* thickness */
		const float	beta[], /* Vs */
		const float	kapa[],	/* Vp/Vs */
		float		p,	/* ray parameter */
		float 		dt,	/* sampling interval */
		float 		gauss,	/* Gaussian filter parameter */
		float 		shft,	/* time shift */
		float 		db,	/* amount of Vs perturbation */
		float 		dh	/* amount of thickness perturbation */
	)
{
  int     	i, j, k, n, nft2, nlyrs;
  unsigned char	vORh;
  float   	w, delw, agg, delta[2];
  complex 	temp, *a;
  Layer		*lyr;
  matrix  	dis;

  nft2 = nft/2;
  delw = 3.1415926/(dt*nft2);
  
  lyr = mdSetup(m, thik, beta, kapa, db, dh, &nlyrs);

  n = 1;
  if (db > 0.) { /* need to comput differential RF */
     delta[0] = delta[1] = 1./db;
     if (dh > 0.) delta[1] = -1./dh;
     n = nlyrs;
  }

  if ( (a = (complex *) malloc(nft2*n*sizeof(complex)))  == NULL ) {
     fprintf(stderr, "failed to allocat memeory for receiver fn\n");
     return NULL;
  }

  ifmat(p, nlyrs, lyr);
  if (db > 0.) delifm(p, nlyrs, lyr);
  for (i=0,w=0.; i<nft2; i++,w+=delw) {
      agg = 0.5*w/gauss;
      temp = cphase(cmplx(-agg*agg, -shft*w));
      dis = rcvrfn(w, nlyrs, lyr);
      if (ps==0) a[i] = cmltp(dis.sp, cinvs(cmltp(IMAGE, dis.pp)));
      else a[i] = cngtv(conjg(cmltp(cmltp(IMAGE,dis.ps), cinvs(dis.ss))));
      a[i] = cmltp(temp, a[i]);
      if (db > 0.) {
         rcvrtd(w, nlyrs, lyr);
         for (k=i+nft2,vORh=0,j=1; j<nlyrs; j++,k+=nft2) {
	     dis = delrcv(w,lyr+j);
	     if (ps==0) a[k] = cmltp(dis.sp, cinvs(cmltp(IMAGE, dis.pp)));
	     else a[k] = cngtv(conjg(cmltp(cmltp(IMAGE,dis.ps), cinvs(dis.ss))));
             a[k] = cmltp(temp, a[k]);
	     a[k] = dmltp(delta[vORh],cplus(a[k], cngtv(a[i])));
	     vORh ^= 1;
         }
      }
  }
  for (k=0,j=0; j<n; j++, k+=nft2) {
      fftr(&a[k], nft2, -dt);
  }

  free(lyr);

  return (float *) a;

}


/* calculate the radial and vertical displacement at free surface for P or S incident wave */
void	respknt(int		ps,	/* 0=incident P; 1=incident S */
		int		nft,	/* number of pts for fft */
		int		m,	/* number of input layers */
		const float	thik[],	/* thickness */
		const float	beta[], /* Vs */
		const float	kapa[],	/* Vp/Vs */
		float		p,	/* ray parameter */
		float 		dt,	/* sampling interval */
		complex		z[],	/* vertical response */
		complex		r[]	/* radial response */
	)
{
  int     	i, nft2, nlyrs;
  float   	w, delw;
  matrix  	aa;
  Layer		*lyr;

  nft2 = nft/2;
  delw = 3.1415926/(dt*nft2);

  lyr = mdSetup(m, thik, beta, kapa, 0., 0., &nlyrs);

  ifmat(p, nlyrs, lyr);
  for (i=0,w=0.; i<nft2; i++,w+=delw) {
      aa = rcvrfn(w, nlyrs, lyr);
      if (ps==0) {
         r[i] = aa.sp;
         z[i] = cmltp(IMAGE, aa.pp);
      } else {
         r[i] = aa.ss;
         z[i] = cmltp(IMAGE, aa.ps);
      }
  }
  fftr(r, nft2, -dt);
  fftr(z, nft2, -dt);
}


/* read in model parameters and return number of layers read */
int	mdin(const char *fo, float thik[], float beta[], float kapa[]) {
  int		m;
  char		line[128];
  FILE		*iopo;

  if ( (iopo = fopen(fo, "r")) == NULL ) {
     fprintf(stderr,"fail to open %s\n",fo);
     return 0;
  }
  m=0;
  while (fgets(line, 128, iopo)) {
     sscanf(line, "%f%f%f", thik+m, beta+m, kapa+m);
     if (m>MAXL) {
        fprintf(stderr,"too many layers in the model >%d\n",m);
	return 0;
     }
     m++;
     if ( thik[m-1] == 0. ) break;
  }
  fclose(iopo);

  return(m);
}
