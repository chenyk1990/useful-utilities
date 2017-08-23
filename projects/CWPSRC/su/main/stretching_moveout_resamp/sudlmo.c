/* Copyright (c) Colorado School of Mines, 1998.*/
/* All rights reserved.                       */

/* SUDLMO: $Revision: 1.1 $ ; $Date: 2016/04/06 17:45:36 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" SUDLMO -- Dynamic Linear Move Out Correction for Surface Waves  	",
"									",
" sudlmo <indata >outdata [optional parameters]				",
"									",
" Optional parameters:							",
"	vnmo=400	Velocity corresponding to  fnmo			",
"	fnmo=0.0	Frequency corresponding to  vnmo		",
"									",
"	dt=(from header)	time sample interval (in seconds)	",
"       cm = 0 			for offset in m   m                     ",
"       		       1 for offset in cm                     ",
"       invert = 0             1 to perform invers DLMO        		",
"									",
" Algorithm:								",
"	Gdlmo(t,x) = Re[INVFTT{ ( (sign) iw*off/Cf*(Gdlmo(f,x)}]	",
"									",
NULL};

/* Credits:
 *	BRGM: Adnand Bitri, oct 1998; based on sufrac
 *
 * Trace header fields accessed: ns, dt, offset
*/
/**************** end self doc ***********************************/


#define	I		cmplx(0.0, 1.0)
#define	PIBY2		0.5 * PI
#define TWOPI		2.0 * PI
#define LOOKFAC		2	/* Look ahead factor for npfao	  */
#define PFA_MAX		720720	/* Largest allowed nfft	          */

static void interpovv (int nt, int ncdp, float *cdp,
        float **ovv,  float cdpt, float *ovvt);


segy tr;

int
main(int argc, char **argv)
{
	float phase;		/* phase shift = phasefac*PI		*/
	float power;		/* phase shift = phasefac*PI		*/
	register float *rt;	/* real trace				*/
	register complex *ct;	/* complex transformed trace		*/
	complex *filt;		/* complex power	 		*/
	int nt;			/* number of points on input trace	*/
	size_t ntsize;		/* nt in bytes				*/
	int ncdp;               /* number of cdps specified */
	int icdp;       	/* index into cdp array */

	long oldoffset;         /* offset of previous trace */
        long oldcdp;    	/* cdp of previous trace */
        int newsloth;   	/* if non-zero, new sloth function was computed */
	int jcdp;       	/* index into cdp array */
	float dt;		/* sample spacing (secs) on input trace	*/
	float tn;		/* sample spacing (secs) on input trace	*/
	float omega;		/* circular frequency			*/
	float domega;		/* circular frequency spacing (from dt)	*/
	int nfft;		/* number of points in nfft		*/
	int ntnmo;      	/* number of tnmos specified            */
	float *cdp;     	/* array[ncdp] of cdps */

	float *vnmo;   		/* array[nvnmo] of vnmos               */
	float *ovvt;   		/* array[nvnmo] of vnmos               */
	int nvnmo;      	/* number of tnmos specified            */
	float *fnmo;   		 /* array[ntnmo] of tnmos               */
	float **ovv;   		 /* array[nf] of fnmos                  */
	float doffs;             /* offset                            */
	float acdp;     	/* temporary used to sort cdp array */
        float *aovv;    	/* temporary used to sort ovv array */

        int invert;              /*  if non-zero, do invers DLMO       */
        int cm;                  /*  if non-zero, the offset in cm     */
        int nf;                 /* number of frequencies (incl Nyq)     */
        int it;                 /* number of frequencies (incl Nyq)     */
	float onfft;		/* 1 / nfft				*/
	float v;                 /* velocity                            */
	size_t nzeros;		/* number of padded zeroes in bytes	*/
	
	
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Set parameters */
	power=0.0;

	/* Get info from first trace*/ 
	if (!gettr(&tr))	err("can't get first trace");

	nt = tr.ns;

	if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
	if (!dt)	err("dt field is zero and not getparred");
	ntsize = nt * FSIZE;

        if (!getparint("invert",&invert)) invert = 0;
        if (!getparint("cm",&cm)) cm = 0;

	/* Set up for fft */
	nfft = npfaro(nt, LOOKFAC * nt);
	if (nfft >= SU_NFLTS || nfft >= PFA_MAX)
		err("Padded nt=%d -- too big", nfft);

        nf = nfft/2 + 1;
        onfft = 1.0 / nfft;
	nzeros = (nfft - nt) * FSIZE;
	domega = TWOPI * onfft / dt;


	/* get velocity functions, linearly interpolated in frequency */

	ncdp = countparval("cdp");

	if (ncdp>0) {
                if (countparname("vnmo")!=ncdp)
                        err("a vnmo array must be specified for each cdp");
                if (countparname("fnmo")!=ncdp)
                        err("a tnmo array must be specified for each cdp");
        } else {
                ncdp = 1;
                if (countparname("vnmo")>1)
                        err("only one (or no) vnmo array must be specified");
                if (countparname("fnmo")>1)
                        err("only one (or no) tnmo array must be specified");
        }

	cdp = ealloc1float(ncdp);
        if (!getparfloat("cdp",cdp)) cdp[0] = tr.cdp;
        ovv = ealloc2float(nf,ncdp);

	for (icdp=0; icdp<ncdp; ++icdp) {
                nvnmo = countnparval(icdp+1,"vnmo");
                ntnmo = countnparval(icdp+1,"fnmo");
                if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
                        err("number of vnmo and tnmo values must be equal");
                if (nvnmo==0) nvnmo = 1;
                if (ntnmo==0) ntnmo = nvnmo;
                /* equal numbers of parameters vnmo, fnmo  */

                vnmo = ealloc1float(nvnmo);
                fnmo = ealloc1float(nvnmo);

                if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 400.0;
                if (!getnparfloat(icdp+1,"fnmo",fnmo)) fnmo[0] = 0.0;
		

		for (it=0; it<ntnmo; ++it)
			fnmo[it]*=TWOPI;

                for (it=1; it<ntnmo; ++it)
                        if (fnmo[it]<=fnmo[it-1])
                                err("tnmo values must increase monotonically");

		for (it=0,tn=0; it<nf; ++it,tn+=domega) {
			intlin(ntnmo,fnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&v);
			ovv[icdp][it] = 1.0/(v); 
		}
                free1float(vnmo);
                free1float(fnmo);
        }



/* sort (by insertion) sloth and anis functions by increasing cdp */

        for (jcdp=1; jcdp<ncdp; ++jcdp) {
                acdp = cdp[jcdp];
                aovv = ovv[jcdp];
                for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
                        cdp[icdp+1] = cdp[icdp];
                        ovv[icdp+1] = ovv[icdp];
                }
                cdp[icdp+1] = acdp;
                ovv[icdp+1] = aovv;
        }

/* allocate workspace */


        ovvt = ealloc1float(nf);

/* interpolate sloth and anis function for first trace */

        interpovv(nf,ncdp,cdp,ovv,(float)tr.cdp,ovvt);

        /* set old cdp and old offset for first trace */
        oldcdp = tr.cdp;
        oldoffset = tr.offset-1;



	/* Allocate fft arrays */
	rt   = ealloc1float(nfft);
	ct   = ealloc1complex(nf);
	filt = ealloc1complex(nf);


		

	/* Loop over traces */
	do {

		/* if necessary, compute new sloth and anis function */
                	if (tr.cdp!=oldcdp && ncdp>1) {
                        interpovv(nt,ncdp,cdp,ovv,(float)tr.cdp,
                                  ovvt);
                        newsloth = 1;
                } else {
                        newsloth = 0;
                }

		/* if sloth and anis function or offset has changed */

                if (newsloth || tr.offset!=oldoffset) {

		doffs = (fabs)((float)(tr.offset));
		if (cm==1) doffs/=100;
		/* Load trace into rt (zero-padded) */
		memcpy( (void *) rt, (const void *) tr.data, ntsize);
		memset((void *) (rt + nt), (int) '\0', nzeros);

		/* FFT */
		pfarc(1, nfft, rt, ct);


		/* Apply filter */
		{ register int i;
		for (i = 0; i < nf; ++i){
			omega = i * domega;
			if (power < 0 && i == 0) omega = FLT_MAX;
			if (invert==0)
			phase = -1.0*omega*ovvt[i]*doffs;
			else
			phase = 1.0*omega*ovvt[i]*doffs;
			
		/*	filt[i] = cmplx(cos(phase),sin(phase)); */
			filt[i] = cwp_cexp(crmul(I,phase)); 
			filt[i] = crmul(filt[i], onfft);
			ct[i] = cmul(ct[i], filt[i]);
			
			}
		}
	  }
		/* Invert */
		pfacr(-1, nfft, ct, rt);


		/* Load traces back in, recall filter had nfft factor */
		{ register int i;
		for (i = 0; i < nt; ++i)  tr.data[i] = rt[i];
		}


		puttr(&tr);

	} while (gettr(&tr));


	return EXIT_SUCCESS;
}


/* linearly interpolate/extrapolate sloth and anis between cdps */
static void interpovv (int nt, int ncdp, float *cdp, float **ovv,
        float cdpt, float *ovvt)
{
        static int indx=0;
        int it;
        float a1,a2;

        /* if before first cdp, constant extrapolate */
        if (cdpt<=cdp[0]) {
                for (it=0; it<nt; ++it) {
                        ovvt[it] = ovv[0][it];
                      };

        /* else if beyond last cdp, constant extrapolate */
        } else if (cdpt>=cdp[ncdp-1]) {
                for (it=0; it<nt; ++it) {
                        ovvt[it] = ovv[ncdp-1][it];
                      };

        /* else, linearly interpolate */
        } else {
                xindex(ncdp,cdp,cdpt,&indx);
                a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
                a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
                for (it=0; it<nt; ++it) {
                        ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
                      };
        }
}

