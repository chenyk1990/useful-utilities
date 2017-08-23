/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* SUGROLL: $Revision: 1.1 $ ; $Date: 2016/04/05 15:37:43 $	*/

#include <cwp.h>
#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SUGROLL - Ground roll supression using Karhunen-Loeve transform	",
"                                                                       ",
"    sugroll <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
" dt=tr.dt (from header) 	time sampling interval (secs)           ",
" nx=ntr   (counted from data)	number of horizontal samples (traces)	",
" sb=0      1 for the graund-roll                                       ",
" verbose=0	verbose = 1 echoes information			",
" nrot=3        the principal components for the m largest eigenvalues ",
" tmpdir= 	 if non-empty, use the value as a directory path",
"		 prefix for storing temporary files; else if the",
"	         the CWP_TMPDIR environment variable is set use	",
"	         its value for the path; else use tmpfile()	",
" 								",
" Notes:                                                                ",
" The method developed here is to extract the ground-roll from the	",
" common-shot gathers using Karhunen-Loeve transform, and then to substract  ",
" it from the original data. The advantage is the ground-roll is suppresed  ",
" with negligible distortion of the signal.	",
"                                                                       ",
NULL};

/*
 * Credits: BRGM: Adnand Bitri, 1999.
 *
 * Reference:       
 *    Xuewei Liu, F., 1999, Ground roll supresion using the Karhunen-Loeve transform 
 *         Geophysics vol. 64 No. 2 pp 564-566
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: dt
 */
/**************** end self doc ********************************/

static void closefiles(void);
static void jacobi(float **a, int n, float w[], float **v,int *nrots);
static void eigsrt(float d[], float **v, int n);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/


segy tr;

int
main(int argc, char **argv)
{
	int ix,it,i,j,k,nrots;		/* loop counters */
	int ntr;		/* number of input traces */
	int nt;			/* number of time samples */
	int nx,sb,nrot;		/* number of horizontal samples */
	float dt,singval;               /* Time sample interval */
	float **in_traces;	/* array[nx][nt] of input traces */	
	float **out_traces;	/* array[nx][nt] of output traces */	
	float **XT,**G,**v,*w,**u;
	int verbose;		/* flag for echoing information */
	char *tmpdir;		/* directory path for tmp files */
	cwp_Bool istmpdir=cwp_false;/* true for user-given path */
	
        /* hook up getpar to handle the parameters */
        initargs(argc,argv);
        requestdoc(1);

	if (!getparint("verbose", &verbose))	verbose = 0;

	/* Look for user-supplied tmpdir */
	if (!getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		err("you can't write in %s (or it doesn't exist)", tmpdir);

        /* get info from first trace */
        if (!gettr(&tr))  err("can't get first trace");
        nt = tr.ns;
        dt = (float) tr.dt/1000000.0;
         if (!getparint("sb",&sb)) sb = 0;
         if (!getparint("nrot",&nrot)) nrot = 3;
        /* Store traces in tmpfile while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}
        ntr = 0;
        do {
                ++ntr;
                efwrite(&tr, 1, HDRBYTES, headerfp);
                efwrite(tr.data, FSIZE, nt, tracefp);
        } while (gettr(&tr));

        /* get general flags and parameters and set defaults */
        if (!getparint("nx",&nx))          	nx = ntr;
	if (!getparfloat("dt",&dt))		dt = dt;


	

	if (dt == 0.0)
		err("header field dt not set, must be getparred");

	/* allocate space */
        in_traces = alloc2float(nt, ntr);
        out_traces = alloc2float(nt, ntr);
        XT = alloc2float(ntr, nt);	
        G = alloc2float(ntr, ntr);	
        v = alloc2float(ntr, ntr);	
        u = alloc2float(nt, ntr);	
	w = alloc1float(ntr);

        /* load traces into an array and close temp file */
	erewind(headerfp);
        erewind(tracefp);
        for (ix=0; ix<ntr; ix++)
                fread (in_traces[ix], FSIZE, nt, tracefp);
        efclose (tracefp);


	/* Transpose matrix of the data XT */

	for(ix=0; ix<ntr; ix++)
		for(it=0; it<nt; it++)
			XT[it][ix] = in_traces[ix][it];

	/* Covariance matrix constaraction G=X*XT */

		for(i=0; i<ntr; i++){
			for(j=0; j<ntr; j++){
				G[i][j] = 0.0;
			for(k=0; k<nt; k++)
				G[i][j] +=in_traces[i][k]*XT[k][j];
			}
		}


	/*     transposition verification */


		for(i=0; i< ntr; i++)
			for(j=0; j< ntr; j++)

				if(G[i][j] != G[j][i]) {
			warn(" G n'est pas symetrique ");
			}
	/* Perform Sisgular Value Decomposition (svd) for a symetric matrix G */

		jacobi(G,ntr,w,v,&nrots);
                 eigsrt(w,v,ntr); 

		for(i=0; i<ntr; i++){
			singval = w[i]/w[0];
			warn(" Singular value number singval=%g", singval);
			}
	
			
	/*Psi matrix construction by multiplication v-transpose with the data */

		for(i=0; i<ntr; i++){
                        for(j=0; j<nt; j++){
                                u[i][j] = 0.0;
                        for(k=0; k<ntr; k++)
                                u[i][j] +=v[k][i]*in_traces[k][j];
                        }
                }

	/*Psi' matrix construction by selecting the m top rows of Psi */
	/* X' matrix construction by X' = R*Psi'                      */



		for( i=nrot; i<nx; i++)
			for(j=0; j<nt; j++)
				u[i][j] = 0.0; 
/***************************************************************/

		for(i=0; i<ntr; i++){
                        for(j=0; j<nt; j++){
                                out_traces[i][j] = 0.0;
                        for(k=0; k<ntr; k++)
                                out_traces[i][j] +=v[i][k]*u[k][j];
                        }
                }

	if (istmpdir) eremove(tracefile);

		
        /* write output traces */
        erewind(headerfp);
	{       register int itr;
		for (itr=0; itr<ntr; itr++) {
			efread(&tr, 1, HDRBYTES, headerfp);
			for (it=0; it<nt; it++){ 
			       if(sb==0)
				tr.data[it]=in_traces[itr][it]-out_traces[itr][it];
				else
				tr.data[it]=out_traces[itr][it];
				
				}
			puttr(&tr);
		}
	}
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);

	/* free allocated space */
	free2float(in_traces);
	free2float(out_traces);
	free2float(G);
	free2float(XT);
	free2float(v);
	free1float(w);

	return EXIT_SUCCESS;

}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

static void eigsrt(float d[], float **v, int n)
{
	int k,j,i;
	float p;
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i+1;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

static void jacobi(float **a, int n, float d[], float **v, int *nrots)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	b=alloc1float(n);
	z=alloc1float(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrots=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip< n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)  
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free1float(z);
			free1float(b);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrots);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
fprintf(stderr,"To many iterations in routine jacobi\n");exit (-2);
}
