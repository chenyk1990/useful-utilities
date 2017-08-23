/* Copyright (c) Durham University, 2013.*/
/* All rights reserved.		       */

/* SUSVDFILT: $Revision: 1.1 $ ; $Date: 2016/04/05 19:44:27 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation *****************************/
char *sdoc[] = {
"								       ",
" SUSVDFILT - SVD (Eigen) filter					",
"								       ",
" susvdfilt <stdin >stdout [optional parameters]			",
"								       ",
" Required parameters:							",
"    none							       ",
"								       ",
" Optional parameters:						  ",
"    ntrw=all traces   number of traces in window		       ",
"    verbose=0	 1 = echo additional information		  ",
"								       ",
"    numpp=1	   number of principal planes to retain	     ",
"    subtract=1	subtarct selected principal planes from i/p      ",
"	     0	output selected principal planes		 ",
"									",
"    tmpdir=	   if non-empty, use the value as a directory path  ",
"		       prefix for storing temporary files; else if the  ",
"		       the CWP_TMPDIR environment variable is set use   ",
"		       its value for the path; else use tmpfile()       ",
" 									",
"								       ",
" Notes:								",
"    Input data is windowed to give an area of twlen x ntrw samples     ",
"    which is then decomposed to into eigen vectors (by SVD) and	",
"    eigenvalues. A percentage of these are selected for inverse	",
"    transform this can be subtracted from input if difference is       ",
"    required.							  ",
"								       ",
" Caveat:							       ",
" The pre-requiste is your target event has to be pre-flattened using NMO,",
" LMO or other static correction. Works well to remove direct wave from",
" marine seismic data for those of you interested in seismic oceanography.",
"								       ",
" The code really needs someone to make the windowing work in both time ",
" and space. However as it stands it works, it has excellent amplitude 	",
" response and exhibits minimal edge effects.			   ",
"								       ",
NULL};

/* 
 * Author: Richard Hobbs
 *	 Durham University
 *	 E-mail: r.w.hobbs@durham.ac.uk
 *
 *
 * Based on sueipofi.c
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 */
/**************** end self doc *******************************************/

static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

segy tr;

int
main(int argc, char **argv)
{	
    
    register int i;       /* indices for components (in loops) */
    register int j,k;     /* loop index */
    int nsp;	      /* number of samples in input trace	*/
    int ntrc;		  /* number of input traces		*/
    int ntrw;		  /* number of traces in window		*/
    float dt;	     /* sampling intervall in seconds */
    int numpp;	    /* number of principal planes to retain */
    int subtract;	 /* flag for subtracing eigen image from input */
    float **data;	 /* input data */
    float **udata;	/* windowed data and matrix U after SVD */
    float **vdata;	/* matrix V */
    float *sval;	  /* array of singular values */
    float **c;	    /* work space */
    float **out=NULL;     /* array[nx][nt] of output traces */	
    int verbose;	  /* flag for echoing information */
    char *tmpdir=NULL;	  /* directory path for tmp files */
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
    nsp = tr.ns;
    dt = (float) tr.dt/1000000.0;

    /* Store traces in tmpfile while getting a counsp */
    if (STREQ(tmpdir,"")) {
	tracefp = etmpfile();
	headerfp = etmpfile();
	if (verbose) warn("using tmpfile() call");
    } 
    else { /* user-supplied tmpdir */
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

    ntrc = 0;
    do {
	++ntrc;
	efwrite(&tr, 1, HDRBYTES, headerfp);
	efwrite(tr.data, FSIZE, nsp, tracefp);
    } while (gettr(&tr));

    /* get general flags and parameters and set defaults */
    if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
    if (!getparint("ntrw", &ntrw))       ntrw = ntrc;
    if (ntrw == 0) {
	ntrw = ntrc;
	if (verbose) warn("ntrw cannot be zero, using all traces");
    }
    if (!getparint("subtract", &subtract)) subtract=1;
    if (!getparint("numpp", &numpp))       NINT(numpp=ntrc/100);
    if (numpp == 0 ) {
	numpp=1;
	warn("numpp=0 set to 1");
    }

    checkpars();

    if (dt == 0.0) err("header field dt not set, must be getparred");

	
    /* data validation */
    if (ntrw<3) err("ntrw=%g must be >= 3 traces. (3 samples)", ntrw);
    if (ntrw>ntrc) err("ntrw=%g larger then number of available traces %g", ntrw, ntrc);
    
    /* echo window information */
    if (verbose) {
	warn("nsp=%d ntrc=%d ntrw=%d\n",nsp,ntrc,ntrw);
    }   

    /* allocate space */
    data = alloc2float(nsp,ntrc);   /* input data */
    udata = alloc2float(ntrw,nsp);  /* SVD matrix U */
    vdata = alloc2float(ntrw,ntrw);  /* SVD matrix V */
    sval = alloc1float(ntrw);       /* singular values */ 
    out = alloc2float(ntrc,nsp);
    c = alloc2float(ntrw,ntrw);      /* work-space */ 
    
    /* load traces into an array */
    rewind(headerfp);
    rewind(tracefp);
    for (i=0; i<ntrc; i++) {
	fread (data[i], FSIZE, nsp, tracefp);
    }
    rewind(headerfp);
    rewind(tracefp);

    if (verbose) fprintf(stderr,"Processing");

   
    /* copy data into matrix for SVD */
    for (i=0; i<ntrw; i++) {
	for (j=0; j<nsp; j++) {
	    udata[j][i]=data[i][j];
	}
    }

    /* perform singular value decomposition (SVD) */
    compute_svd(udata, nsp, ntrw, sval, vdata);
    svd_sort(udata, sval, vdata, ntrw, nsp);
 
    /* reconstruct output */
    for (j=0; j<ntrw; j++) {
	for (i=0; i<ntrw; i++) {
	    c[j][i]=vdata[i][j]*sval[j];
	}
    }

    for (i=0; i<ntrw; i++) {
	for (j=0; j<nsp; j++) {
	    out[j][i]=0;
	    for (k=0; k<numpp; k++) {
		out[j][i]=out[j][i]+(udata[j][k]*c[k][i]);
	    }
	}
    }
	     
    for (i=0; i<ntrc; i++) {
	fread(&tr, 1, HDRBYTES, headerfp);
	if (!subtract) {
	    for (j=0; j<nsp; j++) {
		tr.data[j]=out[j][i];
	    }
	}
	else {
	    for (j=0; j<nsp; j++) { 
		tr.data[j]=data[i][j]-out[j][i];
	    }
	}
	puttr(&tr);
    }


    if (verbose) {
	warn("done");
    }

    fclose (headerfp);
    fclose (tracefp);
    if (istmpdir) {
	remove(tracefile);
	remove(headerfile);
    }

    return(CWP_Exit());
}


/* for graceful interrupt termination */
static void closefiles(void)
{
	fclose(headerfp);
	fclose(tracefp);
	remove(headerfile);
	remove(tracefile);
	exit(EXIT_FAILURE);
}

/* END OF FILE */
