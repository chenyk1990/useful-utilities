#include "su.h"
#include "segy.h"

#include <assert.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"								",
" sur2ai <stdin >sdout 						",
"       Reflectivity to Acoustic Impedance (Inversion)		",
"							        ",
" Required parameters: (none)				        ",
"							        ",
" Optional parameters: (none)				        ",
"							        ",
" Notes:						        ",
" Reflectivity traces in acoustic impedence out		        ",

"							        ",
NULL};

/*
 * Shuki Ronen: 2017
 *
 * Simple theory:
 * Reflectivity = Delta(Impedance) / Impedance
 *
 * Technical Reference:
 * Partial implementation of Aki & Richards, volume II, pages 661-662.
 * It's partial because transmission effects are ignored.
 * 
 *
 */


/**************** end self doc ***********************************/

segy r,ai;

int
main(int argc, char **argv)
{
	int nt;			/* number of points on input trace	*/
	int itr;		/* counter				*/
	int verbose = 0;
	
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	getparint("verbose",&verbose);

	/* Get info from first trace */
	itr = gettr(&r);
	if(itr==0) {
		warn ("can't get first trace");
	}
	assert(itr);

	nt = r.ns;

	bcopy(&r,&ai,240);

	/* Loop over traces */
	itr = 1;			/* First trace already in */
	do {
		int it;

		ai.data[0] = 1.0;
		for(it=0;it<nt-1;it++) {
			/* if(r.data[it]>=1.0 || r.data[it]<= -1.0) */
			if(r.data[it]>1.0 || r.data[it]< -1.0)
				err(" Illegal Reflectivity=%g ... abs value must be <1\n",r.data[it]);
			/* ai.data[it+1] = ai.data[it]*(1.0+r.data[it])/(1.0-r.data[it]); */
			ai.data[it+1] = ai.data[it]*(1.0-r.data[it])/(1.0+r.data[it]);
			/* fprintf(stderr,"ai.data[%d]=%g ai.data[%d]=%g r.data[%d]=%g\n", */
			/*	it+1,ai.data[it+1], it,ai.data[it], it,r.data[it]); */
		}

		puttr(&ai);

		itr++;

	} while (gettr(&r));

	return 0;
}
