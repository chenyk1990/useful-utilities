#include "su.h"
#include "segy.h"

#include <assert.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"								",
" suai2r <stdin >sdout 						",
"	Acoustic Impedance to Reflectivity (Forward Modeling)   ",
"							        ",
"							        ",
" Required parameters: (none)				        ",
"							        ",
" Optional parameters: (none)				        ",
"							        ",
" Notes:						        ",
"  Acoustic impedence traces in Reflectivity out 		",
"							        ",
NULL};

/*
 *   Shuki Ronen: 2017
 *
 * Simple theory:
 * Reflectivity = Delta(Impedance) / Impedance
 *
 * Technical reference:
 * Partial implementation of Aki & Richards, volume II, pages 661-662.
 * It's partial because transmission effects are ignored.
 *
 */

/**************** end self doc ***********************************/

segy tr;

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
	itr = gettr(&tr);
	if(itr==0) {
		warn ("can't get first trace");
	}
	assert(itr);

	nt = tr.ns;

	/* Loop over traces */
	itr = 1;			/* First trace already in*/
	do {
		int it;

		for(it=0;it<nt-1;it++)
			tr.data[it] = -(tr.data[it+1]-tr.data[it])/(tr.data[it+1]+tr.data[it]);
		tr.data[nt-1] = 0.0;

		puttr(&tr);

		itr++;

	} while (gettr(&tr));

	return 0;
}
