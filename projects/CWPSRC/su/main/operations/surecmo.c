/* Copyright (c) Colorado School of Mines, 2017.*/
/* All rights reserved.                       */

/* SURECMO: $Revision: 1.1 $ ; $Date: 2017/05/30 21:41:31 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" SURECMO - compensate for the continuously moving streamer in marine   ",
"           seismic acquisition (assume far offset is first channel)    ",
"           right now valid for data where tmax < 2*dx/vb               ",
"									",
" surecmo <stdin >sdout [vb= ] [dx= ] [dt= ] [fill= ] [tmax= ]          ",
"                                                                       ",
" vb=         boat speed in m/s                                         ",
" dx=         channel spacing in m                                      ",
" dt=         sample rate in microseconds                               ",
" fill=0.0    value to pad                                              ",
"									",
" Examples: 								",
" surecmo <stdin >sdout vb=2.5 dt=25 dt=4000                            ",
"                 (valid for tmax < 20s) 				",
NULL};

/* Author:
 *	CWP: Taylor Goss, May 2017
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/

segy tr;

int
main(int argc, char **argv)
{
	int nmix;		/* number of traces to mix over		*/
	int imix;		/* mixing counter			*/
	int it;			/* sample counter			*/
	int nt;			/* number of time samples per trace	*/
	int itr=0;		/* trace counter			*/
	size_t databytes;	/* number of bytes (nt*FSIZE)		*/
	size_t mixbytes;	/* number of bytes (nt*FSIZE*nmix)	*/
	float vb;               /* boat speed in m/s                    */
	float dx;               /* channel spacing in m                 */
	int dt;                 /* sample-rate in micro seconds         */
	float *mix=NULL;	/* array of mix values			*/
	float *temp=NULL;	/* temp array for mixing 		*/
	float **data=NULL;	/* array for mixing 			*/
	
	
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Get info from first trace */
	if(!gettr(&tr))
		err("can't get first trace");
	nt = tr.ns;
	
        if( !getparint("dt" ,&dt ) ){
           dt=tr.dt;
        };

	if (!getparfloat ("vb", &vb)) vb = 0;
	if (!getparfloat ("dx", &dx)) dx = 0;
	
	/* Get mix weighting values values */

	nmix = 3;
	mix = ealloc1float(nmix);

        checkpars();

	/* Compute databytes per trace and bytes in mixing panel */
	databytes = FSIZE*nt;
	mixbytes = databytes*nmix;

	/* Allocate temporary space for mixing  */
	data = ealloc2float(nt,nmix);
	temp = ealloc1float(nt);

	/* Zero out data array */
	memset((void *) data[0], 0, mixbytes);

	/* Loop over remaining traces */
	do {

		++itr;

		/* Zero out temp */
		memset((void *) temp, 0, databytes);

		/* Read data portion of trace into first column of data[][] */
		memcpy( (void *) data[0], (const void *) tr.data, databytes);
	
		/* Loop over time samples */
		for (it=0; it<nt; ++it) {
		  mix[1] = vb*it*dt/(dx*1e6);
		  mix[0] = 1.0 - mix[1];
		  mix[2] = 0.0;

		  if (mix[1] > 1) {
		    mix[2] = mix[1] - 1;
		    mix[1] = 1.0 - mix[2];
		    mix[0] = 0.0; };

		        /* Weighted moving average (mix) */
			for(imix=0; imix<nmix; ++imix)
				temp[it]+=data[imix][it]*mix[imix];

			/* put mixed data back in seismic trace */
			tr.data[it] = temp[it]; 
		}
		
                /* Bump columns of data[][] over by 1 */
                /* to make space for data from next trace */
                for (imix=nmix-1; 0<imix; --imix)
                        for (it=0; it<nt; ++it) 
                                data[imix][it] = data[imix-1][it];
		
		puttr(&tr);
	} while (gettr(&tr)); 

	return(CWP_Exit());

}
