/* Bayesshrink.
*/
/*
  Copyright (C) 2014 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, n, ifverb;
    float *dat=NULL, *adat=NULL, *diff, t, d, sigma_d2,sigma_n,sigma_s;
    char *type, *thre;
    sf_complex *cdat=NULL, *cdiff;
    sf_file in=NULL, out=NULL, other=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    n = sf_filesize(in);
    adat = sf_floatalloc(n);

    if (NULL == (type=sf_getstring("type"))) type="soft";
    /* [soft,hard] thresholding type, the default is soft  */

    if(!sf_getint("ifverb",&ifverb)) ifverb=0;
    /* 0, not print threshold value; 1, print threshold value. */

    if (NULL != (thre=sf_getstring("other"))){other = sf_output("other");}
    /* If output the difference between the thresholded part and the original one */

    	if (SF_FLOAT == sf_gettype(in)) {
		dat = sf_floatalloc(n);
		sf_floatread(dat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = fabsf(dat[i]);
		}
    	} else if (SF_COMPLEX == sf_gettype(in)) {
		cdat = sf_complexalloc(n);
		sf_complexread(cdat,n,in);
		for (i=0; i < n; i++) {
	    		adat[i] = cabsf(cdat[i]);
		}
    	} else {
		sf_error("Need float or complex input");
    	}
	/* adat[i] is the absolute value of the input sequence */

	
	sigma_n = sf_quantile(floorf(n/2),n,adat) / 0.6745;
	sigma_d2=0.0;
	for(i=0;i<n;i++)
		sigma_d2 += 1.0/n * adat[i]*adat[i];
	sigma_s = sqrtf(SF_MAX(sigma_d2 - sigma_n*sigma_n, 0));
	if(sigma_s!=0) t=2*sigma_n*sigma_n/sigma_s;
	else	t=0;
    
    if(ifverb==1) sf_warning("Threshold=%g",t);   	

    if (NULL != dat) {
	if(thre !=NULL ) {diff=sf_floatalloc(n);}
	for (i=0; i < n; i++) {
	    d = dat[i];
	    if (d < -t) {
		if(type[0]=='s') dat[i] = d+t; 
		if(thre!=NULL) diff[i]=d-dat[i];
	    } else if (d > t) {
		if(type[0]=='s') dat[i] = d-t;
		if(thre!=NULL){ diff[i]=d-dat[i];}
	    } else {
		dat[i] = 0.;
		if(thre!=NULL) diff[i]=d;
	    }
	}
	if(thre != NULL) 
		sf_floatwrite(diff,n,other); /* write the difference */		
	sf_floatwrite(dat,n,out);
    } else {
	if(thre !=NULL ){ cdiff=sf_complexalloc(n);}
	for (i=0; i < n; i++) {
	    d = cabsf(cdat[i]);
	    if (d < -t) {
#ifdef SF_HAS_COMPLEX_H
		if(thre!=NULL)   cdiff[i] = cdat[i]-cdat[i]*(d+t)/d;
		if(type[0]=='s') cdat[i] *= (d+t)/d;
#else
		if(thre!=NULL)   cdiff[i] = sf_csub(cdat[i],sf_crmul(cdat[i],(d+t)/d));
		if(type[0]=='s') cdat[i] = sf_crmul(cdat[i],(d+t)/d);
#endif
	    } else if (d > t) {		
#ifdef SF_HAS_COMPLEX_H
		if(thre!=NULL)   cdiff[i] = cdat[i]-cdat[i]*(d-t)/d;
		if(type[0]=='s') cdat[i] *= (d-t)/d;
#else
		if(thre!=NULL)   cdiff[i] = sf_csub(cdat[i],sf_crmul(cdat[i],(d-t)/d));
		if(type[0]=='s') cdat[i] = sf_crmul(cdat[i],(d-t)/d);
#endif
	    } else {
		if(thre!=NULL)   cdiff[i] = cdat[i];
		cdat[i] = sf_cmplx(0.,0.);
	    }
	}
	if(thre!=NULL) sf_complexwrite(cdiff,n,other);
	sf_complexwrite(cdat,n,out);
    }

    exit(0);
}
