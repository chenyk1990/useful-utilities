/*
 * rcvrfn.c:	calculate receiver functions from velocity model
 *
 * author	Lupei Zhu
 * revision history:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"
#include "knt.h"

int main (int argc, char **argv) {
  int		i, j, n, nt, nft, nlyrs, error=0, response=0, ps=0, com='r';
  float		p, dt, tt, gauss=5., shft=0., db=0., dh=0., *a, *z;
  float		beta[MAXL], thik[MAXL], kapa[MAXL];
  char		temp[64];
  SACHEAD	hd;

  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'G':
        sscanf(&argv[i][2],"%f",&gauss);
	break;
      case 'P':
        sscanf(&argv[i][2],"%f",&p);
        break;
      case 'Q':
	ps = 1;
	com = 'z';
	break;
      case 'R':
	response = 1;
	break;
      case 'S':
        sscanf(&argv[i][2],"%f",&shft);
        break;
      case 'T':
	sscanf(&argv[i][2],"%f/%f",&tt,&dt);
	break;
      case 'V':
        j = sscanf(&argv[i][2],"%f/%f",&db, &dh);
	if (j==1) dh = -1;
	break;
      default:
	error = 1;
	break;
      }
    }
  }

  if (error || argc == 1 ) {
     fprintf(stderr,"usage: %s -Tdura/dt -Pp [-Q] [-R[s] | -Sshift -Ggauss -Vdb[/dh]] models \n\
     It computes theoretical receiver functions for layered models (format is (thickness Vs Vp/Vs)).\n\
	-T: duration and sampling interval in sec. \n\
	-P: ray parameter in s/km \n\
	-Q: incident wave is S wave (P)\n\
	-R: output vertical and radial responses instead of recv. function. \n\
	-S: shift the output trace by shift sec. (0) \n\
	-G: Gaussion filter parameter in Hz. (5) \n\
	-V:  compute derivative of recv. functions using finite differences (no)\n\
	    db -> w.r.t. Vs in each layer: (r(b+db)-r(b))/db \n\
	    dh -> w.r.t. layer thickness: (r(h+dh)-r(h))/dh \n \
	    The order is top layer Vs/h, 2nd layer, ...\n",argv[0]);
     return -1;
  }

  nt = (int) (tt/dt);
  nft = 1; while (nft < nt) nft *= 2;
  hd = sachdr(dt,nt,0.);
  hd.user0 = p;
  hd.user1 = gauss;

  for (j=1; j<argc; j++) {
	if (argv[j][0] == '-') continue;
	fprintf(stderr,"%s\n",argv[j]);
	if ( (nlyrs = mdin(argv[j], thik, beta, kapa)) < 1 ||  nlyrs > MAXL ) continue;

	if (response) {
	   if ( (a=(float *)malloc(nft*sizeof(float))) == NULL ||
		(z=(float *)malloc(nft*sizeof(float))) == NULL ) {
		fprintf(stderr, "failed to allocate memeory\n");
		continue;
	   }
	   respknt(ps, nft, nlyrs, thik, beta, kapa, p, dt, (complex *) z, (complex *) a);
	   sprintf(temp, "%s.%05.3f.z", argv[j], p);
	   write_sac(temp, hd, z);
	   temp[strlen(temp)-1] = 'r';
	   write_sac(temp, hd, a);
	   free(a); free(z);
	}
	else {
  	   if ( (a = partial(ps, nft, nlyrs, thik, beta, kapa, p, dt, gauss, shft, db, dh)) == NULL) continue;
	   sprintf(temp, "%s.%05.3f.%cd", argv[j], p, com);
	   hd.b = -shft;
	   write_sac(temp, hd, a);
	   if (db>0.) {
	      n = nlyrs;
	      if (dh>0) n=2*nlyrs-1;
	      for(i=1;i<n;i++) {
	         sprintf(temp, "%s.%05.3f-%03d.%cd", argv[j], p, i, com);
		 write_sac(temp, hd, a+i*nft);
	      }
	   }
	   free(a);
	}
  }

  return 0;

}
