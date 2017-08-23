/* rayt3d: $Revision: 1.0 $ ; $Date: 2004/05/27 20:26:14 $	*/
#include "usgrid.h"
#include "grid.h"
#include "su.h"
#include "header.h"


char *sdoc[] = {
"rayt3d - traveltime tables calulated by paraxial ray tracing ",
" ",
"rayt3d vfile= tfile= [optional parameters] ",
" ",
"Required Parameters:",
"vfile=stdin            file containing velocitiy v(nz,nx,ny)    	",
"tfile=stdout           file containing traveltime tables  ",
"			t(nzo,nxo,nyo,nxs,nys)	    ",
" ",
"Optional Parameters:",
"dt=0.008  		time sampling interval in ray tracing	",
"nt=401  		number of time sampling in ray tracing	",
"                       (one-way travel time calculation)	",
" ",
"the following nine parameters are from velocity grid header ONLY",
"fz=from-vfile	        first depth sample in velocity 	",
"nz=from-vfile       	Number of depth samples in velocity 	",
"dz=from-vfile       	depth interval in velocity		",
"fx=from-vfile       	first inline coordinate in velocity 	", 
"nx=from-vfile		Number of inline coordinates in velocity	",
"dx=from-vfile		inline interval in velocity			",
"fy=from-vfile       	first crossline coordinate in velocity 	" ,
"ny=from-vfile       	Number of crossline coordinates in velocity	",
"dy=from-vfile        	crossline interval in velocity			",
" ",
"fxo=fx                 first x coordinate in traveltime table",
"fyo=fy                 first y coordinate in traveltime table",
"fzo=fz                 first z coordinate in traveltime table",
"dxo=dx                 x interval in traveltime table",
"dyo=dy                 y interval in traveltime table",
"dzo=dz                 z interval in traveltime output",
"nxo=nx                 number of x samples in traveltime table",
"nyo=ny                 number of y samples in traveltime table",
"nzo=nz                 number of z samples in traveltime table",
"ntline=0		number of target lines (inlinel, crossline or	",
"			slant line). Traveltime outputs on each target 	",
"			line instead of whole volume			",
"when ntline>0, specify 			",
"	v0=5000		reference velocity at the surface		",
"	dvz=0		reference velocity vertical gradient		",
"and for each target line specify					",
"	xbeg		x coordinate at which the line begins		",
"	xend		x coordinate at which the line ends		",
"	ybeg		y coordinate at which the line begins		",
"	yend		y coordinate at which the line ends		",
"	nout		number of output points in the line		",
"	tofile		filename containing traveltime output  		",
"                       tfile will be ignored				",
"                       the above 6 parameters must be given ntline times ",
" ",
"nxs=1                  number of x samples for source locations	",
"nys=1                  number of y samples for source locations	",
"fxs=fx                 x coordinate of first source ",
"fys=fy                 y coordinate of first source ",
"dxs=2*dxo              x coordinate increment of sources ",
"dys=2*dyo              y coordinate increment of sources ",
"aperx=0.5*nx*dx        ray tracing aperature in x-direction ",
"apery=0.5*ny*dy        ray tracing aperature in y-direction ",
" ",
"fa=0            	first take-off polar angle ",
"da=1            	increment of take-off polar angle ",
"na=61            	number of polar angles",
"amin=0            	minimum emergence polar angle  ",
"                       measured from vertical --- moving upwards	",
"amax=90            	maximum emergence polar angle	",
"azhmin=0            	minimum take-off azimuth angle ",
"                       measured from y axis (clockwise)	",
"azhmax=360            	maximum take-off azimuth angle",
" ",
"fac=0.01            	factor to determine radius for extrapolation	",
"ek=1              	=1 to implement eikonal equation in shadow zones",
"ms=1			print verbal information at every ms finished sources",
"ncpu=1			number of cpu to use in the computation ",
"jpfile=stderr          name of job print file; default to standard error",
"restart=n		job is restarted (=y yes; =n no)		",
"isres=                 restart source index (sequential) ",
"                       (default will be determined from the travel ",
"                        file size) ",
" ",
"Note:	",
"1. Each traveltime table is calculated by paraxial ray tracing; then 	",
"   traveltimes in shadow zones are computed by solving eikonal equation.",
"2. A smoothed velocity is prefered. All sampling information of the	",
"   velocity is included in the file grid header.			",
"3. Traveltime table and source locations must be within velocity model. ",
"4. Ray tracing apertures can be choosen as sum of migration apertures ",
"   plus half of maximum offset. ",
"5. When ntline=0, the program outputs the whole traveltime tables.",
"   Otherwise, the outputs are given at the single target lines.	",
"6. When a single target line is slant, the specified parameters for  	",
"   output must be consistent with those used in migration. 	",	
"7. memory requirement for this program is about,		",
"   (7*nz*(nx*ny+ncpu*nxyt)+(1+2*cpu)*nxo*nyo*nzo+28*128*1001)*4 bytes	",
"   where 1001 is the maximum number of nt, 128 is maximum number of	",
"   rays to compute at vector unit, 28 is the number of auxiliary arrays",
"   used in ray tracing, nxyt=(2*aperx/dx)*(2*apery/dy)			",
" ",
" ",
NULL};

void rt3dp_(int *np,int *nzyx,int *nzyxo,int *nz,int *ny,int *nx,
  int *nzo,int *nyo,int *nxo,int *ek,int *nt, int *na,
  float *fa,float *da,float *amin,float *amax,float *azhmin,float *azhmax,
  float *dt,float *tmax,
  float *fxo,float *fyo,float *fzo,float *dxo,float *dyo,float *dzo,
  float *fac,float *fx,float *fy,float *fz,
  float *dx,float *dy,float *dz,float *aperx,float *apery,
  float *v,float *vxx,float *vxy,float *vxz,float *vyy,float *vyz,
  float *vzz,float *ov2,
  float *vt,float *vxxt,float *vxyt,float *vxzt,float *vyyt,
  float *vyzt,float *vzzt,
  float *tt1,float *tt2,float *t,float *s,
  float *xsp,float *ysp,float *azhnp,float *azhxp,float *fxtp,float *fytp,
  int *nxtp,int *nytp,int *nzyxt,
  float *xp,float *yp,float *zp,float *pxp,float *pyp,float *pzp,
  float *e1xp,float *e1yp,float *e1zp,float *e2xp,float *e2yp,float *e2zp,
  float *q111p,float *q112p,float *q121p,float *q122p,float *p211p,
  float *p212p,float *p221p,float *p222p,float *q211p,float *q212p,
  float *q221p,float *q222p,float *vp,float *dvdxp,float *dvdyp,
  float *dvdzp,int *nrsp,float *a0p,float *azh0p,int *n1,int *n2,float *p2p,
  float *q2p,float *hp,float *gradvp,float *d2tp,int *i2,int *i3,int *i6, 
  int *map,float *vs,float *dvdxs,float *dvdys,float *dvdzs,
  float *uxx,float *uxy,float *uxz,float *uyy,float *uyz,float *uzz,float *tzt,
  float *xx,float *yy,float *zz,float *pxs,float *pys,float *pzs,
  float *e1xs,float *e1ys,float *e1zs,float *e2xs,float *e2ys,float *e2zs,
  float *p111s,float *p112s,float *p121s,float *p122s,float *q111s,
  float *q112s,float *q121s,float *q122s,float *p211s,float *p212s,
  float *p221s,float *p222s,float *q211s,float *q212s,float *q221s,float *q222s,
  float *dxx,float *dyy,float *dzz,float *dpxs,float *dpys,float *dpzs,
  float *de1x,float *de1y,float *de1z,float *de2x,float *de2y,float *de2z,
  float *dp111,float *dp112,float *dp121,float *dp122,float *dq111,
  float *dq112,float *dq121,float *dq122,float *dp211,float *dp212,
  float *dp221,float *dp222,float *dq211,float *dq212,float *dq221,float *dq222,
  float *xt,float *yt,float *zt,float *pxt,float *pyt,float *pzt,
  float *e1xt,float *e1yt,float *e1zt,float *e2xt,float *e2yt,float *e2zt,
  float *p111t,float *p112t,float *p121t,float *p122t,float *q111t,
  float *q112t,float *q121t,float *q122t,float *p211t,float *p212t,
  float *p221t,float *p222t,float *q211t,float *q212t,float *q221t,
  float *q222t,float *dxt,float *dyt,float *dzt,float *dpxt,float *dpyt,
  float *dpzt,float *de1xt,float *de1yt,float *de1zt,float *de2xt,
  float *de2yt,float *de2zt,float *dp111t,float *dp112t,float *dp121t,
  float *dp122t,float *dq111t,float *dq112t,float *dq121t,float *dq122t,
  float *dp211t,float *dp212t,float *dp221t,float *dp222t,float *dq211t,
  float *dq212t,float *dq221t,float *dq222t);
				   

void ov2int_(float *v,int *nx,int *ny,int *nz,float *fx,float *fy,float *fz,
  float *dx,float *dy,float *dz,float *ov2,int *nxo,int *nyo,int *nzo,
  float *fxo,float *fyo,float *fzo,float *dxo,float *dyo,float *dzo);
 
void dv2_(int *nx,int *ny,int *nz,float *dx,float *dy,float *dz,float *v,
  float *vxx,float *vxy,float *vxz,float *vyy,float *vyz,float *vzz);
 
void timeb_(int *nr,int *nz,float *dr,float *dz,float *fz,float *a,
  float *v0,float *t);

void resit_(float *t,int *nx,int *ny,int *nz,int *nr,float *dx,float *dy,
  float *dr,float *fx,float *fy,float *x0,float *y0,float *tb);

void interp_(int *nx,int *ny,int *nz,int *nout,float *fx,float *fy,
  float *dx,float *dy,float *x,float *y,float *t,float *tout);

void recot_(float *t,int *nout,int *nz,int *nr,float *dr,float *x,float *y,
  float *x0,float *y0,float *tb);


main(int argc, char **argv)
{
	int 	na,nt,nxs,nys,nxo,nyo,nzo,nx,ny,nz;
	int 	nxot,nyot,nzot,ixot,iyot,ixo,iyo,itemp,ixt,iyt;
	int 	ix,iy,nxtmax,nytmax,nzyxt;
	float   dt,xs,ys,fxs,fys,dxs,dys,exs,eys,
		fxo,fyo,fzo,dxo,dyo,dzo,exo,eyo,
		fa,amin,amax,da,azhmin,azhmax,
		fx,fy,fz,dx,dy,dz,ex,ey,ez,
		fac,tmax,aperx,apery,temp,
		*v,*vxx,*vxy,*vxz,*vyy,*vyz,*vzz,
		*vt,*vxxt,*vxyt,*vxzt,*vyyt,*vyzt,*vzzt,		
 		*t,*s,*ov2,*tt1,*tt2;
	float	fxot,fyot,fzot,dxot,dyot,dzot,xomin,xomax,yomin,yomax,
		xbmin,xbmax,xemin,xemax,yemin,yemax,ybmin,ybmax;
	float 	*ysp,*xsp,*azhnp,*azhxp,*fxtp,*fytp;
	int	*nxtp,*nytp,nzyx,nzyxo;
	int ixs,iys,ixs0,iys0,nsize,ek,ms,ncpu,np,is,is0,isnow,ip;
	int maxno,nev,iev,iout,*nout;
	float *xbeg,*xend,*ybeg,*yend,*dxout,*dyout,*xout,*yout,*tout;
	float v0,dvz,*tb,dr,rmax;
	int nr;
 	char *vfile, *tfile, *jpfile;
	string *tofile;
 	char *restart; 	
	int isres;
 	FILE *vfp, *tfp,**tofp,*jpfp;

	float gmin, gmax;
	int itmp;
  
        usghed ugh;
	int ierr; 
/*	int dtype, n1, n2, n3, n4, n5;
	float d1,d2,d3,d4,d5,o4,o5,dcdp2,dline3,ocdp2,oline3;*/

	long long isize;
	char *envs;
  
  float *xp,*yp,*zp,*pxp,*pyp,*pzp,*e1xp,*e1yp,*e1zp,*e2xp,*e2yp,*e2zp;
  float *q111p,*q112p,*q121p,*q122p,*p211p;
  float *p212p,*p221p,*p222p,*q211p,*q212p;
  float *q221p,*q222p,*vp,*dvdxp,*dvdyp,*dvdzp;
  int *nrsp;
  float *a0p,*azh0p;
  int n1,n2;
  float *p2p,*q2p,*hp,*gradvp,*d2tp;
  int i2,i3,i6,*map;
  float *vs,*dvdxs,*dvdys,*dvdzs,*uxx,*uxy,*uxz,*uyy,*uyz,*uzz,*tzt;
  float *xx,*yy,*zz,*pxs,*pys,*pzs;
  float *e1xs,*e1ys,*e1zs,*e2xs,*e2ys,*e2zs;
  float *p111s,*p112s,*p121s,*p122s,*q111s;
  float *q112s,*q121s,*q122s,*p211s,*p212s;
  float *p221s,*p222s,*q211s,*q212s,*q221s,*q222s;
  float *dxx,*dyy,*dzz,*dpxs,*dpys,*dpzs;
  float *de1x,*de1y,*de1z,*de2x,*de2y,*de2z;
  float *dp111,*dp112,*dp121,*dp122,*dq111;
  float *dq112,*dq121,*dq122,*dp211,*dp212;
  float *dp221,*dp222,*dq211,*dq212,*dq221,*dq222;
  float *xt,*yt,*zt,*pxt,*pyt,*pzt;
  float *e1xt,*e1yt,*e1zt,*e2xt,*e2yt,*e2zt;
  float *p111t,*p112t,*p121t,*p122t,*q111t;
  float *q112t,*q121t,*q122t,*p211t,*p212t;
  float *p221t,*p222t,*q211t,*q212t,*q221t;
  float *q222t,*dxt,*dyt,*dzt,*dpxt,*dpyt;
  float *dpzt,*de1xt,*de1yt,*de1zt,*de2xt;
  float *de2yt,*de2zt,*dp111t,*dp112t,*dp121t;
  float *dp122t,*dq111t,*dq112t,*dq121t,*dq122t;
  float *dp211t,*dp212t,*dp221t,*dp222t,*dq211t;
  float *dq212t,*dq221t,*dq222t;
				   
	
	/* hook up getpar to handle the parameters */
	initargs(argc,argv);
	requestdoc(1);
		
	/* get velocity information from header file */
	if( !getparstring("vfile",&vfile) ) {
		vfp = stdin;
	} else {
		vfp = efopen(vfile,"r"); 
	} 
	ierr = fgetusghdr(vfp,&ugh);
	if(ierr!=0) err("fgetusghdr error");
	nx = ugh.n3;
	ny = ugh.n2;
	nz = ugh.n1;
	fx = ugh.o3;
	fy = ugh.o2;
	fz = ugh.o1;
	dx = ugh.d3;
	dy = ugh.d2;
	dz = ugh.d1;
	if(nx<3 || ny<3 || nz<3 )
err("number of velocity samples in each direction must be greater than 3!\n");

	ex = fx+(nx-1)*dx;
	ey = fy+(ny-1)*dy;
	ez = fz+(nz-1)*dz;

	/* get optional parameters */
	if (!getparstring("jpfile",&jpfile)) {
		jpfp = stderr;
	}else {
		jpfp = fopen(jpfile,"w");
	}
	if (!getparint("nt",&nt)) nt = 401;
	if (!getparint("ncpu",&ncpu)) ncpu = 1;
	/* use ncpu=1 now */
	/*
	if(ncpu>1) { ncpu = 1; fprintf(jpfp,"ncpu reset to 1 \n"); }
	*/
	envs = (char*) emalloc(80*sizeof(char));
    if(!getenv("PARALLEL")) {
                sprintf(envs,"%s=%d","PARALLEL",ncpu);
                putenv(envs);
				/* free(envs); */
    }
	if(nt>1001) err("nt cannot exceed 1001!\n");
	if (!getparfloat("dt",&dt)) dt = 0.008;
	tmax = (nt-1)*dt; 
	if (dt<0.000001) err("dt must be positive!\n");
 	if (!getparint("nyo",&nxot)) nxot = nx;
	if (!getparint("nxo",&nyot)) nyot = ny;
	if (!getparint("nzo",&nzot)) nzot = nz;
	if (!getparfloat("fyo",&fxot)) fxot = fx;
	if (!getparfloat("fxo",&fyot)) fyot = fy;
	if (!getparfloat("fzo",&fzot)) fzot = fz;
	if (!getparfloat("dyo",&dxot)) dxot = dx;
	if (!getparfloat("dxo",&dyot)) dyot = dy;
 	if (!getparfloat("dzo",&dzot)) dzot = dz;

	nzo = nzot;
	fzo = fzot;
	dzo = dzot;

	if (!getparint("ntline",&nev)) nev = 0;
	if(nev){
	    if (countparname("nout")!=nev)
		err("a nout array must be specified for each event");
	    if (countparname("xbeg")!=nev)
		err("a xbeg array must be specified for each event");
	    if (countparname("xend")!=nev)
		err("a xend array must be specified for each event");
	    if (countparname("ybeg")!=nev)
		err("a ybeg array must be specified for each event");
	    if (countparname("yend")!=nev)
		err("a yend array must be specified for each event");
	    if (countparname("tofile")!=nev)
		err("a tofile array must be specified for each event");

	    maxno = 0;
	    nout = alloc1int(nev);
	    dxout = alloc1float(nev);
	    dyout = alloc1float(nev);
	    xbeg = alloc1float(nev);
	    xend = alloc1float(nev);
	    ybeg = alloc1float(nev);
	    yend = alloc1float(nev);
	    tofile = (string *) malloc(nev*sizeof(string));
	    tofp = (FILE**) malloc(sizeof(FILE *)*nev);
  	    for(iev=0; iev<nev; ++iev) {
		if(!getnparint(iev+1,"nout",&nout[iev])) 
			err("cannot get nout for event %i. \n",1+iev);
		if(maxno<nout[iev]) maxno = nout[iev];
		if(!getnparfloat(iev+1,"xbeg",&ybeg[iev])) 
			err("cannot get xbeg for event %i. \n",1+iev);
		if(!getnparfloat(iev+1,"xend",&yend[iev])) 
			err("cannot get xend for event %i. \n",1+iev);
		if(!getnparfloat(iev+1,"ybeg",&xbeg[iev])) 
			err("cannot get ybeg for event %i. \n",1+iev);
		if(!getnparfloat(iev+1,"yend",&xend[iev]))
			err("cannot get yend for event %i. \n",1+iev);
		if(xbeg[iev]<fx || xbeg[iev]>ex ||
		   xend[iev]<fx || xend[iev]>ex ||
		   ybeg[iev]<fy || ybeg[iev]>ey ||
 		   yend[iev]<fy || yend[iev]>ey)
		    	err("event %i is out of output range!\n",1+iev);
		dxout[iev] = (xend[iev]-xbeg[iev])/(nout[iev]-1);
		dyout[iev] = (yend[iev]-ybeg[iev])/(nout[iev]-1);
		if(!getnparstring(iev+1,"tofile",&tofile[iev]))
			err("cannot get output file for event %i. \n",1+iev);
 	    }
	    tout = alloc1float(maxno*nzo);
	    xout = alloc1float(maxno);
	    yout = alloc1float(maxno);
	}
/* compute x-y geometry of extraploation/eikonal time table */
	if(nev==0){
		temp = dxot/dx;
		ixot = NINT(temp);
		ixot = MAX(ixot,1);
		dxo = dxot/ixot;
		temp = dyot/dy;
		iyot = NINT(temp);
		iyot = MAX(iyot,1);
		dyo = dyot/iyot;
		nxo = (nxot-1)*ixot+1;
		nyo = (nyot-1)*iyot+1;
		fxo = fxot;
		fyo = fyot;
	} else {
		fminmax(xbeg,nev,&xbmin,&xbmax);		
		fminmax(xend,nev,&xemin,&xemax);		
		fminmax(ybeg,nev,&ybmin,&ybmax);		
		fminmax(yend,nev,&yemin,&yemax);		
		xomin = MIN(xbmin,xemin);
		xomax = MAX(xbmax,xemax);
		yomin = MIN(ybmin,yemin);
		yomax = MAX(ybmax,yemax);
		dxo = dx;
		dyo = dy;
		nxo = (xomax-xomin)/dxo+1.5;
		nyo = (yomax-yomin)/dyo+1.5;
		fxo = xomin;
		fyo = yomin;
/* added for nev>0  z.li 10/2/96 */
		fxot = xomin;
		fyot = yomin;
	} 

	if (nxo<3) {
		nxo = 3;	
		fxo = fxot - dxo;
	}  
	if (nyo<3) {
		nyo = 3; 
		fyo = fyot - dyo;
	}

	exo = fxo+(nxo-1)*dxo;
	eyo = fyo+(nyo-1)*dyo;
	

 	if (!getparint("nys",&nxs)) nxs = 1;
	if (!getparint("nxs",&nys)) nys = 1;
	if (!getparfloat("fys",&fxs)) fxs = fx;
	if (!getparfloat("fxs",&fys)) fys = fy;
	if (!getparfloat("dys",&dxs)) dxs = dxot * 2.;
	if (!getparfloat("dxs",&dys)) dys = dyot * 2.;
	exs = fxs+(nxs-1)*dxs;
	eys = fys+(nys-1)*dys;
	if (!getparfloat("apery",&aperx)) aperx = 0.5*nx*dx;
	if(nxs>1) aperx += dxs;
	if (!getparfloat("aperx",&apery)) apery = 0.5*ny*dy;
	if(nys>1) apery += dys;
	if (!getparfloat("fa",&fa)) fa = 0;
	if (fa<0) err("fa must be not negative!\n");
	if (!getparfloat("da",&da)) da = 1.;
	if (da<0.0001) err("da must be positive!\n");
	if (!getparint("na",&na)) na = 61;
	if (!getparfloat("amin",&amin)) amin = 0;
	if (!getparfloat("amax",&amax)) amax = 90;
	if (amax>180 || amin<0 ) 
		err("polar angle must be within 0 to 180 degrees!\n");	
	if (!getparfloat("azhmin",&azhmin)) azhmin = 0;
	if (!getparfloat("azhmax",&azhmax)) azhmax = 360;
	if (azhmax>360 || azhmin<0 ) 
		err("azimuth angle must be within 0 and 360 degrees!\n");	
 	fa *= PI/180.;
	da *= PI/180.;
	amin *= PI/180.;
	amax *= PI/180.;
	azhmin *= PI/180.;
	azhmax *= PI/180.;
  	if (!getparfloat("fac",&fac)) fac = 0.01;
	fac *= 1.4141;
	if (!getparint("ek",&ek)) ek = 1;
 	if (!getparstring("restart",&restart)) restart = "n"; 
 	if (!getparint("isres",&isres)) isres = -999; 
	if (!getparint("ms",&ms)) ms = 1;
 

	if(!nev){
 	    ixs0 = 0;
	    iys0 = 0;
	    isize = 0;

	    nsize = nzot*nyot*nxot;
	    if( !getparstring("tfile",&tfile) ) {
		tfp = stdout;
	    } else {
		if((tfp = fopen(tfile,"r"))!=NULL) {
			fclose(tfp);
			tfp = fopen(tfile,"r+");	
		} else {
	   		tfp = fopen(tfile,"w");
		}
		if(restart[0]=='y') { 
			fseek2g(tfp,0,SEEK_END);
			isize = ftell64(tfp);
			isize=isize/sizeof(float)/nsize;
			if(isize>isres && isres>0) isize=isres-1;
			isres = isize + 1;
			ixs0 = isize/nys;
			iys0 = isize-ixs0*nys;
			/* if(isize==ns) return 0; */
		} else {
			fclose(tfp);
			tfp = fopen(tfile,"w");	
			isres = 1; 
		}
		fseek2g(tfp,isize*nsize*sizeof(float),SEEK_SET);
	    }
	} else {
 	    ixs0 = nxs;
	    iys0 = nys;
	    for(iev=0; iev<nev; ++iev){
		ix = 0;
		iy = 0;
		isize = 0;
		nsize = nzo*nout[iev];

 		if((tofp[iev] = fopen(tofile[iev],"r"))!=NULL) {
			fclose(tofp[iev]);
			tofp[iev] = fopen(tofile[iev],"r+");	
		} else {
	   		tofp[iev] = fopen(tofile[iev],"w");
		}
		if(restart[0]=='y') { 
			fseek2g(tofp[iev],0,SEEK_END);
			isize = ftell64(tofp[iev]);
			isize = isize/sizeof(float)/nsize;
			ix = isize/nys;
			iy = isize-ix*nys;
			/* if(isize==ns) return 0; */
		} else {
			fclose(tofp[iev]);
			tofp[iev] = fopen(tofile[iev],"w");	
			isres = 1;
		}
	    	if(ixs0>ix) {
			ixs0=ix;
	    		iys0=iy;
		}
	    }
	    is0 = ixs0*nys + iys0;
	    if(is0>isres && isres>0) is0 = isres - 1;
	    isres = is0 + 1; 
	    ixs0 = is0/nys;
	    iys0 = is0-ixs0*nys;
	    for(iev=0; iev<nev; ++iev){
		nsize = nzo*nout[iev];
		isize = ixs0*nys+iys0;
		fseek2g(tofp[iev],isize*nsize*sizeof(float),0);
	    }
	}


 
	/* ensure sources is in grid */

	if(fx>fxs || ex<exs || fy>fys || ey<eys || fz>0) 
	err("source lies outside of specified v(x,y,z) grid\n");
  
	if(fx>fxo || ex<exo  || fy>fyo || ey<eyo
		|| fz>fzo || ez<fzo+(nzo-1)*dzo) { 
	warn("fxo=%g fx=%g exo=%g ex=%g \n",fyo,fy,eyo,ey);
	warn("fyo=%g fy=%g eyo=%g ey=%g \n",fxo,fx,exo,ex);
	warn("fzo=%g fz=%g ezo=%g ez=%g \n",fzo,fz,fzo+(nzo-1)*dzo,ez);
	err("extrapolation lies outside of specified v(x,y,z) grid \n");
	}
 
	/* compute nxtmax and nytmax for ***t arrays */
	temp = 2.*aperx/dx + 2.5;
	nxtmax = temp;
	if(nxtmax>nx) nxtmax = nx;
	temp = 2.*apery/dy + 2.5;
	nytmax = temp;
	if(nytmax>ny) nytmax = ny;
	nzyxt = nxtmax * nytmax * nz;
/* velocity input */
 	/* allocate space */

	fprintf(stderr,"nxtmax=%d nytmax=%d nz=%d \n",nxtmax,nytmax,nz);
	fprintf(stderr,"before ealloc v nx=%d ny=%d nz=%d \n",nx,ny,nz);
	fprintf(stderr,"nzo=%d nxo=%d nyo=%d \n",nzo,nxo,nyo);
	/*
 	v = ealloc1float(nx*ny*nz);
 	vxx = ealloc1float(nx*ny*nz);
	vxy = ealloc1float(nx*ny*nz);
	vxz = ealloc1float(nx*ny*nz);
	vyy = ealloc1float(nx*ny*nz);
	vyz = ealloc1float(nx*ny*nz);
	*/
 	v = emalloc(nx*ny*nz*sizeof(float));
 	vxx = emalloc(nx*ny*nz*sizeof(float));
 	vxy = emalloc(nx*ny*nz*sizeof(float));
 	vxz = emalloc(nx*ny*nz*sizeof(float));
 	vyy = emalloc(nx*ny*nz*sizeof(float));
 	vyz = emalloc(nx*ny*nz*sizeof(float));
	fprintf(stderr,"before ealloc vzz nzyxt=%d ncpu=%d\n",nzyxt,ncpu);
 	vzz = emalloc(nx*ny*nz*sizeof(float));
	/*
	vzz = ealloc1float(nx*ny*nz);
	*/
  	vt = ealloc1float(nzyxt*ncpu);
 	vxxt = ealloc1float(nzyxt*ncpu);
	vxyt = ealloc1float(nzyxt*ncpu);
	vxzt = ealloc1float(nzyxt*ncpu);
	vyyt = ealloc1float(nzyxt*ncpu);
	vyzt = ealloc1float(nzyxt*ncpu);
	vzzt = ealloc1float(nzyxt*ncpu);
 	ov2 = ealloc1float(nzo*nyo*nxo);

/*
	fprintf(stderr,"before read vfp nx=%d ny=%d nz=%d \n",nx,ny,nz);
*/
 
	/* read velocities */
	fseek(vfp,0,0);
	fread(v,sizeof(float),nx*ny*nz,vfp);
  
	/* compute second derivatives of velocity  */
	dv2_(&nx,&ny,&nz,&dx,&dy,&dz,v,vxx,vxy,vxz,vyy,vyz,vzz);
  
 	/* compute slowness squares  */
	ov2int_(v,&nx,&ny,&nz,&fx,&fy,&fz,&dx,&dy,&dz,ov2,&nxo,&nyo,&nzo,
		&fxo,&fyo,&fzo,&dxo,&dyo,&dzo);

	if(nev){
    	    if (!getparfloat("v0",&v0)) v0 = 5000;
    	    if (!getparfloat("dvz",&dvz)) dvz = 0.;
	    dr = (dxo<dyo)?dxo/3.:dyo/3.;
	    rmax = sqrt((exs-fxo)*(exs-fxo)+(eys-fyo)*(eys-fyo));
	    if(rmax<sqrt((exo-fxs)*(exo-fxs)+(eyo-fys)*(eyo-fys)))
		rmax = sqrt((exo-fxs)*(exo-fxs)+(eyo-fys)*(eyo-fys));
	    nr = 1.5+rmax/dr;
	    tb = alloc1float(nr*nzo);
 	    timeb_(&nr,&nzo,&dr,&dzo,&fzo,&dvz,&v0,tb);
 	}
   
	/* allocate time arrays   */
 	tt1 = ealloc1float(nyo*nxo*ncpu);
 	tt2 = ealloc1float(nyo*nxo*ncpu);
	t = ealloc1float(nzo*nyo*nxo*ncpu);
	s = ealloc1float(nzo*nyo*nxo*ncpu);
	ysp = ealloc1float(ncpu);
	xsp = ealloc1float(ncpu);
	azhnp = ealloc1float(ncpu);
	azhxp = ealloc1float(ncpu);
	fxtp = ealloc1float(ncpu);
	fytp = ealloc1float(ncpu);
	nxtp = ealloc1int(ncpu);
	nytp = ealloc1int(ncpu);

    i2 = 2;
	i3 = 3;
	i6 = 6;
	n1 = 128;
	n2 = 1001*n1; 
	/* allocate other working arrays */
	xp = ealloc1float(n2*ncpu);
	yp = ealloc1float(n2*ncpu);
	zp = ealloc1float(n2*ncpu);
	pxp = ealloc1float(n2*ncpu);
	pyp = ealloc1float(n2*ncpu);
	pzp = ealloc1float(n2*ncpu);
	e1xp = ealloc1float(n2*ncpu);
	e1yp = ealloc1float(n2*ncpu);
	e1zp = ealloc1float(n2*ncpu);
	e2xp = ealloc1float(n2*ncpu);
	e2yp = ealloc1float(n2*ncpu);
	e2zp = ealloc1float(n2*ncpu);
	q111p = ealloc1float(n2*ncpu);
	q112p = ealloc1float(n2*ncpu);
	q121p = ealloc1float(n2*ncpu);
	q122p = ealloc1float(n2*ncpu);
	p211p = ealloc1float(n2*ncpu);
	p212p = ealloc1float(n2*ncpu);
	p221p = ealloc1float(n2*ncpu);
	p222p = ealloc1float(n2*ncpu);
	q211p = ealloc1float(n2*ncpu);
	q212p = ealloc1float(n2*ncpu);
	q221p = ealloc1float(n2*ncpu);
	q222p = ealloc1float(n2*ncpu);
	vp = ealloc1float(n2*ncpu);
	dvdxp = ealloc1float(n2*ncpu);
	dvdyp = ealloc1float(n2*ncpu);
	dvdzp = ealloc1float(n2*ncpu);

	nrsp = ealloc1int(n1*ncpu);
	a0p = ealloc1float(n1*ncpu);
	azh0p = ealloc1float(n1*ncpu);
	p2p = ealloc1float(i2*i2*ncpu);
	q2p = ealloc1float(i2*i2*ncpu);
	hp = ealloc1float(i3*i3*ncpu);
	gradvp = ealloc1float(i3*ncpu);
	d2tp = ealloc1float(i6*ncpu);

	map = ealloc1int(n1*ncpu);
	vs = ealloc1float(n1*ncpu);
	dvdxs = ealloc1float(n1*ncpu);
	dvdys = ealloc1float(n1*ncpu);
	dvdzs = ealloc1float(n1*ncpu);
	uxx = ealloc1float(n1*ncpu);
	uxy = ealloc1float(n1*ncpu);
	uxz = ealloc1float(n1*ncpu);
	uyy = ealloc1float(n1*ncpu);
	uyz = ealloc1float(n1*ncpu);
	uzz = ealloc1float(n1*ncpu);
	tzt = ealloc1float(n1*ncpu);
	xx = ealloc1float(n1*ncpu);
	yy = ealloc1float(n1*ncpu);
	zz = ealloc1float(n1*ncpu);
	pxs = ealloc1float(n1*ncpu);
	pys = ealloc1float(n1*ncpu);
	pzs = ealloc1float(n1*ncpu);
	e1xs = ealloc1float(n1*ncpu);
	e1ys = ealloc1float(n1*ncpu);
	e1zs = ealloc1float(n1*ncpu);
	e2xs = ealloc1float(n1*ncpu);
	e2ys = ealloc1float(n1*ncpu);
	e2zs = ealloc1float(n1*ncpu);
	p111s = ealloc1float(n1*ncpu);
	p112s = ealloc1float(n1*ncpu);
	p121s = ealloc1float(n1*ncpu);
	p122s = ealloc1float(n1*ncpu);
	q111s = ealloc1float(n1*ncpu);
	q112s = ealloc1float(n1*ncpu);
	q121s = ealloc1float(n1*ncpu);
	q122s = ealloc1float(n1*ncpu);
	p211s = ealloc1float(n1*ncpu);
	p212s = ealloc1float(n1*ncpu);
	p221s = ealloc1float(n1*ncpu);
	p222s = ealloc1float(n1*ncpu);
	q211s = ealloc1float(n1*ncpu);
	q212s = ealloc1float(n1*ncpu);
	q221s = ealloc1float(n1*ncpu);
	q222s = ealloc1float(n1*ncpu);
	dxx = ealloc1float(n1*ncpu);
	dyy = ealloc1float(n1*ncpu);
	dzz = ealloc1float(n1*ncpu);
	dpxs = ealloc1float(n1*ncpu);
	dpys = ealloc1float(n1*ncpu);
	dpzs = ealloc1float(n1*ncpu);
	de1x = ealloc1float(n1*ncpu);
	de1y = ealloc1float(n1*ncpu);
	de1z = ealloc1float(n1*ncpu);
	de2x = ealloc1float(n1*ncpu);
	de2y = ealloc1float(n1*ncpu);
	de2z = ealloc1float(n1*ncpu);
	dp111 = ealloc1float(n1*ncpu);
	dp112 = ealloc1float(n1*ncpu);
	dp121 = ealloc1float(n1*ncpu);
	dp122 = ealloc1float(n1*ncpu);
	dq111 = ealloc1float(n1*ncpu);
	dq112 = ealloc1float(n1*ncpu);
	dq121 = ealloc1float(n1*ncpu);
	dq122 = ealloc1float(n1*ncpu);
	dp211 = ealloc1float(n1*ncpu);
	dp212 = ealloc1float(n1*ncpu);
	dp221 = ealloc1float(n1*ncpu);
	dp222 = ealloc1float(n1*ncpu);
	dq211 = ealloc1float(n1*ncpu);
	dq212 = ealloc1float(n1*ncpu);
	dq221 = ealloc1float(n1*ncpu);
	dq222 = ealloc1float(n1*ncpu);

	xt = ealloc1float(n1*ncpu);
	yt = ealloc1float(n1*ncpu);
	zt = ealloc1float(n1*ncpu);
	pxt = ealloc1float(n1*ncpu);
	pyt = ealloc1float(n1*ncpu);
	pzt = ealloc1float(n1*ncpu);
	e1xt = ealloc1float(n1*ncpu);
	e1yt = ealloc1float(n1*ncpu);
	e1zt = ealloc1float(n1*ncpu);
	e2xt = ealloc1float(n1*ncpu);
	e2yt = ealloc1float(n1*ncpu);
	e2zt = ealloc1float(n1*ncpu);
	p111t = ealloc1float(n1*ncpu);
	p112t = ealloc1float(n1*ncpu);
	p121t = ealloc1float(n1*ncpu);
	p122t = ealloc1float(n1*ncpu);
	q111t = ealloc1float(n1*ncpu);
	q112t = ealloc1float(n1*ncpu);
	q121t = ealloc1float(n1*ncpu);
	q122t = ealloc1float(n1*ncpu);
	p211t = ealloc1float(n1*ncpu);
	p212t = ealloc1float(n1*ncpu);
	p221t = ealloc1float(n1*ncpu);
	p222t = ealloc1float(n1*ncpu);
	q211t = ealloc1float(n1*ncpu);
	q212t = ealloc1float(n1*ncpu);
	q221t = ealloc1float(n1*ncpu);
	q222t = ealloc1float(n1*ncpu);
	dxt = ealloc1float(n1*ncpu);
	dyt = ealloc1float(n1*ncpu);
	dzt = ealloc1float(n1*ncpu);
	dpxt = ealloc1float(n1*ncpu);
	dpyt = ealloc1float(n1*ncpu);
	dpzt = ealloc1float(n1*ncpu);
	de1xt = ealloc1float(n1*ncpu);
	de1yt = ealloc1float(n1*ncpu);
	de1zt = ealloc1float(n1*ncpu);
	de2xt = ealloc1float(n1*ncpu);
	de2yt = ealloc1float(n1*ncpu);
	de2zt = ealloc1float(n1*ncpu);
	dp111t = ealloc1float(n1*ncpu);
	dp112t = ealloc1float(n1*ncpu);
	dp121t = ealloc1float(n1*ncpu);
	dp122t = ealloc1float(n1*ncpu);
	dq111t = ealloc1float(n1*ncpu);
	dq112t = ealloc1float(n1*ncpu);
	dq121t = ealloc1float(n1*ncpu);
	dq122t = ealloc1float(n1*ncpu);
	dp211t = ealloc1float(n1*ncpu);
	dp212t = ealloc1float(n1*ncpu);
	dp221t = ealloc1float(n1*ncpu);
	dp222t = ealloc1float(n1*ncpu);
	dq211t = ealloc1float(n1*ncpu);
	dq212t = ealloc1float(n1*ncpu);
	dq221t = ealloc1float(n1*ncpu);
	dq222t = ealloc1float(n1*ncpu);

   	fprintf(jpfp," \n");
   	fprintf(jpfp," RAYT3D parameters \n");
   	fprintf(jpfp,"===================\n");
   	fprintf(jpfp," vfile=%s \n",vfile);
   	fprintf(jpfp," tfile=%s \n",tfile);
   	fprintf(jpfp," one-way time: dt=%g nt=%d \n",dt,nt);
   	fprintf(jpfp," nz=%d nx=%d ny=%d \n",nz,ny,nx);
   	fprintf(jpfp," fz=%g fx=%g fy=%g \n",fz,fy,fx);
   	fprintf(jpfp," dz=%g dx=%g dy=%g \n",dz,dy,dx);
   	fprintf(jpfp," nzo=%d nxo=%d nyo=%d \n",nzot,nyot,nxot);
   	fprintf(jpfp," fzo=%g fxo=%g fyo=%g \n",fzot,fyot,fxot);
   	fprintf(jpfp," dzo=%g dxo=%g dyo=%g \n",dzot,dyot,dxot);
   	fprintf(jpfp," nzoe=%d nxoe=%d nyoe=%d \n",nzo,nyo,nxo);
   	fprintf(jpfp," fzoe=%g fxoe=%g fyoe=%g \n",fzo,fyo,fxo);
   	fprintf(jpfp," dzoe=%g dxoe=%g dyoe=%g \n",dzo,dyo,dxo);
	if(nev){
   	fprintf(jpfp," \n");
   	fprintf(jpfp," ntline=%d v0=%g dvz=%g \n",nev,v0,dvz);
	for(iev=0; iev<nev; ++iev) {
   	fprintf(jpfp," \n");
   	fprintf(jpfp," itline=%d tofile=%s \n",iev+1,tofile[iev]);
   	fprintf(jpfp," xbeg=%g xend=%g ybeg=%g yend=%g \n",
		ybeg[iev],yend[iev],xbeg[iev],xend[iev]);
   	fprintf(jpfp," nout=%d dxout=%g dyout=%g \n",
		nout[iev],dyout[iev],dxout[iev]);
	}
	}
   	fprintf(jpfp," \n");
	fprintf(jpfp," nxs=%d nys=%d fxs=%g fys=%g dxs=%g dys=%g \n",
		nys,nxs,fys,fxs,dys,dxs);			  
	fprintf(jpfp," aperx=%g apery=%g \n",apery,aperx);
	fprintf(jpfp," fa=%g da=%g na=%d \n",fa*180./PI,da*180./PI,na);
	fprintf(jpfp," amin=%g amax=%g azhmin=%g azhmax=%g \n",
		amin*180./PI,amax*180./PI,azhmin*180./PI,azhmax*180./PI);
	fprintf(jpfp," ek=%d ms=%d ncpu=%d jpfile=%s restart=%s isres=%d\n",
		ek,ms,ncpu,jpfile,restart,isres);
	

   	fprintf(jpfp," \n");
   	fprintf(jpfp," finish velocity input ... \n");
   	fprintf(jpfp," begin traveltime calculation ... \n");
	fflush(jpfp);


	is0 = ixs0*nys + iys0;
	nzyx = nz * ny * nx;
	nzyxo = nzo * nyo * nxo;

	/* loop over sources */
	for (is=is0;is<nxs*nys;is=is+ncpu) {
            np = ncpu;
	    if(is+np>nxs*nys) np = nxs*nys - is;
	    for(ip=0;ip<np;ip++) {
		ixs = (is+ip)/nys;
		iys = (is+ip) - ixs*nys; 
		ysp[ip] = fys + iys*dys;		 
		xsp[ip] = fxs + ixs*dxs;		 
	    }
		
/*
	fprintf(stderr,"nzyx=%d nzyxo=%d nz=%d ny=%d nx=%d \n",
		nzyx,nzyxo,nz,ny,nx);
	fprintf(stderr,"nzo=%d nyo=%d nxo=%d ek=%d nt=%d dt=%g tmax=%g\n",
		nzo,nyo,nxo,ek,nt,dt,tmax);
fprintf(stderr,"na=%d fa=%g da=%g amin=%g amax=%g azhmin=%g azhmax=%g \n",
		na,fa,da,amin,amax,azhmin,azhmax);
fprintf(stderr,"fxo=%g fyo=%g fzo=%g dxo=%g dyo=%g dzo=%g fac=%g \n",
		fxo,fyo,fzo,dxo,dyo,dzo,fac);
fprintf(stderr,"fx=%g fy=%g fz=%g dx=%g dy=%g dz=%g aperx=%g apery=%g \n",
		fx,fy,fz,dx,dy,dz,aperx,apery);
		for(ip=0;ip<np;ip++) 
			fprintf(stderr,"ysp=%g xsp=%g \n",xsp[ip],ysp[ip]);
*/


	    /* parallel computation */
    rt3dp_(&np,&nzyx,&nzyxo,&nz,&ny,&nx,&nzo,&nyo,&nxo,&ek,&nt,
       &na,&fa,&da,&amin,&amax,&azhmin,&azhmax,&dt,&tmax,&fxo,&fyo,&fzo,
       &dxo,&dyo,&dzo,&fac,&fx,&fy,&fz,&dx,&dy,&dz,&aperx,&apery,v,vxx,vxy,vxz,
	   vyy,vyz,vzz,ov2,vt,vxxt,vxyt,vxzt,vyyt,vyzt,vzzt,tt1,tt2,t,s,
       xsp,ysp,azhnp,azhxp,fxtp,fytp,nxtp,nytp,&nzyxt,
	   xp,yp,zp,pxp,pyp,pzp,e1xp,e1yp,e1zp,e2xp,e2yp,e2zp,
	   q111p,q112p,q121p,q122p,p211p,p212p,p221p,p222p,
	   q211p,q212p,q221p,q222p,vp,dvdxp,dvdyp,dvdzp,nrsp,
	   a0p,azh0p,&n1,&n2,p2p,q2p,hp,gradvp,d2tp,&i2,&i3,&i6,
	   map,vs,dvdxs,dvdys,dvdzs,uxx,uxy,uxz,uyy,uyz,uzz,tzt,
	   xx,yy,zz,pxs,pys,pzs,e1xs,e1ys,e1zs,e2xs,e2ys,e2zs,
	   p111s,p112s,p121s,p122s,q111s,q112s,q121s,q122s,
   	   p211s,p212s,p221s,p222s,q211s,q212s,q221s,q222s,
	   dxx,dyy,dzz,dpxs,dpys,dpzs,de1x,de1y,de1z,de2x,
	   de2y,de2z,dp111,dp112,dp121,dp122,dq111,dq112,dq121,dq122,
	   dp211,dp212,dp221,dp222,dq211,dq212,dq221,dq222,
	   xt,yt,zt,pxt,pyt,pzt,e1xt,e1yt,e1zt,e2xt,e2yt,e2zt,
	   p111t,p112t,p121t,p122t,q111t,q112t,q121t,q122t,p211t,
	   p212t,p221t,p222t,q211t,q212t,q221t,q222t,dxt,dyt,
	   dzt,dpxt,dpyt,dpzt,de1xt,de1yt,de1zt,de2xt,de2yt,de2zt,
	   dp111t,dp112t,dp121t,dp122t,dq111t,dq112t,dq121t,dq122t,
	   dp211t,dp212t,dp221t,dp222t,dq211t,dq212t,dq221t,dq222t);
				   
	   /* check the travel time calculation */
	   for(ip=0;ip<np;ip++) {
		isnow = is + ip;
		itmp = nxo * nyo * nzo;
	   	fminmax(t+ip*itmp,itmp, &gmin, &gmax);
	   	/* print warning message if needed */
		if(gmin > 9999.) {
 		fprintf(jpfp,
		" NULL travel time table at source is=%d xs=%g ys=%g\n",
		isnow+1,ysp[ip],xsp[ip]);
   	      	fflush(jpfp);
 		warn(
		" NULL travel time table at source is=%d xs=%g ys=%g\n",
		isnow+1,ysp[ip],xsp[ip]);
		}
	    }


	    /* output travel time table */
	    if(!nev){ 
	    	/* write traveltime  	*/
	        for(ip=0;ip<np;ip++) {
		for(ixt=0;ixt<nxot;ixt++) {
			ixo = ixt*ixot;
			if(nxot<3) ixo +=1;
			for(iyt=0;iyt<nyot;iyt++) {
				iyo = iyt*iyot;
				if(nyot<3) iyo +=1;
				itemp = (ixo*nyo+iyo)*nzo;
  	    		fwrite(t+itemp+ip*nzyxo,sizeof(float),nzo,tfp);
			}
		}
		}
 	    } else { 
	        for(ip=0;ip<np;ip++) {
		xs = xsp[ip];
		ys = ysp[ip]; 
 		resit_(t+ip*nzyxo,&nyo,&nxo,&nzo,&nr,&dyo,&dxo,&dr,&fyo,&fxo,
			&ys,&xs,tb);
 	    	for(iev=0; iev<nev; ++iev){

		    for(iout=0; iout<nout[iev]; ++iout){
			xout[iout] = xbeg[iev]+iout*dxout[iev];
			yout[iout] = ybeg[iev]+iout*dyout[iev];
		    }
 
		    interp_(&nyo,&nxo,&nzo,&nout[iev],&fyo,&fxo,&dyo,&dxo,
		    	yout,xout,t+ip*nzyxo,tout);

  		    recot_(tout,&nout[iev],&nzo,&nr,&dr,yout,xout,&ys,&xs,tb);
 
   	    	    fwrite(tout,sizeof(float),nout[iev]*nzo,tofp[iev]);

/*
		itmp = nout[iev]*nzo;
		fminmax(tout,itmp, &gmin, &gmax);
		fprintf(stderr," ixs=%d gmin=%g gmax=%g ip=%d \n",
			ixs,gmin,gmax,ip);
*/

	        }
	        }
	    }

	    for(ip=0;ip<np;ip++) {
		isnow = is + ip;
	    if((isnow+1)%ms==0 && (isnow+1)>=ms ) {
fprintf(jpfp," travel time computed at source is=%d xs=%g ys=%g\n",
		isnow+1,ysp[ip],xsp[ip]);
   	      fflush(jpfp);
	    }
	    }
	}

	fprintf(jpfp," RAYT3D job done\n");

	ugh.scale = 1.e-6;
	ugh.dtype = 4;
 	ugh.n1 = nzot;
	ugh.n2 = nyot;
 	ugh.n3 = nxot;
 	ugh.n4 = nys;
	ugh.n5 = nxs;
 	ugh.d1 = dzot;
	ugh.d2 = dyot;
 	ugh.d3 = dxot;
 	ugh.d4 = dys;
	ugh.d5 = dxs;
	ugh.o1 = fzot; 
	ugh.o2 = fyot; 
	ugh.o3 = fxot; 
	ugh.o4 = fys;
	ugh.o5 = fxs;
	ugh.ocdp2 = 0.;
	ugh.oline3 = 0.;
	ugh.dcdp2 = 0.;
	ugh.dline3 = 0.;
	ugh.gmin = 0.;
	ugh.gmax = tmax;
	ugh.orient = 1;
	ugh.gtype = 0;


	for(iev=0; iev<nev; ++iev){
	    ugh.o2 = ybeg[iev];
	    ugh.o3 = xbeg[iev];

	    if(dxout[iev]==0){
		ugh.n3 = 1;
		ugh.n2 = nout[iev];
		ugh.d2 = dyout[iev];
		ugh.d3 = dxo;
	    } else if(dyout[iev]==0){
		ugh.n2 = 1;
		ugh.n3 = nout[iev];
		ugh.d3 = dxout[iev];
		ugh.d2 = dyo;
 	    } else {
		ugh.n2 = nout[iev];
		ugh.n3 = 1;
		ugh.d2 = sqrt(dxout[iev]*dxout[iev]+dyout[iev]*dyout[iev]);
		ugh.d3 = 999999.;
 	    }
 	    fputusghdr(tofp[iev],&ugh); 
 	}
 		 

 	if(!nev) {
		fputusghdr(tfp,&ugh); 
 		fclose(tfp);
	} else {
		for(iev=0; iev<nev; ++iev) fclose(tofp[iev]);
	}
 
 	/* free space */
 	free1float(v);
  	free1float(vxx);
 	free1float(vxy);
 	free1float(vxz);
 	free1float(vyy);
 	free1float(vyz);
 	free1float(vzz);
 	free1float(vt);
   	free1float(vxxt);
 	free1float(vxyt);
 	free1float(vxzt);
 	free1float(vyyt);
 	free1float(vyzt);
 	free1float(vzzt);
  	free1float(t);
	free1float(s);
	free1float(xsp);
	free1float(ysp);
	free1float(azhnp);
	free1float(azhxp);
	free1float(fxtp);
	free1float(fytp);
	free1int(nxtp);
	free1int(nytp);
	if(nev){
		free1int(nout);
		free1float(dxout);
		free1float(dyout);
		free1float(xbeg);
		free1float(ybeg);
		free1float(xend);
		free1float(xend);
		free1float(tout);
	}

	free1float(xp);
	free1float(yp);
	free1float(zp);
	free1float(pxp);
	free1float(pyp);
	free1float(pzp);
	free1float(e1xp);
	free1float(e1yp);
	free1float(e1zp);
	free1float(e2xp);
	free1float(e2yp);
	free1float(e2zp);
	free1float(q111p);
	free1float(q112p);
	free1float(q121p);
	free1float(q122p);
	free1float(p211p);
	free1float(p212p);
	free1float(p221p);
	free1float(p222p);
	free1float(q211p);
	free1float(q212p);
	free1float(q221p);
	free1float(q222p);
	free1float(vp);
	free1float(dvdxp);
	free1float(dvdyp);
	free1float(dvdzp);

	free1int(nrsp);
	free1float(a0p);
	free1float(azh0p);
	free1float(p2p);
	free1float(q2p);
	free1float(hp);
	free1float(gradvp);
	free1float(d2tp);

	free1int(map);
	free1float(vs);
	free1float(dvdxs);
	free1float(dvdys);
	free1float(dvdzs);
	free1float(uxx);
	free1float(uxy);
	free1float(uxz);
	free1float(uyy);
	free1float(uyz);
	free1float(uzz);
	free1float(tzt);
	free1float(xx);
	free1float(yy);
	free1float(zz);
	free1float(pxs);
	free1float(pys);
	free1float(pzs);
	free1float(e1xs);
	free1float(e1ys);
	free1float(e1zs);
	free1float(e2xs);
	free1float(e2ys);
	free1float(e2zs);
	free1float(p111s);
	free1float(p112s);
	free1float(p121s);
	free1float(p122s);
	free1float(q111s);
	free1float(q112s);
	free1float(q121s);
	free1float(q122s);
	free1float(p211s);
	free1float(p212s);
	free1float(p221s);
	free1float(p222s);
	free1float(q211s);
	free1float(q212s);
	free1float(q221s);
	free1float(q222s);
	free1float(dxx);
	free1float(dyy);
	free1float(dzz);
	free1float(dpxs);
	free1float(dpys);
	free1float(dpzs);
	free1float(de1x);
	free1float(de1y);
	free1float(de1z);
	free1float(de2x);
	free1float(de2y);
	free1float(de2z);
	free1float(dp111);
	free1float(dp112);
	free1float(dp121);
	free1float(dp122);
	free1float(dq111);
	free1float(dq112);
	free1float(dq121);
	free1float(dq122);
	free1float(dp211);
	free1float(dp212);
	free1float(dp221);
	free1float(dp222);
	free1float(dq211);
	free1float(dq212);
	free1float(dq221);
	free1float(dq222);

	free1float(xt);
	free1float(yt);
	free1float(zt);
	free1float(pxt);
	free1float(pyt);
	free1float(pzt);
	free1float(e1xt);
	free1float(e1yt);
	free1float(e1zt);
	free1float(e2xt);
	free1float(e2yt);
	free1float(e2zt);
	free1float(p111t);
	free1float(p112t);
	free1float(p121t);
	free1float(p122t);
	free1float(q111t);
	free1float(q112t);
	free1float(q121t);
	free1float(q122t);
	free1float(p211t);
	free1float(p212t);
	free1float(p221t);
	free1float(p222t);
	free1float(q211t);
	free1float(q212t);
	free1float(q221t);
	free1float(q222t);
	free1float(dxt);
	free1float(dyt);
	free1float(dzt);
	free1float(dpxt);
	free1float(dpyt);
	free1float(dpzt);
	free1float(de1xt);
	free1float(de1yt);
	free1float(de1zt);
	free1float(de2xt);
	free1float(de2yt);
	free1float(de2zt);
	free1float(dp111t);
	free1float(dp112t);
	free1float(dp121t);
	free1float(dp122t);
	free1float(dq111t);
	free1float(dq112t);
	free1float(dq121t);
	free1float(dq122t);
	free1float(dp211t);
	free1float(dp212t);
	free1float(dp221t);
	free1float(dp222t);
	free1float(dq211t);
	free1float(dq212t);
	free1float(dq221t);
	free1float(dq222t);

	return 0;
}
