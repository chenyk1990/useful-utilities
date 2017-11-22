function [phs,delay,cc,obj,kx,ky]=phaseslice(seis,t,x,y,xref,yref,t1s,twin,delx,dely)
% PHASESLICE: Compute a relative phase slice across a 3D seismic volume
%
% [phs,delay,cc,obj,kx,ky]=phaseslice(seis,t,x,y,xref,yref,t1s,twin,delx,dely)
%
% Method: A reference trace is chosen from the seismic volume and compared to traces at chosen
% locations. At each location, a local ensemble of traces, restricted in time to a chosen window, is
% compared to the reference trace (in the window) to determine the best time shift (delay) and
% constant phase rotation required to match the ensemble to the reference trace. The comparison is
% done by calling maxcorr_ephs . Phase measurments are made at every delx,dely locations beginning
% at 1,1 and extending to length(x) and length(y). For a particular location, the ensemble uses a
% small rectangular patch centered at the location and of size 2*delx+1 by 2*dely+1 of
% (2*delx+1)*(2*dely+1) traces. Thus the phase and delay are least-squares estimates for the
% ensemble and this tends to given spatially smooth estimates.  The returned phs and delay are
% interpolated (spline) to the input x,ylocations.
%
% seis ... seismic matrix as a 3D volume. t is the first dimension, x the second, and y the third
% t ... time coordinate for seis
% x ... x (dimension 2) coordinate for seis
% y ... y (dimension 3) coordinate for seis
% xref ... x location of reference trace
% yref ... y location of reference trace
% t1s ... matrix of times of the center of the analysis window. Must be either
%       length(x)-by-length(y) or a single scalar (for a horizontal window)
% NOTE: It is best to specify t1s over an area the encompases the seismic volume. For example,
% choose a rectangle that includes all the seismic traces and specify t1s at the four corners of
% this rectangle.
% twin ... time width of the phase analysis window in seconds. Window is centered on t1s and is of
%       the same width all along the profile. Don't make this too small. .5 sec or 1 second are good
%       choices.
% delx ... measurments will be made every this many x locations. Should be an integer>=1.
% ************* default =10 *************
% dely ... measurments will be made every this many y locations. Should be an integer>=1.
% ************* default =10 *************
%
% phs ... measured phases, one for each x,y
% delay ... measured delays, one for each x,y
% cc ... max cc after phase and delay correction, one for each x,y
% obj ... the phase objective function. This is a 3D matrix with the 2 and 3 dimensions being the
%       x-y measurment grid and the number of rows is the length of theta=-180:179. The rows are
%       thus phases at 1 degree increments. The estimated phase is the minimum of obj in each
%       column.
% kx ... vector of indices of x locations where phase was measured. Phase measurments were at x(kx)
% ky ... vector of indices of y locations where phase was measured. Phase measurments were at y(ky)
%
% The phase at location k is that required in phsrot to rotate location k to look like the reference
% That is, let sk be the trace at location k and sr be the reference trace, the
% skr=phsrot(sk,phs(k)) will look similar to sr except for a time shift. 
%

if(nargin<10)
    error('not enough input parameters');
end
if(nargin<10)
    dely=10;
end
if(nargin<9)
    delx=10;
end

[nt,nx,ny]=size(seis);
if(length(t)~=nt)
    error('t and seis have incompatable sizes');
end

if(length(x)~=nx)
    error('x and seis have incompatable sizes');
end

if(length(y)~=ny)
    error('y and seis have incompatable sizes');
end

if(length(t1s)==1)
    t1s=t1s*ones(nx,ny);
end

if(size(t1s,1)~=nx)
    error('row dimension of t1s must equal length(x)'); 
end
if(size(t1s,2)~=ny)
    error('column dimension of t1s must equal length(y)'); 
end


delx=round(delx);
if(delx<1 || delx>.25*nx)
    error('unreasonable value for delx');
end
dely=round(dely);
if(dely<1 || dely>.25*ny)
    error('unreasonable value for dely');
end

if(~between(x(1),x(end),xref,2))
    error('xref must lie between first and last x')
end
if(~between(y(1),y(end),yref,2))
    error('yref must lie between first and last y')
end


%adjust twin if needed
tmin=min(t1s);
tmax=max(t1s);
if(tmin-.505*twin<t(1))
    twin=(tmin-t(1))*1.95;
end
if(tmax+.505*twin>t(end))
    twin=(t(end)-tmax)*1.95;
end

%build reference vol
dt=t(2)-t(1);
iwin2=round(twin/(2*dt));
iref=near(x,xref);
jref=near(y,yref);
kref=near(t,t1s(iref(1),jref(1)));

sref=seis(kref(1)-iwin2:kref(1)+iwin2,iref(1),jref(1));

nlag=round(iwin2/2);

kx=1:delx:nx;
ky=1:dely:ny;
% ky=dely+1;
if(kx(end)~=nx)
    kx=[kx nx];
end
if(ky(end)~=ny)
    ky=[ky ny];
end

%find an inc
inc=1;
test=-inc*floor(dely/inc):inc:inc*floor(dely/inc);
while(length(test)>7)
    inc=inc+1;
    test=-inc*floor(dely/inc):inc:inc*floor(dely/inc);
end
iymax=inc*floor(dely/inc);
ixmax=inc*floor(delx/inc);

phsx=zeros(length(kx),length(ky));
delayx=phsx;
ccx=phsx;
my=0;
theta=-180:179;
objx=zeros(length(theta),length(kx),length(ky));
t0=clock;
for k=ky
    my=my+1;
    mx=0;
    for j=kx
        mx=mx+1;
%         kk=max([k-dely, 1]):min([k+dely,ny]);
%         jj=max([j-delx, 1]):min([j+delx,nx]);
        kk=max([k-iymax,1]):inc:min([k+iymax,ny]);
        jj=max([j-ixmax,1]):inc:min([j+ixmax,nx]);
        jt=near(t,t1s(k));
        s=seis(jt-iwin2:jt+iwin2,jj,kk);
        sr=sref(:,ones(1,length(jj)),ones(1,length(kk)));
        %xx=maxcorr_ephs(s(:),sr(:),nlag);
        [xx,~,~,tmp]=maxcorr_ephs2(s(:),sr(:),nlag);
        objx(:,mx,my)=flipud(tmp');%flipud to be consistent with next line
        phsx(mx,my)=-xx(3);%minus sign for consistency with previous
        ccx(mx,my)=xx(4);
        delayx(mx,my)=-xx(2)*dt;
    end
    tnow=clock;
    timeused=etime(tnow,t0);
    timeleft=(length(ky)-my)*timeused/my;
    disp(['Completed inline ' int2str(my) ' of ' int2str(length(ky))])
    disp(['time used= ' num2str(timeused/60) ', time left= ' num2str(timeleft/60)])
end


xm=x(kx);
ym=y(ky);
[YM,XM]=meshgrid(ym,xm);
[Y,X]=meshgrid(y,x);
phs=interp2(YM,XM,phsx,Y,X,'spline');
delay=interp2(YM,XM,delayx,Y,X,'spline');
cc=interp2(YM,XM,ccx,Y,X,'spline');
om=max(abs(objx(:)));
obj=zeros(length(theta),length(x),length(y));

for k=1:length(theta)
    obj(k,:,:)=interp2(YM,XM,squeeze(objx(k,:,:))/om,Y,X,'spline');
end




