function [phs,delay,cc,obj,kx]=phaseprofile(seis,t,x,xref,t1s,xs,twin,delx)
% PHASEPROFILE: Compute a relative phase profile across a seismic gather
%
% [phs,delay,c,obj,kx]=phaseprofile(seis,t,x,xref,t1s,xs,twin,delx)
%
% Method: A reference trace is chosen from the seismic matrix and compared to traces at chosen
% locations along the x axis. At each location, a local ensemble of traces, restricted in time to a
% chosen window, is compared to the reference trace (in the window) to determine the best time shift
% (delay) and constant phase rotation required to match the ensemble to the reference trace. The
% comparison is done by calling maxcorr_ephs . Phase measurments are made at every delx locations
% beginning at 1 and extending to length(x). For a particular location, the ensemble takes delx
% traces to the left and delx traces to the right for a total of 2*delx+1 traces. Thus the phase and
% delay are least-squares estimates for the ensemble and this tends to give spatially smooth
% estimates.  The returned phs, delay, and cc are interpolated (spline) to the input x locations.
%
% seis ... seismic matrix
% t ... time coordinate for seis
% x ... space coordinate for seis
% xref ... location of reference trace
% t1s ... vector of times of the center of the analysis window
% xs ... vector of x coordinates at which the t1s are specified. Must be the same size as t1s.
% twin ... time width of the phase analysis window in seconds. Window is centered on t1s and is of
%       the same width all along the profile. Don't make this too small. .5 sec or 1 second are good
%       choices.
% delx ... measurments will be made every this many x locations. Should be an integer>=1. Local
%       trace ensembles will contain 2*delx+1 traces at each location except for the ends. This is
%       specified in traces, not in x units.
% ************* default =10 *************
%
% phs ... measured phases, one for each x
% delay ... measured delays, one for each x
% cc ... max cc after phase and delay correction, one for each x
% obj ... the phase objective function. This is a matrix with length(x) columns and the number of
%       rows is the length of theta=-180:179. The rows are thus phases at 1 degree increments. The
%       estimated phase is the minimum of obj in each column.
% kx ... vector of indices of x locations where phase was measured. Phase measurments were at x(kx)
%
% The phase at location k is that required in phsrot to rotate location k to look like the reference
% That is, let sk be the trace at location k and sr be the reference trace, then
% skr=phsrot(sk,phs(k)) will look similar to sr except for a time shift. 
%

if(nargin<7)
    error('not enough input parameters');
end

if(nargin<8)
    delx=10;%specified in traces
end

if(length(t1s)~=length(xs))
    error('t1s and xs must have the same length'); 
end

[nt,nx]=size(seis);
if(length(t)~=nt)
    error('t and seis have incompatable sizes');
end

if(length(x)~=nx)
    error('x and seis have incompatable sizes');
end

delx=round(delx);
if(delx<1 || delx>.25*nx)
    error('unreasonable value for delx');
end

if(~between(x(1),x(end),xref,2))
    error('xref must lie between first and last x')
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

t1x=interpextrap(xs,t1s,x,0);

dt=t(2)-t(1);
iwin2=round(twin/(2*dt));
iref=near(x,xref);
jref=near(t,t1x(iref(1)));
sref=seis(jref(1)-iwin2:jref(1)+iwin2,iref(1));

nlag=round(iwin2/2);
% kx1=delx+1;
% kxn=nx-delx;
% kx=kx1:incx:kxn;
kx=1:delx:nx;
if(kx(end)~=nx)
    kx=[kx nx];
end
phsx=zeros(size(kx));
delayx=phsx;
ccx=phsx;
n=0;
theta=-180:179;
objx=zeros(length(theta),length(kx));
for k=kx
    n=n+1;
    kk=max([k-delx, 1]):min([k+delx,nx]);
    jx=near(t,t1x(k));
    s=seis(jx-iwin2:jx+iwin2,kk);
    sr=sref(:,ones(size(kk)));
    %xx=maxcorr_ephs(s(:),sr(:),nlag);
    [xx,~,~,tmp]=maxcorr_ephs2(s(:),sr(:),nlag);
%     [xx,tmp]=maxcorr_ephspro(s,sref,nlag);
    objx(:,n)=flipud(tmp');%flipud to be consistent with next line
    phsx(n)=-xx(3);%minus sign for consistency with previous
    ccx(n)=xx(4);
    delayx(n)=-xx(2)*dt;
end

% phs=interpextrap(x(kx),phsx,x,0);
% delay=interpextrap(x(kx),delayx,x,0);
% cc=interpextrap(x(kx),ccx,x,0);
xm=x(kx);
phs=interp1(xm,phsx,x,'spline');
delay=interp1(xm,delayx,x,'spline');
cc=interp1(xm,ccx,x,'spline');
om=max(abs(objx(:)));
obj=zeros(length(theta),length(x));

for k=1:length(theta)
    obj(k,:)=interp1(xm,objx(k,:)/om,x,'spline');
end




