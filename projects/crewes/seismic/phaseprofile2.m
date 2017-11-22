function [phs,delay,cc]=phaseprofile(seis,t,x,xref,t1s,xs,twin,delx)
% PHASEPROFILE: Compute a relative phase profile across a seismic gather
% [phs,delay]=phaseprofile(seis,t,x,xref,t1s,xs,twin)
%
% Method: A reference trace is chosen from the seismic matrix and compared to traces at chosen
% locations along the x axis. At each location, a local ensemble of traces, restricted in time to a
% chosen window, is compared to the reference trace (in the window) to determine the best time shift
% (delay) and constant phase rotation required to match the ensemble to the reference trace. The
% comparison is done by calling maxcorr_ephs . Phase measurments are made at every delx locations
% beginning at 1 and extending to length(x). For a particular location, the ensemble taks delx
% traces to the lect and delx traces to the right for a total of 2*delx+1 traces. Thus the phase and
% delay are least--squares estimates for the ensemble and this tensd to given spatially smooth
% estimates.  The returned phs and delay are interpolated (spline) to the input x locations.
%
% seis ... seismic matrix
% t ... time coordinate for seis
% x ... space coordinate for seis
% xref ... location of reference trace
% t1s ... vector of times of the center of the analysis window
% xs ... vector of x coordinates at which the t1s are specified. Must be the same size as t1s.
% twin ... time width of the phase analysis window. Window is centered on t1s and is of the same
%       width all along the profile.
% delx ... measurments will be made every this many x locations. Should be in integer>=1.
% ************* default =10 *************
%
% phs ... measured phases, one for each x
% delay ... measured delays, one for each x
% cc ... max cc after phase and delay correction
%

if(nargin<7)
    error('not enough input parameters');
end

if(nargin<8)
    delx=10;
end

if(length(t1s)~=length(xs))
    error('t1s and xs must have the same length'); 
end

[nt,nx]=size(seis);
if(length(t)~=nt)
    error('t and seis have incompaible sizes');
end

if(length(x)~=nx)
    error('x and seis have incompaible sizes');
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
for k=kx
    n=n+1;
    kk=max([k-delx, 1]):min([k+delx,nx]);
    jx=near(t,t1x(k));
    s=seis(jx-iwin2:jx+iwin2,kk);
    sr=sref(:,ones(size(kk)));
    xx=maxcorr_ephs(s(:),sr(:),nlag);
    phsx(n)=-xx(3);
    ccx(n)=xx(4);
    delayx(n)=-xx(2)*dt;
end

% phs=interpextrap(x(kx),phsx,x,0);
% delay=interpextrap(x(kx),delayx,x,0);
% cc=interpextrap(x(kx),ccx,x,0);
phs=interp1(x(kx),phsx,x,'spline','extrap');
delay=interp1(x(kx),delayx,x,'spline','extrap');
cc=interp1(x(kx),ccx,x,'spline','extrap');

end

function sref=referencegather(seis,t,x,xrefs,n)
%xrefs ... coordinates of the desired reference traces
%the first xrefs will be the location to which all the others are aligned
% n= 2*n +1 lags will be computed
% ******* default= round(length(t)/10) *********

if(nargin<5)
    n=round(length(t)/10);
end
dt=t(2)-t(1);
sref=zeros(length(t),length(xrefs));
for k=1:length(xrefs)
    ix=near(x,xrefs(k));
    if(k==1)
        sref(:,k)=seis(:,ix(1));
    else
        s2=seis(:,ix(1));
        s1=sref(:,1);
        x=maxcorr_ephs(s1,s2,n);
        sref(:,k)=phsrot(stat(s2,t,dt*x(2)),x(3));
    end
    
end


end



