function [seisf,masks,kx,ky]=fkevfilter3Dani(seis,t,x,y,vax,vay,dvx,dvy,tpad,xpad,ypad,fmask)
%FKEVFILTER3Dani ... apply an f-k evanescent filter to a 3D seismic dataset
%
% [seisf,masks,kx,ky]=fkevfilter3Dani(seis,t,x,y,vax,vay,dvx,dvy,tpad,xpad,ypad,fmask)
%
% FKEVFILTER3Dani designs and applies an f-k (frequency-wavenumber) evanescent reject filter in 3D.
% All apparent velocities slower thatn the specified filter boundary are rejected. In cross section, the reject
% region is elliptical where the ellipse is oriented with axes parallel to the x and y axes.
% This requires specifying the boundary apparent
% velocities for both the x and y directions. These are expressed as positive numbers and the
% rejection region in any radial direction is determined by the elliptical "radius" in that
% direction. A raised cosine taper is applied to the filter edge. The filter is applied to each
% frequency as a real-valued multiplier (or mask) whose values lie between 0 and 1. This mask can be
% examined by returning all four values and then plotting it.
%
% seis ... 3D seismic matrix to be filtered. The first dimension is time,
%       the second x and the third y. All three coordinates must be
%       regularly sampled.
% t ... time coordinate vector for seis. The length(t) must equal
%       size(seis,1)
% x ... first space coordinate vector for seis. The length(x) must equal
%       size(seis,2)
% y ... second space coordinate vector for seis. The length(y) must equal
%       size(seis,3)
% vax ... apparent velocity defining the rejection fan in the x direction. 
% vay ... apparent velocity defining the rejection fan in the x direction. 
%   Requirement: Both vax and vay must be positive
% dvx  ... with of the taper in the x direction for the edge of the rejection fan in velocity units
% dvy  ... with of the taper in the y direction for the edge of the rejection fan in velocity units
%   Note: Very small dvx and dvy will likely cause visible filter artefacts. A reasonable choice might
%   be 10% of the vax and vay values.
%
% tpad ... size (in t units) of temporal zero pad to be afixed to seis in
%       the first dimension
% ********* default = 0.1*(max(t)-min(t)) **********
% xpad ... size (in x units) of spatial zero pad to be afixed to seis in
%       the second dimension.
% ********* default = 0.1*(max(x)-min(x))***********
% ypad ... size (in y units) of spatial zero pad to be afixed to seis in
%       the third dimension.
% ********* default = 0.1*(max(x)-min(x))***********
% fmask ... vector of frequencies at which the filter mask is to be output
% ********* default = [] (no mask output) **********
%
% NOTE: The values supplied for xpad and tpad are minimal pads because,
% after afixing these pads, the maxtrix is further extended to the next
% power of 2 in all dimensions. THe purpose of this pad is to minimize
% wraparound of the filter impulse response.
%
% seisf ... the f-k filtered result with all pads removed. It will be the
%       same size as seis.
% masks ...a cell array of length fmask of filter multipliers. Each cell
%       contains a 2-D mask with kx as the row coordinate and ky as the
%       column coordinate.
% kx ... x wavenumber coordinate for mask.
% ky ... y wavenumber coordinate for mask.
% 
% The mask can be displayed by: plotimage(mask{ifreq},kx,ky) where ifreq is
% an integer corresponding to one of the input values of fmask.
% 
% G.F. Margrave, Devon Energy, 2017
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE
disp('fkfanfilter3Dani')
if(nargin<12)
    fmask=[];
end
if(nargin<9)
    tpad=0.1*(max(t)-min(t));
end
if(nargin<11)
    ypad=0.1*(max(y)-min(y));
end
if(nargin<10)
    xpad=0.1*(max(x)-min(x));
end

[nt,nx,ny]=size(seis);
if(length(t)~=nt)
    error('seis and t have incompatible sizes')
end
if(length(x)~=nx)
    error('seis and x have incompatible sizes')
end
if(length(y)~=ny)
    error('seis and y have incompatible sizes')
end

if(dvy<=0 || dvx<=0 || vax<0 || vay<=0)
    error('dvx, dvy, vax,vay must not be negative')
end

dx=abs(x(2)-x(1));
dy=abs(y(2)-y(1));
dt=t(2)-t(1);
small=.00001*dx;
if(sum(abs(diff(diff(x))))>small)
    error('x coordinates must be regularly spaced')
end
small=.00001*dy;
if(sum(abs(diff(diff(y))))>small)
    error('y coordinates must be regularly spaced')
end
small=.00001*dt;
if(sum(abs(diff(diff(t))))>small)
    error('t coordinates must be regularly spaced')
end

%attach t pad
nt2=nt;
if(tpad>0)
    ntpad=round(tpad/dt);
    nt2=nt+ntpad;
    t=(0:nt2-1)*dt;
    seis=[seis;zeros(ntpad,nx,ny)];
end
%attach x and y pads
if(xpad>0 || ypad>0)
    nxpad=round(xpad/dx);
    nypad=round(ypad/dy);
    seis2=zeros(nt2,nx+nxpad,ny+nypad);
    seis2(:,1:nx,1:ny)=seis;
    nx2=nx+nxpad;
    x=(0:nx2-1)*dx;
    ny2=ny+nypad;
    y=(0:ny2-1)*dy;
else
    seis2=seis;
end
clear seis;

%fk transform
t0=clock;
disp('Beginning forward transform');
[seisfk,f,kx,ky]=fkktran(seis2,t,x,y);
tnow=clock;
timeused=etime(tnow,t0);
disp(['forward transform complete in ' int2str(timeused) ' seconds'])
df=f(2);
nf=length(f);
kx=kx(:);
ky=ky(:)';
dkx=kx(2)-kx(1);
dky=ky(2)-ky(1);
nkx=length(kx);
nky=length(ky);
kkx=kx(:,ones(1,nky));
kky=ky(ones(nkx,1),:);
%kr=sqrt(kkx.^2+kky.^2);%radial wavenumber

%express the filter reject regions in slowness
sax=1/vax;%x slowness of filter edge
say=1/vay;%y slowness of filter edge
dsx=1/(vax-dvx)-sax;
dsy=1/(vay-dvy)-say;
% If sr is radial slowness on a constant f plane, then filter is all-pass for 
% s2<sa2p, is all reject for sa2<sr<sa1, and tapers from all-pass to all
% reject for sa2p<sr<sa2 and tapers from all-reject to all-pass for
% sa1<sr<sa1p. 
%


ifmask=round(fmask/df)+1;
nfmask=1;
masks=cell(size(ifmask));
theta=0:360;
ievery=50;
for k=2:nf%skip 0 hz
% ind=near(f,30);
% for k=ind
   mask=zeros(nkx,nky);
   %sr=kr/f(k);%radial slowness at all points on the grid
   %define an ellipse with the x and y slownesses
   a=f(k)*sax;
   b=f(k)*say;
   r=a*b./sqrt((b*cosd(theta)).^2+(a*sind(theta)).^2);
   kxe=r.*cosd(theta);
   kye=r.*sind(theta);
   %kxe and kye are the coordinates of a 360 point polygonal approximation to the ellipse
   ind=inside(kkx,kky,kxe,kye);
   %ind=inpolygon(kkx,kky,kxe,kye);
   mask(ind)=1;
   %now convolve with a radial gaussian to get a taper
   sigmax=max([3*dkx f(k)*dsx]);
   sigmay=max([3*dky f(k)*dsy]);
   g=gaus_radial(dkx,sigmax,dky,sigmay);
   mask=conv2(mask,g,'same')/sum(g(:));
   
   %apply the mask
   tmp=squeeze(seisfk(k,:,:));
   seisfk(k,:,:)=mask.*tmp;
   if(~isempty(fmask))
       if(nfmask<=length(fmask))
           if(k==ifmask(nfmask))
               %the transpose is here because we want kx to be the column coordinate. However, the
               %squeeze command eliminates the first dimension (f) and promotes kx to dimension 1
               %and ky to dimension 2. This makes kx the row coordinate so we transpose.
               masks{nfmask}=mask';
               nfmask=nfmask+1;
           end
       end
   end
   if(rem(k,ievery)==0)
       tnow1=clock;
       timeused=etime(tnow1,tnow);
       timeperf=timeused/(k-1);
       timeleft=timeperf*(nf-k);
       disp(['completed frequency number ' int2str(k) ' of ' int2str(nf)]);
       disp(['time used = ' int2str(timeused) ' s, time remaining= ' int2str(timeleft) ' s'])
   end
end
tnow2=clock;
timeused=etime(tnow2,tnow);
disp(['Filter applied in ' int2str(timeused) ' seconds'])
%inverse transform
disp('Beginning inverse transform');
seisf=ifkktran(seisfk,f,kx,ky,nt,nx,ny);
tnow3=clock;
timeused=etime(tnow3,tnow2);
disp(['Inverse transform completed in ' int2str(timeused) ' seconds'])
timeused=etime(tnow3,t0);
disp(['3D fan filter total time ' int2str(timeused) ' seconds'])
 