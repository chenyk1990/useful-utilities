function [seisf,masks,kx,ky]=fkfanfilter3D(seis,t,x,y,va1,va2,dv,tpad,xpad,ypad,fmask)
%FKFANFILTER3D ... apply an f-k fan filter to a 3D seismic dataset
%
% [seisf,masks,kx,ky]=fkfanfilter3D(seis,t,x,y,va1,va2,dv,tpad,xpad,ypad,fmask)
%
% FKFANFILTER3D designs and applies an f-k (frequency-wavenumber) fan reject filter in 3D. The
% reject region is fan-shaped (when viewed in frequency and radial wavenumber) and defined by two
% bounding apparent velocities, va1 and va2. These are expressed as positive numbers and the filter
% is rotationally invariant about the vertical (t) axis. A raised cosine taper is applied to the
% filter edge. The filter is applied to each frequency as a radially symmetric multiplier (or mask)
% whose values lie between 0 and 1. This mask can be examined by returning all four values and then
% plotting it.
%
% seis ... 3D seismic matrix to be filtered. The first dimension is time, the second x and the third
%       y. All three coordinates must be regularly sampled.
% t ... time coordinate vector for seis. The length(t) must equal size(seis,1)
% x ... first space coordinate vector for seis. The length(x) must equal size(seis,2)
% y ... second space coordinate vector for seis. The length(y) must equal size(seis,3)
% va1 ... minumum apparent velocity defining the rejection fan. 
%       Enter 0 to reject everything slower than va2.
% va2 ... maximum apparent velocity defining the rejection fan.
% Requirement: va2>va1. neither value can be negative.
% dv  ... width of the taper on the edge of the rejection fan in velocity units
% REQUIREMENT: 0<dv<va1<=va2.  Setting va1=va2 gives a very narrow reject
% region. Better rejection of a specific velocity, vn, follows from making
% va1 slightly less than vn and va2 slightly greater.
%
% tpad ... size (in t units) of temporal zero pad to be afixed to seis in the first dimension
% ********* default = 0.1*(max(t)-min(t)) **********
% xpad ... size (in x units) of spatial zero pad to be afixed to seis in the second dimension.
% ********* default = 0.1*(max(x)-min(x))***********
% ypad ... size (in y units) of spatial zero pad to be afixed to seis in the third dimension.
% ********* default = 0.1*(max(x)-min(x))***********
% fmask ... vector of frequencies at which the filter mask is to be output
% ********* default = [] (no mask output) **********
%
% NOTE: The values supplied for xpad, ypad and tpad are minimal pads because, after afixing these pads,
% the matrix is further extended to the next power of 2 in all dimensions. The purpose of this pad
% is to minimize wraparound of the filter impulse response.
%
% seisf ... the f-k filtered result with all pads removed. It will be the same size as seis.
% masks ...a cell array of length fmask of filter multipliers. Each cell contains a 2-D mask with 
%       kx as the row coordinate and ky as the column coordinate.
% kx ... x wavenumber coordinate for mask
% ky ... y wavenumber coordinate for mask
% 
% The mask can be displayed by: plotimage(mask{ifreq},kx,ky) where ifreq is an integer corresponding
% to one of the input values of fmask.
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
disp('fkfanfilter3D')
if(nargin<11)
    fmask=[];
end
if(nargin<8)
    tpad=0.1*(max(t)-min(t));
end
if(nargin<10)
    ypad=0.1*(max(y)-min(y));
end
if(nargin<9)
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
if(dv<=0 || va1<0 || va2<=0)
    error('dv, va1,va2 must not be negative')
end
if(va1>va2)
    error('va1 must be less than va2')
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
kx=kx(:);
ky=ky(:)';
nkx=length(kx);
nky=length(ky);
kkx=kx(:,ones(1,nky));
kky=ky(ones(nkx,1),:);
kr=sqrt(kkx.^2+kky.^2);%radial wavenumber

%express the filter reject regions in slowness
sa1=1/va1;%might be infinite
sa2=1/va2;
if(~isinf(sa1))
    sa1p=sa1*(1+dv/va1);
else
    sa1p=inf;
end
sa2p=sa2*(1-dv/va2);
% If sr is radial slowness on a constant f plane, then filter is all-pass for 
% s2<sa2p, is all reject for sa2<sr<sa1, and tapers from all-pass to all
% reject for sa2p<sr<sa2 and tapers from all-reject to all-pass for
% sa1<sr<sa1p. 
%


ifmask=round(fmask/df)+1;
nfmask=1;
masks=cell(size(ifmask));

for k=2:length(f)%skip 0 hz
   mask=ones(nkx,nky);
   sr=kr/f(k);%radial slowness at all points on the grid
   ind1=find(sr>sa2p);
   ind2=find(sr(ind1)<=sa2);
   mask(ind1(ind2))=.5+.5*cos((sr(ind1(ind2))-sa2p)*pi/(sa2-sa2p));%first taper
   if(~isinf(sa1))
       ind1=find(sr>sa2);
       ind2= sr(ind1)<=sa1;
       mask(ind1(ind2))=0;%all reject
       ind1=find(sr>sa1);
       ind2=find(sr(ind1)<=sa1p);
       mask(ind1(ind2))=.5+.5*cos((sr(ind1(ind2))-sa1p)*pi/(sa1-sa1p));%second taper
   else
       %this is an evanescent filter with va2 defining the edge
       ind1= sr>sa2;
       mask(ind1)=0;
   end
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
 