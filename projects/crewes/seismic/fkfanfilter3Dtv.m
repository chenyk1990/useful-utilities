function seisf=fkfanfilter3Dtv(seis,t,x,y,twins,toverlap,va1,va2,dv,tpad,xpad,ypad)
%FKFANFILTER3Dtv ... time-variant 3D f-k fan filter
%
% seisf=fkfanfilter3Dtv(seis,t,x,y,va1,va2,dv,tpad,xpad,ypad)
%
% This function breaks a 3D seismic dataset into overlapping temporal windows and applies
% fkfanfilter3D to each window with window-dependent parameters. Any number of windows can be
% specified and the window size and overlap are user prescribed. FKFANFILTER3D designs and applies
% an f-k (frequency-wavenumber) fan reject filter in 3D. The reject region is fan-shaped (when
% viewed in frequency and radial wavenumber) and defined by two bounding apparent velocities, va1
% and va2. These are expressed as positive numbers and the filter is rotationally invariant about
% the vertical (t) axis. A raised cosine taper is applied to the filter edge. The filter is applied
% to each frequency as a radially symmetric multiplier (or mask) whose values lie between 0 and 1.
% This mask can be examined by returning all four values and then plotting it.
%
%
% seis ... 3D seismic matrix to be filtered. The first dimension is time, the second x and the third
%       y. All three coordinates must be regularly sampled.
% t ... time coordinate vector for seis. The length(t) must equal size(seis,1)
% x ... first space coordinate vector for seis. The length(x) must equal size(seis,2)
% y ... second space coordinate vector for seis. The length(y) must equal size(seis,3)
% twins ... vector of window boundary times. Do not include t(1) and t(end) as they are automatic.
%       The first window extends from t(1) to twins(1), the second from twins(1) to twins(2), and so
%       on until the last which is twins(end) to t(end). May be specified as [] for a single window
%       but in that case you should just use fkfanfilter3D directly.
% toverlap ... vector of window overlap times for each window. May be a single value if all windows
%       are to have the same overlap. This means that window 1 actually extends from t(1) to
%       twins(1)+toverlap(1) where the last toverlap(1) samples are multiplied by a raised cosine
%       taper. Similarly window 2 extends from twins(1)-toverlap(1) to twins(2)+toverlap(2) (with
%       tapering) and so on.
% va1 ... vector of minumum apparent velocities defining the rejection fan in each window. Enter 0 to
%       reject everything slower than va2. The length of va1 must be 1 greater than the
%       length(twins).
% va2 ... vector of maximum apparent velocities defining the rejection fan in each window. The
%       length(va2) must equal length(twins)+1.
% Requirement: va2>va1. Neither value can be negative.
%
% dv  ... vector of widths of the taper on the edge of the rejection fan in velocity
%         units (for each window).
%
% REQUIREMENT: 0<dv<va1<=va2.  Setting va1=va2 gives a very narrow reject region. Better rejection
% of a specific velocity, vn, follows from making va1 slightly less than vn and va2 slightly
% greater.
%
% tpad ... size (in t units) of temporal zero pad to be afixed to each windowed seismic section in
%           the first dimension.
% ********* default = 0.1*size_of_window **********
% xpad ... size (in x units) of spatial zero pad to be afixed to seis in the second dimension.
% ********* default = 0.1*(max(x)-min(x))***********
% ypad ... size (in y units) of spatial zero pad to be afixed to seis in the third dimension.
% ********* default = 0.1*(max(x)-min(x))***********
%
% NOTE: The values supplied for xpad, ypad, and tpad are minimal pads because, after afixing these pads,
% the matrix is further extended to the next power of 2 in all dimensions. The purpose of this pad
% is to minimize wraparound of the filter impulse response.
%
% seisf ... the f-k filtered result with all pads removed. It will be the same size as seis.
% 
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

if(nargin<10)
    tpad=0.1*(max(t)-min(t));
end
if(nargin<12)
    ypad=0.1*(max(y)-min(y));
end
if(nargin<11)
    xpad=0.1*(max(x)-min(x));
end

twins=twins(:);
va1=va1(:);
va2=va2(:);
dv=dv(:);
small=1000*eps;
ind=find(abs(twins-t(1))<small, 1);
if(isempty(ind))
    twins=[t(1);twins];
end
ind=find(abs(twins-t(end))<small, 1);
if(isempty(ind))
    twins=[twins;t(end)];
end
nwins=length(twins)-1;
if(length(va1)~=nwins)
    error('va1 is the wrong size');
end
if(length(va2)~=nwins)
    error('va2 is the wrong size');
end
if(length(dv)~=nwins)
    error('dv is the wrong size');
end
if(length(toverlap)==1)
    toverlap=toverlap*ones(size(va1));
elseif(length(toverlap)~=nwins)
    error('toverlap is the wrong size');
end
toverlap=toverlap(:);
%force twins to fall on a sample
dt=t(2)-t(1);
twins=dt*round(twins/dt);
nhoverlap=floor(.5*toverlap/dt);
nhoverlap=[0;nhoverlap];
nhoverlap(end)=0;
seisf=zeros(size(seis),'like',seis);
disp('fkfanfilter3Dtv');
tnot=clock;
for k=1:nwins
    tw1=twins(k);%beginning of window
    tw2=twins(k+1);%end of window
    v1=va1(k);
    v2=va2(k);
    dvk=dv(k);
    iwin1=round((tw1-t(1))/dt)+1-nhoverlap(k);%first sample in window
    iwin2=round((tw2-t(1))/dt)+1+nhoverlap(k+1);%last sample in window
    iwin=(iwin1:iwin2)';
    mw=ones(size(iwin));
    %front taper
    i1=1;
    i2=i1+2*nhoverlap(k);
    ii=i1:i2;
    if(i1~=i2)
        mw(ii)=.5+.5*cos((ii-i2)*pi/(i1-i2));
    end
    %end taper
    nw=length(mw);
    i2=nw;
    i1=nw-2*nhoverlap(k+1);
    ii=i1:i2;
    if(i1~=i2)
        mw(ii)=.5+.5*cos((ii-i1)*pi/(i2-i1));
    else
        mw(ii)=0;
    end
    mw3d=mw(:,ones(size(x)),ones(size(y)));
    seisw=seis(iwin,:,:).*mw3d;%windowed seismic
    seiswf=fkfanfilter3D(seisw,t(iwin),x,y,v1,v2,dvk,tpad,xpad,ypad);
    seisf(iwin,:,:)=seisf(iwin,:,:)+seiswf;
    tnow=clock;
    timeused=etime(tnow,tnot);
    %assume run time is proportional to the number of time samples
    timeremaining=(length(t)-iwin(end))*timeused/iwin(end);
    disp(['finished layer ' int2str(k) ' after ' int2str(timeused) ' s'])
    disp(['estimated time remaining ' int2str(timeremaining) ' s'])
end

    
    