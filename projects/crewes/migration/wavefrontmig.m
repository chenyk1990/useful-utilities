function seismig=wavefrontmig(seis,t,x,v,aper,filtparms,increments)
% WAVEFRONTMIG: perform a constant velocity migration by direct superposition of wavefronts
%
% seismig=wavefrontmig(seis,t,x,v,aper,filtparms,increments)
%
% This is a post-stack (exploding reflector) migration. Its mostly for teaching/learning.
%
% seis ... input seismic matrix
% t ... time coordinate for seis. Length(t) must equal size(seis,1)
% x ... x coordinate for seis. Length(x) must equal size(seis,2)
% v ... velocity (must be a scalar). Will be divided by 2 internally for exploding reflector
% aper ... maximum allowed radius for a wavefront circle. 
% filtparms ... length 4 vector giving f1,f2,f3,f4 Ormsby filter parameters. Set the
%           values to nan for no filter.
% *********** default = nan*ones(1,4) which means no filter ***********
% increments ... length 2 vector giveing the time and space increments (in samples) for
%           wavefronts. For example increments=[10 5] means that every 10th time sample (row)
%           on every 5th spatial coordinate (column) will have a wavefront.
% *********** default =[1 1] ***********
%
% G.F. Margrave, CREWES, 2014
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
%
%
if( nargin <5 )
    aper = inf;
end
if(nargin<6)
    filtparms=nan*1:4; %filter specification Ormsby style
else
    filtparms=filtparms(1:4);
end
if(nargin<7)
    increments=[1 1];
end

v=v/2;
dt=t(2)-t(1);
tmin=t(1);

[nsamps,nc]=size(seis);
seismig=zeros(nsamps,nc);

tstart=clock;

for k=1:increments(2):nc
    xnot=x(k);
    for k2=1:increments(1):nsamps
        tnot=t(k2);
        r=v*tnot;%radius of the circle
        r=min([r aper]);
        x1=xnot-r;
        x2=xnot+r;
        ix=near(x,x1,x2);%these are the columns of interest
        %loop over wavefront columns
        for k3=ix
            xoff=x(k3)-xnot;
            tk = sqrt(r^2-xoff^2)/v;
            a=seis(k2,k);
            ik=(tk-tmin)/dt+1;
            if( between(1,nsamps,ik) )
                ik1=floor(ik);
                ik2=ceil(ik);
                if(ik1==ik2)
                    %exactly on a sample
                    seismig(ik1,k3)=seismig(ik1,k3)+a;
                else
                    %a simple interpolation
                    seismig(ik1,k3)=seismig(ik1,k3)+a*(ik-ik2)/(ik1-ik2);
                    seismig(ik2,k3)=seismig(ik2,k3)+a*(ik-ik1)/(ik2-ik1);
                end
            end
        end
    end
    if(rem(k,100)==0)
        tnow=clock;
        timeused=etime(tnow,tstart);
        time_per_col=timeused/k;
%         timetotal=nc*time_per_col;
        timeremaining=(nc-k)*time_per_col;
        disp(['finished column ' int2str(k) ' of ' int2str(nc) ', estimated time remaining ' int2str(timeremaining) 's'])
    end
end

%filter
if(~isnan(filtparms(1)))
    seismig=filtf(seismig,t,[filtparms(2) filtparms(2)-filtparms(1)],[filtparms(3) filtparms(4)-filtparms(3)]);
end