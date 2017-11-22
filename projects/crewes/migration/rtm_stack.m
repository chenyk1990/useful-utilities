function [seism,xm,zm,snapshots,xms,zms]=rtm_stack(seis,t,x,dtstep,vel,tout,laplacian,boundary)
% RTM_STACK ... 2D post-stack reverse-time migration (RTM) (a depth migration)
%
% [seism,xm,zm,snapshots,xms,zms]=rtm_stack(seis,t,x,dtstep,vel,tout,laplacian,boundary)
%
% RTM_STACK does post-stack reverse-time depth migration. An input zos
% (zero-offset section) is injected into the subsurface and propagated
% backwards in time using second-order time stepping. Time-stepping is
% initiallized with snapshots at Tmax (=max(t)) and Tmax+dtstep.
% The latter snapshot is all zeros while the former is all zeros except for
% the z=0 samples which are injected from the stack at Tmax. These are
% sufficient to calculate the snapshot at Tmax-dtstep. This process is then
% iterated until the snapshot at t=0 is calculated which is the migrated
% section. At each step, a new row of samples are injected from the
% stacked section. The spatial derivatives can be calculated with a five or
% nine point approximation to the Laplacian operator.  The five point
% approximation is faster, but the nine point results in a broader
% bandwidth.  The two laplacian options have different stability conditions
% (see below).
%
% seis ... the input seismic stack to be migrated. Each column is a trace
%       and each row is a time sample.
% t ... time coordinate vector for seis. length(t) must equal size(seis,1)
% x ... space coordinate vecotr for seis. length(x) must equal size(seis,2)
% NOTE: Both t and x must be regularly sampled
% dtstep = size of time step (in seconds). dtstep must be no greater than
%       dt=t(2)-t(1). For best results it should be about dt/10 .
% vel = the input velocity matrix. This must be instantaneous velocity as a
%       function of x (columns) and z (rows). 
% NOTE: Vel must be on a square regular grid of spatial sample size dx=x(2)-x(1);
% NOTE: do not divide by two to compensate for one way travel time. This is built into the program.
% tout = array of times at which output snapshots are desired
% ********** default = [] (no snapshots) ************
% laplacian - an option between two approximation to the laplacian operator
%           - 1 is a 5 point approximation
%          Stability condition: max(velocity)*dtstep/dx MUST BE < sqrt(2)
%           - 2 is a nine point approximation
%          Stability condition: max(velocity)*dtstep/dx MUST BE < 2*sqrt(3/8)
% ********* default = 2 ************
% boundary = indicate whether all sides of the matrix are absorbing
%          = 0 indicates that no absorbing boundaries are desired
%          = 1 indicates all four sides are absorbing
%          = 2 choses three sides to be absorbing, and the top one not to be
%             this enables sources to be put on the surface
% ********** default = 2 *********
%
% seism = the output migrated depth section
% xm = the x coordinate of the migrated section (same as x)
% zm = the depth coordinate of the migrated section
% snapshots = cell array of snapshots of the reverse-time propagated wavefield at the times tout
% xms = x coordinate for the snapshots
% zms = z coordinates for the snapshots
%
% NOTE: see test_rtm_stack for an example of running this and viewing the snapshots
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

small=10000*eps;

disp('Beginning post-stack reverse time migration');

if(nargin<8)
    boundary=2;% means all top is not absorbing
end
if(nargin<7)
    laplacian=2;
end
if(nargin<6)
    tout=[];
end

if(length(x)~=size(seis,2))
    error('length(x) must equal size(seis,2)')
end
if(length(x)~=size(vel,2))
    error('length(x) must equal size(vel,2)')
end
if(length(t)~=size(seis,1))
    error('length(t) must equal size(seis,1)')
end
if(sum(abs(diff(diff(t))))>small)
    error('t must be regularly sampled')
end
if(sum(abs(diff(diff(x))))>small)
    error('x must be regularly sampled');
end

if(laplacian==1)
    disp('using second order Laplacian')
elseif(laplacian==2)
    disp('using fourth order Laplacian')
else
    error('Invalid value for laplacian')
end


vel=vel/2; %exploding reflector velocity

unstable=0;
dxin=abs(x(2)-x(1));
dt=t(2)-t(1);

if(dt<dtstep)
	error('dt cannot be less than dtstep')
end

vmax=max(vel(:));
if laplacian ==1 
    if vmax*dtstep/dxin > 1/sqrt(2) %note vel has been halved
       unstable=1;
    end
else
    if vmax*dtstep/dxin > sqrt(3/8)
       unstable=1;
    end
end

z=dxin*(0:size(vel,2)-1)';
tstep=(t(1):dtstep:t(end))';
if(unstable)
    error('input parameters are unstable');
    %need to intepolate vel to a coarser spatial grid
%     if(laplacian==1)
%         dxnew=vmax*dtstep*sqrt(2);
%     else
%         dxnew=vmax*dtstep*sqrt(8/3);
%     end
%     %want (x(end)-x(1))/dxnew to be an integer
%     nx=ceil(x(end)-x(1))/dxnew;
%     dxnew=(x(end)-x(1))/nx;
%     disp(['Resampling spatial grid from ' num2str(dxin) ' to ' num2str(dxnew)])
%     xnew=x(1):dxnew:x(end);
%     znew=(z(1):dxnew:z(end))';
%     velnew=interp2(x(ones(size(z)),:),z(:,ones(size(x))),vel,xnew,znew);
%     vel=velnew;
%     seisnew=zeros(length(t),length(xnew));
%     for k=1:length(t)
%         seisnew(k,:)=interp1(x,seis(k,:),xnew,'spline');
%     end
%     seis=seisnew;
%     x=xnew;
%     clear velnew seisnew xnew znew
end

if(dtstep~=dt)
    %interpolate seismic to dtstep
    seisnew=zeros(length(tstep),length(x));
    for k=1:length(x)
        seisnew(:,k)=interpbl(t,seis(:,k),tstep);
    end
    seis=seisnew;
    clear seisnew
end

% xmax=max(x);zmax=max(z);
snap1=zeros(size(vel));
snap2=snap1;
snap2(1,:)=seis(end,:);%inject last row of seis at top of snapshot

nsteps=length(tstep);
disp(['There are ' int2str(nsteps) ' steps to complete']);
time0=clock;

nwrite=2*round(nsteps/50)+1;
dx=abs(x(2)-x(1));

if(~isempty(tout))
    ind=find(tout<=tstep(end));
    snapshots=cell(size(tout(ind)));
    iout=nsteps-round(tout(ind)/dtstep);
    kout=0;
    iout2=zeros(size(iout));
end

for k=1:nsteps
    
    snap3=afd_snap(dx,dtstep,vel,snap1,snap2,laplacian,boundary);
    %snap1 is at time t+dtstep, snap2 is at time t, and snap3 is at time t-dtstep 
    snap1=snap2;%snap2 moves to t+dtstep
    snap2=snap3;%snap3 moves to t
    
    if(k<nsteps)
        snap2(1,:)=snap2(1,:)+seis(end-k,:);%inject the next seismic sample
    end
    
    if rem(k,nwrite) == 0
        timenow=clock;
        timeused=etime(timenow,time0);
        timeremaining=(nsteps-k)*timeused/k;
        
        disp(['Completed ' num2str(k) ' steps in ' int2str(timeused) ' seconds'])
        disp(['Time remaining ' int2str(timeremaining) ' seconds or ' num2str(timeremaining/60) ' minutes'])
    end
    
    if(k==nsteps)
        disp('hey')
    end
    
    ind=find(k==iout, 1);
    if(~isempty(ind))
        kout=kout+1;
        snapshots{kout}=snap3;
        iout2(kout)=iout(ind);
    end
    
end

if(~isempty(tout))
    xms=x(1)+(0:size(snap3,2)-1)*dx;
    zms=(0:size(snap3,1)-1)*dx;
else
    snapshots=[];
    xms=[];
    zms=[];
end

disp('RTM completed')

%downsample to output x and z
if(dxin>dx)
    ievery=round(dxin/dx);%should be an integer
    seism=snap3(1:ievery:end,1:ievery:end);
else
    seism=snap3;
end
xm=x(1)+(0:size(seism,2)-1)*dxin;
zm=(0:size(seism,1)-1)*dxin;


