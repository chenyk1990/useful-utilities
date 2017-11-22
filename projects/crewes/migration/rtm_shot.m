function [shotm,shotmf,xm,zm,recsnaps,sousnaps,refsnaps,xms,zms]=rtm_shot(shot,t,xrec,zrec,xshot,zshot,w,tw,vel,xv,dtstep,itcorr,tout,laplacian,illume,boundary)
% RTM_SHOT ... 2D pre-stack reverse-time migration (RTM) of a shot record (a depth migration)
%
% [shotm,shotmf,xm,zm,recsnaps,sousnaps,refsnaps,xms,zms]=rtm_shot(shot,t,xrec,zrec,xshot,zshot,w,tw,vel,xv,dtstep,itcorr,tout,laplacian,illume,boundary)
%
% RTM_SHOT does pre-stack reverse-time depth migration of a shot record. The migration computes the
% forward-in-time propagation of the source field and the reverse-time propagation of the input shot
% (receiver field). Thus two wavefields (source wavefield and receiver wavefield) are simulated
% by finite-difference time stepping. The migrated shot (or reflectivity estimate) comes from
% crosscorrelation of these two fields at equal times. The crosscorrelation is a 2D (x,z) dot-star
% (.*) operation. This requires that the source wavefield must be available at each time step of the
% receiver field and this is a major complexity. Here this is addressed by fully propagating the
% source field before the reverse-time propagation is begun and saving the needed snapshots of the
% source field in a 3D array. Rather than correlation at each time step (dtstep) the correlation is
% done at time intervals of itcorr*dtstep where itcorr>=1. This allows a control on the memory
% requirements and itcorr must be an integer. A good choice for itcorr is round(dt/dtstep). Also
% required is a wavelet for the source simulation. Ideally this should be the same wavelet that is
% embedded in the receiver field (shot record). All time stepping is second order. At each step of
% the reverse-time propagation, a new row of samples from the shot record is injected into the
% receiver field at the receiver locations. The spatial derivatives can be calculated with a five or
% nine point approximation to the Laplacian operator. The five point approximation is faster, but
% the nine point results in a broader bandwidth.  The two laplacian options have different stability
% conditions (see below).
%
% shot ... the input seismic stack to be migrated. Each column is a trace
%       and each row is a time sample.
% t ... time coordinate vector for shot. length(t) must equal size(shot,1)
% xrec ... space coordinate vector for shot. length(x) must equal size(shot,2)
% zrec ... depth coordinates of the receivers. Enter a scalar 0 to place all rreceivers at the z=0,
%       otherwise, it must be a vector the same size as xrec.
% NOTE: t must be regularly sampled but xrec need not be. The propagation grid size is determined by xv.
% xshot, zshot ... x and z coordinates of the shot (scalars for point souce, same-length vectors for an array)
% w,tw ... the wavelet and its time coordinate. The wavelet may be noncausal. Must be sampled at dtstep.
% vel = the input velocity model. This must be instantaneous velocity as a function of x (columns) 
%       and z (rows). Vel must be on a square regular grid of spatial sample size dx=xv(2)-xv1);
% xv = horizontal coordinate for the velocity model. The velocity model can span a longer x range
%       than the shot record but must be no smaller. The depth coordinate of vel is generated
%       internally as dx*(0:size(vel,1)-1)' . The dx of vel determines the square propagation grid
%       on thich the computation occurs.
% dtstep = size of time step (in seconds). dtstep must be no greater than dt=t(2)-t(1). 
%       For best results dtstep should be about dt/10 .
% itcorr = the source and receiver fields will be correlated (.*) every itcorr*dtstep seconds. itcorr
%       must be an integer no smaller than 1. Note that setting this equal to 1 will require large
%       memory. The .* operation gets the instantaneous reflectivity estimate while the final
%       reflectivity is the sum over time of the instantaneous reflectivites.
% *********** default = round(dt/dtstep) where dt=t(2)-t(1) *********
% tout = array of times at which output snapshots are desired
% ********** default = [] (no snapshots) ************
% laplacian - an option between two approximations to the laplacian operator
%           - 1 is a 5 point approximation
%          Stability condition: max(velocity)*dtstep/dx MUST BE < sqrt(2)
%           - 2 is a nine point approximation
%          Stability condition: max(velocity)*dtstep/dx MUST BE < sqrt(3/8)
% ********* default = 2 ************
% NOTE: dx in the stability formulae is dx=abs(xv(2)-xv(1));
% illume ... if 1 then the illumination is computed at each correlation step and used to normalize
%           the correlation. Illumination is the square (.^2) of the source field. Otherwise, the
%           reflectivity estimate at each time is simply the correlation (.*) of the source and
%           receiver fields. Division of the correlation by the illumination (./) is effectivlty a
%           correction for the amplitude loss caused by wavefront spreading.
% ********* default = 1 ************
% boundary = indicate which sides of the velocity matrix are absorbing
%          = 0 indicates that no absorbing boundaries are desired
%          = 1 indicates all four sides are absorbing
%          = 2 choses sides and bottom to be absorbing, and the top not to be. (Top is a free
%          surface). This enables sources to be put on the surface because placing them on an
%          absorbing boundary tends to be defeating.
% ********** default = 2 *********
%
% shotm = the migrated shot in depth. This will be the same size as the velocity model.
% shotmf = the migrated shot with Laplacian filter applied
% xm = the x coordinate of the migrated shot (same as xv)
% zm = the depth coordinate of the migrated shot (same as depth axis of the velocity model)
% recsnaps = cell array of snapshots of the reverse-time propagated receiver wavefield at the times tout
% sousnaps = cell array of snapshots of the forward-time propagated source wavefield at the times tout
% refsnaps = cell array of snapshots of the reflectivity image (the migrated shot) at the times tout
% xms = x coordinate for the snapshots
% zms = z coordinates for the snapshots
%
% NOTE: The depth migrated shot will be the same size as the velocity model. Its up to the user
% to ensure this is large enough to prevent unwanted boundary interactions.
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

disp('Beginning pre-stack reverse time migration');

if(nargin<16)
    boundary=2;% means all top is not absorbing
end
if(nargin<15)
    illume=1;% means we will divide by the illumination
end
if(nargin<14)
    laplacian=2;
end
if(nargin<13)
    tout=[];
end
dt=t(2)-t(1);
if(nargin<12)
    itcorr=round(dt/dtstep);
end

if(length(xrec)~=size(shot,2))
    error('length(x) must equal size(shot,2)')
end
if(length(xv)~=size(vel,2))
    error('length(xv) must equal size(vel,2)')
end
if(length(t)~=size(shot,1))
    error('length(t) must equal size(shot,1)')
end
if(sum(abs(diff(diff(t))))>small)
    error('t must be regularly sampled')
end
if(sum(abs(diff(diff(xv))))>small)
    error('xv must be regularly sampled');
end
if(max(xv)<max(xrec) || min(xv)>min(xrec))
    error('x range of vel must encompass x range of shot')
end
if(abs(abs(tw(2)-tw(1))-dtstep)>dtstep)
    error('wavelet must be sampled at dtstep');
end

if(floor(itcorr)~=itcorr)
    error('itcorr must be an integer');
end

if(laplacian==1)
    disp('using second order Laplacian')
elseif(laplacian==2)
    disp('using fourth order Laplacian')
else
    error('Invalid value for laplacian')
end

if(length(zrec)==1 && zrec==0)
    zrec=zeros(size(xrec));
end

if(length(zrec)~=length(xrec))
    error('xrec and zrec must be the same size');
end

if(length(xshot)~=length(zshot))
    error('xshot and zshot must be the same size');
end

unstable=0;
dxin=abs(xv(2)-xv(1));

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

%adjust x coordinates so that xv starts at 0
xmin=min(xv);
xv=xv-xmin;
xshot=xshot-xmin;
xrec=xrec-xmin;

zv=dxin*(0:size(vel,1)-1)';
tstep=(t(1):dtstep:t(end))';
% if the input parameters are unstable
if(unstable)
    error('Input parameters are unstable')
%     %need to intepolate vel to a coarser spatial grid
%     if(laplacian==1)
%         dxnew=ceil(vmax*dtstep*sqrt(2));
%     else
%         dxnew=ceil(vmax*dtstep*sqrt(8/3));
%     end
%     %want (x(end)-x(1))/dxnew to be an integer
%     nx=ceil(abs((xv(end)-xv(1)))/dxnew);
%     dxnew=(xv(end)-xv(1))/nx;
%     disp(['Resampling spatial grid from ' num2str(dxin) ' to ' num2str(dxnew)])
%     xnew=xv(1):dxnew:xv(end);
%     znew=(zv(1):dxnew:zv(end))';
%     velnew=interp2(xv(ones(size(zv)),:),zv(:,ones(size(xv))),vel,xnew,znew);
%     vel=velnew;
%     xv=xnew;
%     clear velnew xnew znew
end

%interpolate the input shot record in time to the sample rate of the time stepping
if(dtstep~=dt)
    %interpolate seismic to dtstep
    shotnew=zeros(length(tstep),length(xrec));
    for k=1:length(xrec)
        shotnew(:,k)=interpbl(t,shot(:,k),tstep);
    end
    shot=shotnew;
    clear shotnew
end


if(any(zshot<zv(1)) || any(zshot>zv(end)))
   error('shot depth is outside the model') 
end

if(any(xshot<xv(1)) || any(xshot>xv(end)))
   error('shot location is outside the model') 
end

if(any(zrec<zv(1)) || any(zrec>zv(end)))
   error('receiver depths are outside the model') 
end

%define some constants
dx=abs(xv(2)-xv(1));
nsteps=length(tstep);
nx=length(xv);
nz=length(zv);
tmax=t(end);
nw=length(w);

%set up cell arrays for output snapshots
if(~isempty(tout))
    tout=sort(tout,'descend');
    tout=dtstep*round(tout/dtstep);
    recsnaps=cell(size(tout));
    sousnaps=recsnaps;
    refsnaps=recsnaps;
    iout=nsteps-round(tout/dtstep);
    kout=0;
    iout2=zeros(size(iout));
else
    iout=[];
end

% prepare for forward time propagation loop
disp(['Forward time loop: There are ' int2str(nsteps) ' steps to complete']);

%set up 3D matrix for saveed snapshots to be used in correlation. These are different from the
%output snapshots
icorr=1:itcorr:nsteps;%these are the times at which we will crosscorrelate
source=zeros(nz,nx,length(icorr));

%transform receiver locations to bin locations
ixrec = floor(xrec./dx)+1;
izrec = floor(zrec./dx)+1;

%determine linear addresses for receivers
irec=(ixrec-1)*nz + izrec;

%transform source locations to bin locations
ixshot = floor(xshot./dx)+1;
izshot = floor(zshot./dx)+1;

%determine linear addresses for shots
ishot=(ixshot-1)*nz + izshot;

%determine time shift to allow for non-causal wavelets
inot=find(tw==0);
if(isempty(inot))
    error('wavelet does not have a sample at time 0');
end

%build initial snapshots
snap1=zeros(size(vel));
snap2=snap1;
snap2(ishot)=w(1);

%save time zero from snap2
source(:,:,1)=snap2;

maxstep=round(tmax/dtstep)-1+(inot-1);
disp(['There are ' int2str(maxstep) ' steps to complete']);
time0=clock;

% begin forward time propagation loop
nwrite=2*round(maxstep/100)+1;
isource=2;%counter for saved source snaps
soumax=0;
for k=1:2:maxstep% each loop does two time steps
	
	%first time step
	snap1=afd_snap(dx,dtstep,vel,snap1,snap2,laplacian,boundary);
	%shot(k+1,:)=snap1(irec);%pull out samples at receivers
    if((k+1)<=nw)
        snap1(ishot)=w(k+1);%inject wavelet
    end
    
    if(rem(k+1,itcorr)==0)
        source(:,:,isource)=snap1;
        isource=isource+1;
        tmp=max(abs(snap1(:)));
        if(tmp>soumax);soumax=tmp;end
    end
	
    %second time step
	snap2=afd_snap(dx,dtstep,vel,snap2,snap1,laplacian,boundary);
	%shot(k+2,:)=snap2(irec);%pull out samples at receivers
    if((k+2)<=nw)
        snap2(ishot)=w(k+2);%inject wavelet
    end
    
    if(rem(k+2,itcorr)==0)
        source(:,:,isource)=snap2;
        isource=isource+1;
        tmp=max(abs(snap2(:)));
        if(tmp>soumax);soumax=tmp;end
    end
    
    
    if rem(k,nwrite) == 0
        timenow=clock;
        tottime=etime(timenow,time0);
        
        disp(['wavefield propagated to ' num2str(k*dtstep) ...
            ' s; computation time left ' ...
            num2str((tottime*maxstep/k)-tottime) ' s']);
    end

end
disp('completed source modelling');

%begin reverse-time propagation loop
nsteps=length(tstep);%number of reverse time steps needed
disp(['Reverse time loop: There are ' int2str(nsteps) ' steps to complete']);
time0=clock;
snap1=zeros(size(vel));%snapshot at time tmax+dtstep
snap2=snap1;%snapshot at time tmax
snap2(irec)=shot(end,:);%inject last row of receiver samples at irec at top of snap2
nsources=size(source,3);%number of saved source snapshots
isource=nsources-1;%this will be the first correlation point
shotm=zeros(nz,nx);
nwrite=2*nwrite;
stab=.01;
for k=1:nsteps
    
    snap3=afd_snap(dx,dtstep,vel,snap1,snap2,laplacian,boundary);
    %snap1 is at time t+dtstep, snap2 is at time t, and snap3 is at time t-dtstep 
    snap1=snap2;%snap2 moves to t+dtstep
    snap2=snap3;%snap3 moves to t
    
    if(k<nsteps)
        %inject the next seismic samples from the shot record
        snap2(irec)=snap2(irec)+shot(end-k,:);
    end
    
    if(rem(k,itcorr)==0 && isource>0)
        %correlate course and receiver fields and add to reflectivity
        if(illume==1)
            %gained by illumination
            shotm=shotm+snap2.*source(:,:,isource)./(source(:,:,isource).^2+stab*soumax^2);
        else
            %simple correlation
            shotm=shotm+snap2.*source(:,:,isource);
        end
        isource=isource-1;   
    end
    
    if rem(k,nwrite) == 0
        timenow=clock;
        timeused=etime(timenow,time0);
        timeremaining=(nsteps-k)*timeused/k;
        
        disp(['Completed ' num2str(k) ' reverse-time steps in ' int2str(timeused) ' s. Time remaining ' int2str(timeremaining) ' s'])
    end
    
    
    %save snapshots if requested
    ind=find(k==iout, 1);
    if(~isempty(ind) && isource>-1)
        kout=kout+1;
        recsnaps{kout}=snap3;
        sousnaps{kout}=source(:,:,isource+1);
        refsnaps{kout}=shotm;
        iout2(kout)=iout(ind);
    end
    
end

if(~isempty(tout))
    xms=xv;
    zms=zv;
else
    recsnaps=[];
    sousnaps=[];
    refsnaps=[];
    xms=[];
    zms=[];
end

disp('RTM completed')

%downsample to output x and z
if(dxin>dx)
    ievery=round(dxin/dx);%should be an integer
    shotm=shotm(1:ievery:end,1:ievery:end);
end
xm=xv(1)+(0:size(shotm,2)-1)*dxin+xmin;
zm=(0:size(shotm,1)-1)*dxin;
%apply laplacian filter
shotmf=del2(shotm);


