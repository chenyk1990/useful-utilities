function [shotmigdec,shotmigcc,illumination,chkpnts]=pspi_shot(shot,t,x,vel,xv,zv,xshot,frange,stab,zout,w,tw)
% PSPI_SHOT: shot record migration by the PSPI algorithm
%
% [shotmigdec,shotmigcc,illumination,chkpnts]=pspi_shot(shot,t,x,vel,xv,zv,xshot,frange,stab,zout)
%
% PSPI_SHOT performs 2D depth migration for PP events on a single shot record using the PSPI
% algorithm (Gazdag and Squazerro 1984) for wavefield extrapolation. Required inputs are the shot
% record, which must be regularly sampled in both x and t, and a p-wave instantaneous velocity
% model. The velocity model is a depth matrix and the migrated shot will have the same dimensions as
% the velocity model. Two migrated shot records are returned corresponding to a cross-correlation
% imaging condition and a deconvolution imaging condition. If the former is to be used, then the
% shot record should either be gained before migration or after (see GAINCC). The technical content
% of this code is the same as PSPI_MIG but it is repackaged here to be easier to use.
%
% Prestack migration in this case is accomplished by the simultaneous depth stepping (downward
% continuation) of two wavefields: (1) the input shot record and (2) a model (generated internally)
% of the source wavefield. The downward continuation is done by one-way phase shift so no multiples
% are generated or handled. At each depth the two wavefields are compared in the x-f
% (space-frequency) domain using either the correlation imaging condition (CCIC) or the stabilized
% deconvolution imaging condition (DECIC). Let wrec be the "receiver wavefield" (e.g. the input
% shot) in the x-f domain and wsou be the modelled source wavefield also in the x-f domain and both
% are at some specific depth z. Then the reflectivity estimated at that depth through CCIC is
% rcc=sum(wrec.*conj(wsou)) where conj means complex conjugate and the sum is over frequency.
% Similarly, for DECIC the reflectivity estimate is rdec=sum((wrec.*conj(wsou))./(illum+stab*imax))
% where illum=wsou.*conj(wsou) is the illumination, imax is its maximum absolute value, and stab is
% a small positive number between 0 and 1. The source simulation is done by seeding the 2D
% free-space Green's function into the computation at the first depth step (see greenseed2) and
% optionally applying an input wavelet.
%
% Please see 'Ferguson_Margrave_2005.pdf' or 'Ferguson and Margrave, 2005, Planned seismic imaging
% using explicit one-way operators, Geophysics, V70, NO5, S101 - S109' for details.
%
% shot ... shot record stored as a matrix, one trace per column
% t ... time coordinate vector for shot
% x ... x coordinate vector for shot (x should be regularly sampled).
% NOTE: size(shot) must equal [length(t), length(x)]
% NOTE: the x sample interval for data and velocity model should be the
% same. Ideally, this means that both should be sampled at 1/2 the geophone
% spacing. This will usually require trace interpolation for the shot.
% vel...velocity model. A matrix in consistent units.
% xv ... x (column) coordinate vector for velocity model
% zv ... z (row) coordinate vector for velocity model
% NOTE: size(vel) must equal [length(zv), length(xv)]
% NOTE: The velocity model depth coodinate defines the depth step:
%      dz=zv(2)-zv(1). Currently zv must be regularly sampled.
%      Hence the depth coordinate of the migrated shot will be zv. The x
%      coordinates of the shot record must be contained within the x
%      coordinate span of the velocity model. If the shot does not span the
%      entire velocity model then the shot will be automatically padded
%      with zero traces. If a larger xpad is needed, then do it before migration. 
%      Traces are automatically padded in time to minimize operator
%      wrap-around and the pad is re-zero'd every 10 steps.
% xshot ... x coordinate of shot (a scalar). 
% frange... two element vector giving the frequency range (min and max) to
% be migrated. Example, migrate all frequencies up to 60 Hz: frange=[0 60];
%  ****** default is all frequencies ******
% NOTE: runtime is linearly proportional to the number of frequencies to be
% migrated. There is nothing to be gained by migrating noise or extremely
% low amplitude frequencies.
% stab ... stability factor used in stabilized deconvolution imaging
%   condition. Must be contained in [0,1] .
%  ****** default is .0001 ******
% zout ... array of checkpoint depths at which intermediate results will be output.
%  ****** default is [] (i.e. no checkpoints) ******
% w ... source wavelet to be used in the forward modelling.
% tw ... time coordinate vector for w
% NOTE: w can be causal or acausal but must be sampled at the same rate as the shot record. If not
% provided, the the wavelet is modelled as an impulse.
%
% shotmigdec...depth migrated output using a stabilized deconvolution imaging condition.
% shotmigcc...depth migrated output using crosscorrelation imaging condition.
% illumination ... illum or shot strength at each image point.
% NOTE: Roughly, shotmigdec=shotmigcc./(illumination+stab*max(illumination(:))).
%   This is approximate because the actual calculation is done in the
%   frequency domain (search for variables rcc and rdec in the code to see
%   how its done). Division of rcc by the stabilized illumonation in the
%   frequency domain is not quite the same as doing it in the time domain,
%   but the results are similar. In doubt? Try it yourself.
% NOTE: Both shotmigcc and shotmigdec are the same size as the input
%   velocity model.
% chkpnt ... cell array of length 5 with the following content
%       chkpnt{1} ... array zout of the depths at which checkpoint were written. May not equal zout on input.
%       chkpnt{2} ... cell array of receiver fields extrapolated to the depths in chkpnt{1}
%       chkpnt{3} ... cell array of source fields extrapolated to the depths in chkpnt{1}
%       chkpnt{4} ... cell array of cc depth computations
%       chkpnt{5} ... cell array of dec depth computations
% NOTE chkpnt{2} and chkpnt{3} will be the size of velocity model in columns and the input shot in rows
%      chkpnt{4} and chkpnt{5} will be the size of velocity in rows and columns
% NOTE the chkpnt information can be conveniently examined with plotchkpnts_pspi.
%
% R. J. Ferguson, 2009
% G.F. Margrave 2011-2017
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

%assign defaults
if(nargin<10)
    zout=[];
end
if(nargin>8)
    if(stab<0 || stab>1)
        error('stab must lie in [0,1]');
    end
else
    stab=.0001;
end
if(nargin<7)
    frange=[0 inf];
end

%***get sizes of things***
[Nz,Nx]=size(vel);
dz=zv(2)-zv(1);
[nt,nx]=size(shot);
dt=t(2)-t(1);
small=1e-04;
dx=xv(2)-xv(1);
dxs=x(2)-x(1);
%*************************
%test shot dimensions
if(length(x)~=nx)
    error('shot x axis incorrect');
end
if(length(t)~=nt)
    error('shot t axis incorrect');
end
%test vel dimensions
if(length(zv)~=Nz)
    error('velocity z axis incorrect');
end
if(length(xv)~=Nx)
    error('velocity x axis incorrect');
end
if(abs(dxs-dx)>small)
    error('velocity model and shot must have same x sample size')
end
%determine if the velocity model spans the shot
xmins=min(x);xmaxs=max(x);
xminv=min(xv);xmaxv=max(xv);
if(xmins<xminv || xmaxs>xmaxv)
    error(['Shot x coordinates fall outside the span of the velocity model.'...
        ' Velocity model must be extended']);
end


%test to see if shot and velcoty are on same grid
if(abs(dx*floor((xmaxs-xmins)/dx)-(xmaxs-xmins))>small)
    error('Velocity model and shot record are not on the same x grid');
end
%test for regular shot sampling
if(sum(abs(diff(diff(x))))>small)
    error('Shot record must be regularly sampled in x');
end


%now, pad the shot with zero traces if needed
npadmin=0;
npadmax=0;
if(xmins>xminv)
    npadmin=round((xmins-xminv)/dx);
end
if(xmaxs<xmaxv)
    npadmax=round((xmaxv-xmaxs)/dx);
end
if(npadmin+npadmax>0)
    shot=[zeros(nt,npadmin) shot zeros(nt,npadmax)];
end
x=xv;
nx=length(x);
%pad the shot out so that number of traces is a power of 2
nx2=2^nextpow2(nx);

if(nx<nx2)
    shot=[shot zeros(nt,nx2-nx)];%pad shot with zeros
    vel=[vel vel(:,end)*ones(1,nx2-nx)];%extend vel with last trace
end
x2=(0:nx2-1)*dx;
        
%determine the temporal pad
%pad by enough to hold 10 dz steps
tmax=max(t);
vmin=min(vel(:));
tpad=2*10*dz/vmin;%twice the vertical travel time for 10 steps at slow velocity
npad=round(tpad/dt);
npow=nextpow2(nt+npad);
ntnew=2^npow;%make the new length a power of 2
npad=ntnew-nt;
%pad the traces in time
shot=[shot;zeros(npad,nx2)];
t=dt*(0:ntnew-1);

%fk transform the shot record
[shotfk,f,k]=fktran(shot,t,x2); %#ok<ASGLU>
shotfk=fftshift(shotfk,2);%pspi_ips wants a wrapped wavenumber spectrum

if(frange(1)==0)
    frange(1)=f(2);%don't use 0 Hz
end
if(frange(2)>f(end))
    frange(2)=f(end);
end

%[nf,cd]=size(shotfk);
nf=length(f);
indf=near(f,frange(1),frange(2));%frequencies to migrate
nf2=length(indf);

%***build the source***
%design a secret 85 degree dip limit on energy from the shot. This is done
%to attenuate high angle noise from the shot model. It is a very good
%thing.
xlim=dz*tand(85);
nwindow=round(1.2*2*xlim/dx)+1;%size of a spatial window to be applied to the source
nx0=round(xshot/dx)+1;
nwin2=round((nwindow-1)/2);
mw=mwindow(nwindow,10)';
if((nx0-nwin2)<1)
    %here the window hangs off the left edge
    nskip=nwin2-nx0;%number of samples to skip at beginning of window
    window=[mw((nskip+1):end) zeros(1,nx2-(nwindow-nskip))];
elseif(nx0+nwin2>nx2)
    %here the window hangs off the right edge
    nlast=nwindow-(nx0+nwin2-nx2);
    window=[zeros(1,nx2-nlast) mw(1:nlast)];
else
    %here the window is fully within the span of the model
    nbegin=nx0-nwin2;
    window=[zeros(1,nbegin-1) mw zeros(1,nx2-nwindow-nbegin+1)];
end
temp=zeros(nf,nx2);
%build wavelet
if(nargin<11)
    W=ones(size(f));
else
    if(nargin<12)
        error('tw must be provided when w is provided')
    end
    if(abs(abs(tw(2)-tw(1))-dt)>small)
       error('wavelet must be sampled at the same rate as the shot record'); 
    end
    w2=pad_trace(w(:),t);
    npadw=length(w2)-length(w);
    tw2=[tw(:);tw(end)+dt*(1:npadw)'];
    izero=near(tw2,0);
    w2=[w2(izero:length(w2)); w2(1:izero-1)];
    W=fftrl(w2,tw2);
end
for j=indf
	temp(j,:)=window.*greenseed2(W(j),dx*(0:nx2-1),xshot,f(j),f(end),vel(1,:),dz,1);
end
sourcefk=ifft(temp,[],2);
%**********************
%build the piecwise constant velocity model by the Bagaini method
%the first two input parameters are a mystery but seem to work
vel_blocked=Bagaini(length(x)-1,10,vel);

%prepare for checkdepths
if(~isempty(zout))
    zout=sort(zout);
    ind=find(zout==0);
    if(~isempty(ind))
        zout(ind)=[];
    end
    izout=round(zout/dz);
    reccheck=cell(size(izout));
    soucheck=reccheck;
    deccheck=reccheck;
    cccheck=reccheck;
    izout2=zeros(size(zout));
    icheck=1;
else
    izout=[];
end

%allocate arrays
shotmigcc=zeros(Nz,nx2);
shotmigdec=zeros(Nz,nx2);
illumination=zeros(Nz,nx2);

time1=clock;
timeused=0.0;
ievery=25;
for j=1:Nz-1
    if((rem(j,ievery)-2)==0)
        disp([' pspi prestack mig working on depth ',num2str(j),' of ',num2str(Nz),...
            ' time left ~ ' int2str(timeremaining) '(s)'])
%     else
%         disp([' pspi prestack mig working on depth ',num2str(j),' of ',num2str(Nz)])
    end
    if(rem(j,10)==0)
        shotfk=ps_rezero(shotfk,f,dx,tmax);
        sourcefk=ps_rezero(sourcefk,f,dx,tmax);
    end
    %step the data down. pspi_ips does a domain change from (k,f) to (x,f)
	ftemp=pspi_ips(shotfk(indf,:),f(indf),dx,vel(j,:),vel_blocked(j,:),dz); %ftemp is in (x,f) domain
    %step the source model down
	stemp=pspi_ips(sourcefk(indf,:),f(indf),dx,vel(j,:),vel_blocked(j,:),-dz); %stemp is in (x,f) domain
    %ftemp and stemp are in the (x,f) domain.
%     if(j==141)
%         disp('Break')
%     end
    %imaging conditions
	rcc=ftemp.*conj(stemp);%trivial reflectivity estimate
    illum=stemp.*conj(stemp);%illumination is the shot power
    rdec=rcc./(illum+stab*max(abs(illum(:))));%stabilized decon reflectivity estimate
    %At this point, rcc and rdec are frequency and wavenumber dependent
    %Sum rcc and rdec over temporal frequencyn to get the final estimates
	shotmigcc(j+1,:)=real(sum(rcc)+sum(rcc(1:nf2-1,:)))/(2*nf2-1)/2/pi;
    shotmigdec(j+1,:)=real(sum(rdec)+sum(rdec(1:nf2-1,:)))/(2*nf2-1)/2/pi;
    illumination(j+1,:)=real(sum(illum)+sum(illum(1:nf2-1,:)))/(2*nf2-1);
    %save checkpoints
    if(~isempty(zout))
        if(j==izout(icheck))
            izout2(icheck)=1;%flag to indicate an output has happened
            [rtmp,tout]=back2time(ftemp,f,indf);
            ind=near(tout,0,tmax);
            stmp=back2time(stemp,f,indf);
            reccheck{icheck}=rtmp(ind,1:nx);
            soucheck{icheck}=stmp(ind,1:nx);
            deccheck{icheck}=shotmigdec(:,1:nx);
            cccheck{icheck}=shotmigcc(:,1:nx);
            icheck=min([icheck+1,length(izout)]);
        end
    end
    %transform back to (k,f) domain in preparation for next step
	shotfk(indf,:)=ifft(ftemp,[],2);
	sourcefk(indf,:)=ifft(stemp,[],2);
    timenow=clock;
    timeused=etime(timenow,time1)+timeused;
    time1=timenow;
    timeremaining=(Nz-1)*timeused/j-timeused;
    %disp([' elapsed time ',int2str(timeused),' (s), estimated time remaining '...
    %    ,int2str(timeremaining),' (s)']);
end
shotmigcc=shotmigcc(:,1:nx);
shotmigdec=shotmigdec(:,1:nx);
illumination=illumination(:,1:nx);
disp(['shot migrated in ' int2str(timeused) '(s)'])

if(~isempty(zout))
    chkpnts=cell(1,5);
    chkpnts{1}=zout;
    chkpnts{2}=reccheck;
    chkpnts{3}=soucheck;
    chkpnts{4}=cccheck;
    chkpnts{5}=deccheck;
end
    

function [seist,t]=back2time(seisf,f,indf)
nx=size(seisf,2);
nf=length(f);
tmp=zeros(nf,nx);
tmp(indf,:)=seisf;
[seist,t]=ifftrl(tmp,f);