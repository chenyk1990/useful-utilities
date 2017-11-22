function [shotmig,tmig,xmig,timetables,shottest]=kirk_shot_alt(shotrec,t,x,xshot,velmod,tv,xv,params,timeoutinfo,itest)
% KIRK_SHOT: Kirchhoff prestack shot-record time migration
% 
% [shotmig,tmig,xmig,timetables,shottest]=kirk_shot(shotrec,t,x,xshot,velmod,tv,xv,params,timeoutinfo,itest)
%
% KIRK_SHOT migrates a single shot record using prestack time
% migration. The algorithm is simple traveltime path summation with a few
% options as described in the parameters vector. Traveltime from source to
% scatterpoint (i.e the image point) is approximated by a Dix equation
% using the rms velocity from the model at the lateral position halfway
% between source and receiver and at the vertical traveltime of the 
% scatterpoint. Similarly, from the scatterpoint to a receiver, a Dix 
% equation using the rms velocity at halfway between scatterpoint and 
% receiver is used. There is no topographic compensation at present. The
% source and all receivers are assumed to be on the same horizontal plane.
%
% shotrec ... matrix containing the shot record. One trace per column.
% t ... time coordinate vector for shotrec. 
%           Requirement: length(t)=size(shotrec,1);
% x ... space coordinate vector for shotrec.
%           Requirement: length(x)=size(shotrec,2);
% velmod ... velocity model. This is a matrix of RMS velocities as a
%           function of lateral position and time.
% tv ... time coordinate vector for velmod
%           Requirement: length(tv)=size(velmod,1);
% xv ... space coordinate vector for velmod
%           Requirement: length(xv)=size(velmod,2);
% REQUIREMENT: the span of tv and xv must equal or exceed that of t and x.
% xshot ... lateral position of the shot. It must be in the same coordinate
%       system as vector x, the receiver locations. If x was specified as a
%       scaler
%
% params ... vector of length 13 containing migration parameters. An entry of nan gets the
% default.
%
%    params(1--3) : migration aperture and its taper
%    params(1) ... physical aperture in meters. This is the largest lateral
%               distance a trace may have from the output location and still be
%               allowed to contribute.
%	**********   default is the length of the velocity model
%    params(2) ... width of the aperture taper
%   **********   default is 0.05*params(1)
%    params(3) ... = 0, linear taper
%                     = 1, cosine taper   
%   **********   default is 1 (cosine taper)
%    params(4-6) : scattering angle limit (degrees)
%    params(4) ... maximum scattering angle limit (degrees)
%   **********   default = 60
%    params(5) ... width of angle limit taper
%   **********   default = min([0.15*params(4), 90-params(4)])
%    params(6) ... taper type:
%                     = 0: linear taper;
%                     = 1: cosine taper.
%   **********   default = 1.
%    params(7) : relative to sample interpolation
%    params(7) ... = 1, linear interpolation
%                     = 2, cubic interpolation
%                     = 3, spline interpolation
%                     = 4, sinc interpolation
%   **********   default = 1
%    params(8--11) : relative to migration target window
%	 params(8) ... tmin of migration target window
%	**********   default = min(tv) 
%	 params(9) ... tmax of migration target window
%	**********   default = max(tv)
%    params(10) ... xmin of target window
%	**********   default = min(xv)
%    params(11) ... xmax of migration target window
%	**********   default = max(xv) 
%
%    params(12) : box-car anti-aliasing filter
%    params(12) ... = 0, no box-car filter used;
%                      = 1, box-car filter will be used.
%   **********   default is 0.
%    params(13) ... dxmig: spatial sample rate for migrated traces
%   **********   default ... mean(diff(x))/2
%
%    timeoutinfo ... cell array defining travletime output. This is designed to show the summation
%           curves required to produce a given output trace.
%           timeoutinfo={xout,tinc,xinc} xout is the migrated trace location. A traveltime table
%           will be output for this location with traveltince every tinc seconds and xinc distance
%           across the input shot. xout may be a vector in which case tables will be output for each
%           location.
%   **********  default for timeoutinfo is [] ************
%
%    itest ... a trace number for which the impulse response of the migration is desired. If
%           programmed, then shottest will contain the migration of this single trace. 
%   ********** default = [] ********** 
%
% OUTPUT arguments
%
%    shotmig ... the output migrated shot record
%    tmig ... t coordinates of migrated data
%    xmig ... x coordinates of migrated data.
% NOTE: xmig spans a range of x coordinates controlled by params(10:11). By default xmig is the same
% range as the velocity model. This is usually too large so use params(10:11) to limit. The
% increment between migrated traces is params(13). It defaults to 1/2 the input spacing.
%
%    timetables ... summation curve time tables as requested by timeoutinfo. For each xout, a table 
%           of double square root traveltimes is output describing the summation curve for that
%           output location. Each table row is a constant time. Each table column is a constant x
%           position in the input shot. Time and x increments are as defined on input. timetables{1}
%           is a cell array of traveltimes tables, all identically sized. timetables{2} is the row
%           coordinates of the tables, timetables{3} is the column coordintes, timetables{4} is the
%           spatial locations of the tables. Note that timetables{4} will usually be slightly
%           different from the input timeoutinfo{1}. This is because the output locations are
%           shifted to the nearest output trace location.
%
% By G.F. Margrave and J.K. Cooper
% CREWES Project, U of Calgary, 2007-2008
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

if(nargin<9)
    timeoutinfo=[];
    timetables=[];
else
    if(~iscell(timeoutinfo))
        error('timeoutinfo must be a cell array');
    end
    if(length(timeoutinfo)~=3)
        error('timeoutinfo must be a cell array of length 3');
    end
end

if(nargin<10)
    itest=[];
else
    if(~between(1,length(x),itest))
        error('itest is not a valid trace number');
    end
end
if(isempty(itest))
    itest=-1;
end

%tstart=clock;
[nsamp,ntr]=size(shotrec);
[nvsamp,nvtr]=size(velmod);

% check the validity input arguments

if(length(t)~=nsamp)
    error('Incorrect time specification')
end
t=t(:);
dt=t(2)-t(1);

if(length(x)~=ntr)
    error('Incorrect x specification')
end
dx=x(2)-x(1);

%  ---- test velocity info ----
if(length(tv)~=nvsamp)
    error('Time vector for velocity model is incorrect')
end
if(length(xv)~=nvtr)
    error('Space coordinate vector for velocity model is incorrect')
end
if(min(tv)>min(t) || max(tv)<max(t))
    error('Time vector for velocity model must span that for data')
end
if(min(xv)>min(x) || max(xv)<max(x))
    error('Space vector for velocity model must span that for data')
end

%  ---- examine parameters ----
nparams=13; 				% number of defined parameters
					
if(nargin<8) 				% no parameters inputted
	params= nan*ones(1,nparams); 
end 	 	

if(length(params)<nparams) 
	params = [params nan*ones(1,nparams-length(params))];
end

%assign parameter defaults

if( isnan(params(1)) ) 
	aper = abs(max(xv)-min(xv));
else
    aper = params(1);
end

if( isnan(params(2)) )
    width1 = aper/20;
else
    width1 = params(2);
end

if( isnan(params(3)) )
    itaper1 = 1;
else
    itaper1 = params(3);
end

if( isnan(params(4)) )
	ang_limit = pi/3;
else
	ang_limit = params(4)*pi/180;
end

if( isnan(params(5)) )
	width2 = 0.15*ang_limit;
    if(width2>pi/2-ang_limit)
        width2=pi/2-ang_limit;
    end
else
	width2 = params(5)*pi/180;
end
angle1 = ang_limit + width2;

if( isnan(params(6)) )
	itaper2 = 1;
else
	itaper2 = params(6);
end
if itaper2 ~= 1 && itaper2 ~= 0
	error('the angle limit taper type: params(6) should be 0 and 1 !');
end

if( isnan(params(7)) )
	interp_type = 1;
else
	interp_type = params(7);
end
if interp_type < 1 || interp_type > 4
	error('the interpolation indexx paarams(7) should be 1, 2, 3 and 4 !');
end

if( isnan(params(8)) ) 
		tmig1 = min(tv);
else
		tmig1 = params(8);
end

if( isnan(params(9)) ) 
		tmig2 = max(tv);
else
		tmig2 = params(9);
end
if tmig2 < tmig1
	error(['the target time window start time should be smaller than the end time !'...
         ' i.e. paraams(8) < params(9)']);
end

if( isnan(params(10)) ) 
		xmig1 = min(xv);
else
		xmig1 = params(10);
        ix=near(xv,xmig1);
        xmig1=xv(ix);
end
if( isnan(params(11)) ) 
		xmig2 = max(xv);
else
		xmig2 = params(11);
        ix=near(xv,xmig2);
        xmig2=xv(ix);
end
if xmig2 < xmig1
 	error(['the start location of target trace range should be less than the end location'...
        ' i.e. params(10) < params(11)']);
end

if( isnan(params(12)) )
		ibcfilter = 0;
else
		ibcfilter = params(12);
end

if( isnan(params(13)) )
		dxmig=mean(diff(x))/2;
else
		dxmig = params(13);
end

if ibcfilter 
	% get a cumulative array from shotrec
	arycum=cumsum(shotrec);
end

%one way time
dt1=.5*dt;
t1=t/2;

%compute maximum time needed. This is the traveltime for a scatterpoint
%when the source and receiver are colocated a distance aper away.
vmin=min(velmod(:));
tmax=sqrt(tmig2^2 + (2*aper/vmin)^2);

%pad input to tmaxin
npad=ceil(tmax/dt1)-nsamp+5;
if( npad > 0)
	shotrec= [shotrec; zeros(npad,ntr)];
	t1 = [t1',(nsamp+1:nsamp+npad)*dt1]';
    if ibcfilter
        arycum=[arycum; ones(npad,1)*arycum(nsamp,:)];
    end
end
t2= t1.^2;

%output samples targeted
samptarget=near(t,tmig1,tmig2);
tmig=t(samptarget);

%output traces desired
%trtarget= near(xv,xmig1,xmig2);
xmig=xmig1:dxmig:xmig2;

%initialize output array
shotmig=zeros(length(samptarget),length(xmig));
if(itest>0)
    shottest=shotmig;
end

%loop over migrated traces
%
kmig=0;

disp(' ');
disp([' --- Total number of migrated traces : ' int2str(length(xmig)) ' ---']);
disp(' ');

t0=clock;
% steptimes=nan*ones(size(xmig));
% ntimes=0;
ievery=20;%print a progress message every this many traces
%save1=zeros(size(xmig));
% loop over traces in aperture
if(~isempty(timeoutinfo))
    timetables=cell(1,4);
    xout=timeoutinfo{1};
    tinc=timeoutinfo{2};
    xinc=timeoutinfo{3};
    ttable=t(1):tinc:t(end);
    xtable=x(1):xinc:x(end);
    tables=cell(size(xout));
    timetables{2}=ttable;
    timetables{3}=xtable;
    ixout=nan*zeros(size(xout));
    for k=1:length(ixout)
        ind=near(xmig,xout(k));
        ixout(k)=ind(1);
        xout(k)=xmig(ixout(k));
    end
    ittab=round((ttable-t(1))/dt)+1;
    smallx=abs(x(2)-x(1))*.001;
    timetables{4}=xout;
end
for ktr=1:length(xmig) 			% ktr--the index of the output trace
    xtr=xmig(ktr);%x coordinate of target
	kmig=kmig+1; 			% trace counter

	%determine input traces in aperture
    inaper=near(x,xtr-aper,xtr+aper);

	%shot offset and velocity 
    offsetshot2=(xtr-xshot)^2;
    xvshotside=(xshot+xtr)/2;
    ivshot=near(xv,xvshotside);
	vshot2 = velmod(:,ivshot(1)).^2;
    %gather=zeros(length(tmig),length(inaper));
    %check for table output location
    itable=find(abs(xout-xtr)<smallx);
    if(~isempty(itable))
        disp(['table output at xout=' num2str(xout(itable))])
        thistable=nan*zeros(length(ttable),length(xtable));
    end
    for kaper=1:length(inaper)
        xnow=x(inaper(kaper));
        %receiver offset and velocity
        offsetrec2=(xtr-xnow)^2;
        xvrecside=(xnow+xtr)/2;
        ivrec=near(xv,xvrecside);
        vrec2 = velmod(:,ivrec(1)).^2;
        
        % source-receiver travel time via double square root equation
        tsr= sqrt(offsetshot2./vshot2(samptarget) + t2(samptarget))+...
            sqrt(offsetrec2./vrec2(samptarget) + t2(samptarget))+100*eps;
        
        %check for table output
        if(~isempty(itable))
            ixtable=find(abs(xtable-xnow)<smallx);
            if(~isempty(ixtable))
                thistable(:,ixtable)=tsr(ittab);
            end
        end
        
        %cosine theta amplitude correction
        costheta = tmig./tsr;
        tanalpha = sqrt(1-costheta.^2);
        
        %angle limit and the taper
        
        ind = find( costheta < cos(angle1) );
        i1=1;
        if(~isempty(ind))
            i1 = ind(end);
        end
        ind = find( costheta < cos(ang_limit) );
        i2=1;
        if(~isempty(ind))
            i2 = ind(end);
        end
        
        if i1 < i2
            if itaper2  ==  0
                coef2 = lin_taper(i2,i1);
            else
                coef2 = cos_taper(i2,i1);
            end
            costheta(1:i1) = zeros(i1,1);
            costheta(i1+1:i2) = coef2(i2-i1:-1:1)'.*costheta(i1+1:i2);
        end
        
        % boxcar anti-aliasing filter
        if ibcfilter
            lt0=round((dx*tanalpha./velmod(samptarget,ktr)/dt));
            indt = round((tsr/dt))+1;
            lentr = nsamp+npad;
            lt = ones(lentr,1)*max(lt0);
            lt(indt)=lt0;
            lt(max(indt)+1:lentr) = ones(lentr-max(indt),1)*min(lt0);
            it = (1:lentr)';
            l1=it-lt-1;
            l2=it+lt;
            ind = find(l1 < 1);
            l1(ind) = ones(length(ind),1);
            ind = find(l2> lentr);
            l2(ind)=ones(length(ind),1)*lentr;
            tmp0=t;
            tmp0(1) = arycum(1,inaper(kaper));
            ind = 2:lentr;
            tmp0(ind) = (arycum(l2(ind),inaper(kaper))-arycum(l1(ind),inaper(kaper)))...
                ./(l2(ind)-l1(ind));
        else
            tmp0 = shotrec(:,inaper(kaper));
        end
        
        %interpolation
        % Linear
        if interp_type == 1
            tnumber = tsr/dt;
            it0 = floor( tnumber ) + 1;
            it1 = it0+1;
            xt0 = tnumber - it0 + 1;
            xt1 = it0-tnumber;
            tmp = xt1.*tmp0(it0)+xt0.*tmp0(it1);
        end
        % Spline
        if interp_type == 2
            tmp = interp1(t,tmp0,tsr,'spline');
        end
        % Cubic
        if interp_type == 3
            tmp = interp1(t,tmp0,tsr,'cubic');
        end
        % Sinc
        if interp_type == 4
            tmp = sinci(tmp0,t,tsr);
        end
        
        % aperture taper
        aper_weight = 1.0/length(inaper);
        xtest=abs(aper-abs(xtr-xnow));%distance of trace from edge of aper
        if xtest < width1
            if(itaper1==1)
                aper_weight=(.5+.5*cos(pi*(xtest-width1)/(180*width1)))/length(inaper);
            else
                aper_weight=(xtest-width1)/(width1*length(inaper));
            end
        end
        
        tmp = tmp .* aper_weight;
        
        
        tmp = tmp.* sqrt(costheta.^3);
        %save1(kmig)=save1(kmig)+sum(costheta);
        %gather(:,kaper)=tmp;
        
        shotmig(:,kmig)= shotmig(:,kmig)+tmp;
        if(inaper(kaper)==itest)
            shottest(:,kmig)=shottest(:,kmig)+tmp;
        end
        
    end
    
    %check for table output
    if(~isempty(itable))
        tables{itable}=thistable;
    end
	
	% scaling and 45 degree phase shift
    ivmig=near(xv,xtr);
	scalemig = velmod(samptarget,ivmig(1)).*sqrt(pi.*(tmig+0.0001));
	shotmig(:,kmig) = phsrot(shotmig(:,kmig),-45)./scalemig;
    if(itest>0)
        shottest(:,kmig) = phsrot(shottest(:,kmig),-45)./scalemig;
    end
    
    if(rem(kmig,ievery)==0)
	    disp([' Completed migrated trace no. ' ,int2str(kmig) ,' of ' int2str(length(xmig)) ]);
        timenow=clock;
        timeused=etime(timenow,t0);
        time_per_trace=timeused/ktr;
        time_left=(length(xmig)-ktr)*time_per_trace;
        
        disp([' time so far ' int2str(timeused) ' (s) estimated remaining ' int2str(time_left) ' (s)' ]);
    end

end
if(~isempty(timeoutinfo))
    timetables{1}=tables;
end

totaltime=etime(clock,t0);
disp(['Total time required ' num2str(totaltime)])