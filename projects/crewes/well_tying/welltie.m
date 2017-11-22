function [stie,sreftie,tstretch,tphase]=welltie(s,t,sref,tref,twin,tinc,flag,envflag,padflag,bwflag)
% WELLTIE: perform well tying by time-variant stretch and phase rotation
%
% [stie,sreftie,tstretch,tphase]=welltie(s,t,sref,tref,twin,tinc,flag,envflag,padflag,bwflag)
%
% Given a fully processed seismic trace and a reference trace that is either
% a synthetic seismogram or a reflectivity produced from well logs, the
% seismic trace is "tied" to the reference by measuring time-variant delays
% and time-variant phase rotations. It is important to do a rough
% time-alignment of the reference trace prior to running this. This usually
% amounts to simply changing the start time of the reference trace by a
% constant time shift that represents the unlogged overburden. The well
% tying process here consists of mesuring time-variant time-shifts and
% time-variant phase rotations. Once measured, the default action is to
% apply the time-shifts to the reference trace (since they are time
% variant, this is a stretch) and the phase rotations to the seismic trace.
% Optionally, both can be applied to the seismic trace.
% Method: 
%   (1) Transfer the seismic bandwidth to the reference trace (see bandwidth_xfer)
%   (2) Check for polarity and flip s if necessary
%   (3) Measure time-variant time shift by tvmaxcorr using trace envelopes.
%       Using trace envelopes here helps to avoid the usual coupling
%       between time-shifts and phase rotations.
%   (4) Apply the measured time shifts to either the reference trace (the
%       default) or the seismic trace.
%   (5) Measure time-variant phase rotations using tvconstphase.
%   (6) Apply phase rotations to the seismic trace.
% It is expected that the reference trace will be shorter than the seismic
% trace. However, the longer the reference trace the better. If the
% reference trace is less than half the length of the seismic, consider
% creating a composite log using several wells with different depth ranges.
% Also be aware that shifts and rotations for times outside the range of
% the reference trace are simply constant extrapolations.
%
% s ... input seismic trace to be "tied"
% t ... time coordinate for s (same size as s)
% sref ... reference trace (seismogram at well with zero phase wavelet)
% tref ... time coordinate for sref (same size as sref)
% twin ... half width of Gaussian window used for time variant localization. The maximum 
%       expected time shift should be less than twin. A common value would be 0.1 or 0.2.
%       Optionally, this can be a vector of length 2 where the second entry is the length of a
%       boxcar smoother applied to the trace envelopes. If not specified, then the smoother is
%       calculated as 0.1*twin(1).
% tinc ... increment between windows. Typically tinc would be about twin/4.
%       A smaller tinc can detect more rapid time variance. Optionally, this can be a vector of
%       length 2 where the second value gives the maximum allowed cc lag. If not specified this
%       value is 0.4*twin(1).
% flag ... determines how phase rotations and times shifts are applied
%           0: Don't apply phase rotations. Time shifts still applied to sref in
%               order to measure phase.
%           1: Apply phase rotations to s and time shifts to sref
%           2: Apply both phase rotations and time shifts to s
% ******** default flag=1 ********
% envflag ... 0 means the delays are estimate by cross correlating traces
%             1 means the delays are estimate by cross correlating envelopes
% ******** default envflag=1 *********
% padflag ... 0 means the adjusted reference trace is the same size as sref
%             1 means the adjusted reference trace is the same size as s
% ******** default padflag=0 *********
% bwflag ... 0 means do not transfer the bandwidth from s to sref
%            1 means transfer the bandwidth from s to sref
% ******** default bwflag=1 *******
% stie ... adjusted seismic trace
% sreftie ... adjusted reference trace
% tstretch ... measured and applied time shifts (same size as t)
% tphase ... measured and applied phase rotations (same size as t)
%
% To apply the time shifts to another trace, say s2, use
%   s2s=stretcht(s2,t,-tstretch); if shifts are applied to the reference trace and s2 is another reference trace
%   s2s=stretcht(s2,t,tstretch); if shifts are applied to the seismic trace and s2 is another seismic trace
% To apply the phase rotations use
%   s2sr=tvphsrot(s2s,t,-tphase,t);
% Note the minus sign in the phase rotations and that the rotations are applied after shifting.
% This is the case if flag=2 was used. For flag=1, then the reference trace has been shifted
% and the phase rotations should be applied directly to s2 with the opposite sign.
%  s2sr=tvphsrot(s2,t,tphase,t);
%
%
% by: G.F. Margrave, Devon Canada, 2016
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
    bwflag=1;
end

if(nargin<9)
    padflag=0;
end

if(nargin<8)
    envflag=1;
end

if(nargin<7)
    flag=1;
end

dt=t(2)-t(1);
dtr=tref(2)-tref(1);
if(abs(dt-dtr)>1000*eps)
    error('s and sref have different sample rates')
end

if(length(twin)==1)
    twin=[twin .1*twin];
end

if(length(tinc)==1)
    tinc=[tinc 0.4*twin(1)];
end

extrapbeg=0;
extrapend=0;
if(tref(1)>t(1))
    extrapbeg=1;
end
if(tref(end)<t(end));
    extrapend=1;
end
ind=near(t,tref(1),tref(end));

%transfer the seismic bandwidth to the reference
if(bwflag==1)
    df=1/(t(end)-t(1));
    n=round(5/df);
    sref=bandwidth_xfer(s(ind),sref,n);
end

%check for an overall polarity flip
cc=maxcorr(sref,s(ind));
if(cc(1)<0)
    s=-s;
end

%measure time variant cc of envelopes
if(envflag==1)
    %note the use of env here. it fails with traces themselves
    nsmo=round(twin(2)/dt)+1;
    es=convz(env(s(ind)),ones(nsmo,1))/nsmo;
    er=convz(env(sref),ones(nsmo,1))/nsmo;
    %[cc,tcc]=tvmaxcorr(er,es,t(ind),twin,tinc,.4*twin);
    [cc,tcc]=tvmaxcorr(er,es,t(ind),twin(1),tinc(1),tinc(2));
else
    %note the use of traces themselves
    [cc,tcc]=tvmaxcorr(sref,s(ind),t(ind),twin(1),tinc(1),tinc(2));
end
if(extrapbeg==1)
    cc=[cc(1,:);cc];
    tcc=[t(1);tcc];
end
if(extrapend==1)
    cc=[cc;cc(end,:)];
    tcc=[tcc;t(end)];
end
dt=t(2)-t(1);
delt2=cc(:,2)*dt;
delt2=findmonostretch(tcc,delt2,cc(:,1));
tstretch=interp1(tcc,delt2,t);
% remove estimated delay
t1=tref(1);%first live reference sample
t2=tref(end);%las live reference sample
if(padflag)
    srefp=zeros(size(s));
    srefp(ind)=sref;
    trefp=t;
else
    srefp=sref;
    trefp=tref;
end
ind=near(t,trefp(1),trefp(end));
if(flag==2)
    s2=stretcht(s,t,tstretch);
    sreftie=srefp;
elseif(flag==1 || flag==0)
    sreftie=stretcht(srefp,trefp,-tstretch(ind));
    s2=s;
else
    error('invalid value for flag');
end
%adjust t1 and t2
if(padflag==1)
   i1=near(t,t1);
   t1=t1-tstretch(i1);
   i2=near(t,t2);
   t2=t2-tstretch(i2);
end

ind=near(t,t1,t2);
ind2=near(trefp,t1,t2);

%estimate the phase
[phs,tphs]=tvconstphase(sreftie(ind2),s2(ind),t(ind),twin(1),tinc(1));
extrapbeg=0;
extrapend=0;
if(tphs(1)>t(1))
    extrapbeg=1;
end
if(tphs(end)<t(end));
    extrapend=1;
end
if(extrapbeg==1)
    phs=[phs(1);phs];
    tphs=[t(1);tphs];
end
if(extrapend==1)
    phs=[phs;phs(end)];
    tphs=[tphs;t(end)];
end
tphase=interp1(tphs,phs,t);
%remove estimated phase
if(flag~=0)
    stie=tvphsrot(s2,t,-tphase,t);
else
    stie=s2;
end

function deltnew=findmonostretch(tcc,delt,ccmax)
% 
% test and possibly modify the delt to ensure monotonic solution
%
% tcc ... vector of measurment times
% delt ... vector of estimted shifts
% ccmax ... vector of estimated max cc's
% all three are the same length
%
% Plan: The stretch maps t to tp=t+delt. The solution is monotonic if diff(tp)>0 everywhere.
% So, I hunt for those places where diff(tp)<0 and fudge them. Fudging is done by looking at
% the problem points in conparison to the point before. The ccmax values at each
% problem point and the one before it are compared and the point with the lower ccmax is
% discarded and a new value interpolated in.
%

%first test to see if anything needs to be done
tp=tcc+delt;
ind=find(diff(tp)<0);%ideally diff(tp) should always be positive
deltnew=delt;
if(isempty(ind))
    return;
end
npts=length(tcc);
count=0;
while(~isempty(ind))
    count=count+1;
    if(count>npts)
        error('failed to converge to single-valued stretch solution in welltie');
    end
    for k=1:length(ind)
        % compare the points at ind(k) and ind(k)+1; The second is the low one
        i1=ind(k);
        i2=i1+1;
        if(ccmax(i1)>ccmax(i2))%we favour the point with the higher ccmax
            %discard all points where tp<tp(i1) and that have tcc>tcc(i1) and interpolate a replacement
            ii=find(tp<tp(i1));%first condition
            %iii=find(tp>tp(i1));
            %i3=iii(1);
            i3=ii(end)+1;
            for kk=1:length(ii) %check each point
                if(tcc(ii(kk))>tcc(i1))%second condition
                    i2=ii(kk);%point to be replaced
                    if(i3<=npts)
                        %linear interpolation between i1 and i3
                        deltnew(i2)=deltnew(i1)*(tcc(i2)-tcc(i3))/(tcc(i1)-tcc(i3))+deltnew(i3)*(tcc(i2)-tcc(i1))/(tcc(i3)-tcc(i1));
                    else
                        %there is no i3
                        deltnew(i2)=deltnew(i1);
                    end
                end
            end
        else
            %discard all points where tp>tp(i2) and that have tcc<tcc(i2) and interpolate a replacement
            ii=find(tp>tp(i2));
            i0=ii(1)-1;
            for kk=1:length(ii)
                if(tcc(ii(kk))<tcc(i2))%second condition
                    i1=ii(kk);
                    if(i0>0)
                        %linear interpolation between i0 and i2
                        deltnew(i1)=deltnew(i0)*(tcc(i1)-tcc(i2))/(tcc(i0)-tcc(i2))+deltnew(i2)*(tcc(i1)-tcc(i0))/(tcc(i2)-tcc(i0));
                    else
                        %there is no i0
                        delt(i1)=delt(i2);
                    end
                end
            end
        end
    end
    tp=tcc+deltnew;
    ind=find(diff(tp)<0);
end
                
    