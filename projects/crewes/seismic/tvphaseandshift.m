function [phs,tphs,delay,ccs,ccs0]=tvphaseandshift(s1,s2,t,twin,tinc)
% TVPHASEANDSHIFT: estimates apparent temporally local constant phase rotations and time shifts 
%
% [phs,tphs,delay,ccs,ccs0]=tvphaseandshift(s1,s2,t,twin,tinc)
% 
% The input parameters t,twin, and tinc define a discrete set of Gaussian windows of half-width twin
% and separated by tinc and that span the length of t. Input parameters are altered slightly so that
% the first and last windows are centered on t(1) and t(end). Then, for each window, s1 and s2 are
% localized (by multiplication by the window) to give s1w and s2w and compared.  The first step of
% the comparison is to crosscorrelate the envelopes of s1w and s2w to determine an apparent shift.
% Then this shift is applied to s2w. Then the constant phase rotation that best matches the shifted
% s2s to s1 is determined. Finally crosscorrelations before and after the process are calculated.
%
% s1= input trace to be analyzed
% s2= reference trace. Constant phase rotations are w.r.t. this trace
% t= time coordinate vector for s1
% twin= width (seconds) of the Gaussian window (standard deviation)
% tinc= temporal shift (seconds) between windows
%
% phs= vector of apparent constant phase rotations in each window
% tphs= vector of window center times. Same size as phs
% delay= vector of time shifts of s2 relative to s1. same size as phs
% ccs= vector of maximum cc values calculated after removing delay and phase. same size as phs
% ccs0= vector of maximum cc values calculated before removing delay and phase. same size as phs
%
% NOTE: To correct s1 to look like s2, use 
%         s1d=stretcht(s1,t,delay,tphs);%remove the delays
%         s1r=tvphsrot(s1d,t,phs,tphs);%remove the phase rotations
%
% NOTE2: There is something strange about the ccs and ccs0 values. Don't trust them. Use tvmaxcorr
% and measure them independently.
%
% by G.F. Margrave, Sept 2005-2016
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


if(length(s1)~=length(s2))
    error('s1 and s2 must have the same length');
end
if(length(t)~=length(s1))
    error('s1 and t must have the same length');
end
dt=abs(t(2)-t(1));
tmin=t(1);
t=t(:);
% determine number of windows. tinc will be adjusted to make the
% last window precisely centered on tmax
tmax=t(end);
nwin=(tmax-tmin)/tinc+1; %this will generally be fractional
nwin=round(nwin);
tinc=(tmax-tmin)/(nwin-1); %redefine tinc
tphs=zeros(nwin,1);
phs=tphs;
ccs=tphs;
ccs0=tphs;
delay=tphs;
n=round(twin/dt);%max allowed shift
for k=1:nwin
    %build the gaussian
    tnot=(k-1)*tinc+tmin;
    tphs(k)=tnot;
    gwin=exp(-((t-tnot)/twin).^2)/(sqrt(pi)*twin/tinc);
    %window and measure phase
    s1w=s1.*gwin;
    s2w=s2.*gwin;
    %phs(k)=constphase2(s1w,s2w,flag);
    x=maxcorr_ephs(s1w/max(abs(s1w)),s2w/max(abs(s1w)),n);
    
    ccs0(k)=x(1);
end
