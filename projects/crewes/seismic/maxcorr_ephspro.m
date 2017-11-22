function [x,obj]=maxcorr_ephspro(seis,sr,n)
% MAXCORR_EPHS: determine best shift and phase rotation between two traces
% 
% [x,str1,str2]=maxcorr_ephspro(seis,sr,n)
%
% MAXCORR_EPHSPRO compares a profiles of traces with a reference trace. Each trace in the profiles
% is compared individually to the references trace. The comparison cross correlates the envelopes of two given
% traces to find the maximum correlation and its lag. Then trace2 (the reference) is shifted by the discovered lag
% and the best constant phase rotatation is estimated between trace1 and the shifted trace2. The
% returned values are averages for the profile. The value of the maximum correlation, its lag, and
% the phase rotation (in degrees) are returned as x(1) x(2) and x(3). The maximum correlation and
% lag after shifting and rotating are returned in x(4) and x(5).
%
% See also: maxcorr_phs and maxcorr
% 
% NOTE: to adjust sr to look like the average of seis, do
%       sradj=phsrot(stat(sr,t,dt*x(2)),x(3))
%
% seis= gather of traces to be compared to reference trace
% sr= reference trace 
% NOTE: seis and sr should be the same length (required to determine
%       the phase rotation)
% n= 2*n +1 lags will be computed
% ******* default= round(length(sr)/10) *********
%
% x= output: x(1)-> interpolated maximum cross correlation (between signals)
%            x(2)-> interpolated lag of maximum correlation (between envelopes)
%            x(3)-> best constant phase rotation (degrees) after shifting
%            x(4)-> max cc (between signals) after shift and rotation
%            x(5)-> lag of max cc after shift and rotation
% Note: a negative result for x(2) indicates trace2 is delayed
%       relative to trace 1
%
% by G.F. Margrave, Aug 2016
%

nsamps=length(seis(:,1));
if(length(sr)~=nsamps)
    error('seis and sr must be the same length');
end
if(nargin<3)
    n=round(length(sr/10));
end

ntraces=size(seis,2);
cc=zeros(ntraces,2);
cce=cc;
cc2=cc;
phs=zeros(1,ntraces);
obj=zeros(1,360);
nobs=0;
for k=1:ntraces
    test=sum(abs(seis(:,k)));
    if(test>0)
        e1=env(seis(:,k));
        e2=env(sr);
        cc(k,:)=maxcorr(seis(:,k),sr,n);
        cce(k,:)=maxcorr(e1,e2,n);
        
        dt=.001;%pretend we are at 1 mil. It does not matter.
        t=dt*(0:nsamps-1)';
        sr2=stat(sr,t,dt*cce(2));
        
        [phs(k),tmp]=constphase3(sr2,seis(:,k));
        nobs=nobs+1;
        obj=obj+tmp;
        sr2r=phsrot(sr2,phs(k));
        
        cc2(k,:)=maxcorr(seis(:,k),sr2r);
    end
end
obj=obj/nobs;
ind=find(cc(:,1)~=0);
x=[mean(cc(ind,1)) mean(cce(ind,2)) mean(phs(ind)) mean(cc2(ind,1)) mean(cc2(ind,2))];

