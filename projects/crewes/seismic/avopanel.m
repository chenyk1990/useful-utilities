function [avopanel,xoff]=avopanel(seis,t,xrec,xshot,frange,mute,xoffmax,cmpmin,cmpmax,dxoff,dxinc,flag1,flag2)
% AVOPANEL ... Sums along midpoint at constant offset to create a gather showing average AVO
%
% [avopanel,xoff]=avopanel(seis,t,xrec,xshot,frange,mute,xoffmax,cmpmin,cmpmax,dxoff,flag1,flag2)
% 
% seis ... cell array of shots along a 2D line
% t ... time coordinate vector
% xrec ... cell array of receiver coordinates for seis. Same size as seis.
% xshot ... ordinary array of shot coordinates for seis. Same size as seis.
% frange ... two element vector giving min and max frequencies (no default). Only frequencies in
%           this range (with a slight boundary taper) will contribute to the panel. 
% mute ... (offset,t) values of a mute to be applied to each trace. Only times larger than the mute
%           time at each offset will contribute to the panel. Specified like mute=[off1 t1;off2
%           t2;....]. Mute must not be double valued (offset must always increase).
% ********** default is [0, 0;xofffmax, 0] (or no mute).
% xoffmax ... maximum offset to examine
% ********** default is all offsets (nan gets default) ***********
% cmpmin ... minimum cmp to contribute to the panel
%  ********* default is all cmp's are used (nan gets default) *************
% cmpmax ... maximum cmp to contribute to the panel
%  ********* default is all cmp's are used (nan gets default) *************
% dxoff ... size of offset bins
%  ********* default is the median geophone spacing on shot#1 (nan gets default) *********
% dxinc ... increment between offset bins (nan gets default)
%  ********* default dxinc = dxoff/2 *********
% flag1 ... 1 means the panel shows the average absolute value, 0 means it shows the average trace
%           (for each offset)
%  ********* default is 1 *************
% flag2 ... 1 means compute offsets as xrec-xshot while 2 means compute offsets as 2*(xrec-xshot).
%           The first option is appropriate for unmigrated data and the second for migrated data.
% ********** default is 1 ***************
%
if(nargin<13)
    flag2=1;
end
if(nargin<12)
    flag1=1;
end
if(nargin<11)
    dxinc=nan;
end
if(nargin<10)
    dxoff=nan;
end
if(nargin<9)
    cmpmax=nan;
end
if(nargin<8)
    cmpmin=nan;
end
if(nargin<7)
    xoffmax=nan; 
end
if(nargin<6)
    mute=nan;
end
if(isnan(cmpmin))
    cmpmin=-inf;
end
if(isnan(cmpmax))
    cmpmax=inf;
end
if(isnan(xoffmax))
    xoffmax=inf;
end
if(isnan(dxoff))
    dxoff=median(abs(diff(xrec{1})));
end
if(isnan(dxinc))
    dxinc=dxoff/2;
end


nshots=length(seis);
if(length(xrec)~=nshots)
    error('xrec is the wrong size');
end
if(length(xshot)~=nshots)
    error('nshot is the wrong size');
end
nt=length(t);
if(size(seis{1},1)~=nt)
    error('t is the wrong size');
end

%determine offset range
xoffm=0;
for k=1:nshots
    xoffs=flag2*abs(xrec{k}-xshot(k));
    xm=max(xoffs);
    if(xm>xoffm); xoffm=xm; end
end
if(xoffmax>xoffm)
    xoffmax=xoffm;
end
xoff=0:dxinc:xoffmax;
nbins=length(xoff);

%interpolate mute to get a value for each bin
if(isnan(mute))
    mute=[0 0;xoffmax 0];
end
tmute=interp1(mute(:,1),mute(:,2),xoff);

%allocate space for panel
avopanel=zeros(nt,nbins);
fold=avopanel;
tmax=t(end);
%loop over shots and sum to form panel
fmin=[frange(1) frange(1)/2];
fmax=[frange(2) 10];
ievery=5;
t0=clock;
for k=1:nshots
    xs=xshot(k);
    xr=xrec{k};
    xoff2=flag2*abs(xr-xs);
    [xoff2,isort]=sort(xoff2);
    shot=seis{k}(:,isort);
    xr=xr(isort);
    xcmp=.5*(xr+xs);
    %loop over bins
    for j=1:nbins
        bin1=xoff(j)-dxoff;
        bin2=xoff(j)+dxoff;
        inbin=find(xoff2>=bin1 & xoff2<=bin2);
        if(flag1==1)
            %             if(j==100)
            %                 disp('hey')
            %             end
            for jj=1:length(inbin)
                if(xcmp(inbin(jj))>cmpmin && xcmp(inbin(jj))<cmpmax)
                    indt=near(t,tmute(j),tmax);
                    trace=filtf(shot(:,inbin(jj)),t,fmin,fmax,0);
                    avopanel(indt,j)=avopanel(indt,j)+abs(trace(indt));
                    fold(:,j)=fold(:,j)+ones(nt,1);
                end
            end
        else
            for jj=1:length(inbin)
                indt=near(t,tmute(j),tmax);
                trace=filtf(shot(:,inbin(jj)),t,fmin,fmax,0);
                avopanel(indt,j)=avopanel(indt,j)+trace(indt);
                fold(:,j)=fold(:,j)+ones(nt,1);
            end
            
        end
    end
    %     ioff=xoff2<(xoffmax+.01*dxoff);
    %     joff=floor(xoff2(ioff)/dxoff)+1;
    %     if(flag2==1)
    %         xcmp=.5*(xr(ioff)+xs);
    %     else
    %         xcmp=xr(ioff);
    %     end
    %     for j=1:length(joff)
    %         if(xcmp(j)>=cmpmin && xcmp(k)<=cmpmax)
    %             trace=filtf(shot(:,j),t,fmin,fmax,0);
    %             indt=near(t,tmute(joff(j)),tmax);
    %             if(flag1==1)
    %                 avopanel(indt,joff(j))=avopanel(indt,joff(j))+abs(trace(indt));
    %             else
    %                 avopanel(indt,joff(j))=avopanel(indt,joff(j))+trace(indt);
    %             end
    %             fold(:,joff(j))=fold(:,joff(j))+ones(nt,1);
    %         end
    %     end
    if(rem(k,ievery)==0)
        tnow=clock;
        timeused=etime(tnow,t0);
        timeleft=(timeused/k)*(nshots-k);
        disp(['AVOPANEL completed shot ' int2str(k) ' after ' num2str(timeused/60) ' m, time left ' num2str(timeleft/60) ' m'])
    end
end
ind= fold==0;
fold(ind)=1;
avopanel=avopanel./fold;
