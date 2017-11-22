function [stack,xr]=recstack(shots,xrecs,t,xshots,mute)
% RECSTACK: common receiver stack
%
% [stack,xr]=recstack(shots,xrecs,t,xshots,mute)
%
% RECSTACK does a common receiver stack (no moveout correction)
% 
% shots ... cell array of shots to be stacked. Each shot is an (x,t) matrix
% xrecs ... cell array of receiver coordinates for each shot. May be a single vector if all shots
%       have the same size.
% t ... t coordinate vector of each shot and of the resulting stack
% xshots ... vector of x coordinates for each shot (one number per shot
%   giving the shot location)
% mute ... (offset,t) values of a mute to be applied to each shot. Offsets
%   greater than the maximum offset in the mute will be zero'd. Specified like mute=[off1 t1;off2
%   t2;....]. Mute must not be double valued (offset must always increase). Here offset means
%   abs(xrecs{j}(k)-xshot(j)) where k is trace number and j is shot number.
%
% stack ... stack of gathers normalized by live (non-zero) fold
% xr ... receiver coordinate vector for the stack
%
%G.F. Margrave 2017
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


if(length(shots)~=length(xshots))
    error('shots and xshots must be the same length');
end
if(~iscell(xrecs))
    if(size(shots{1},2)==length(xrecs))
        %copy into cell array
        tmp=cell(size(shots));
        for k=1:length(shots)
            tmp{k}=xrecs;
        end
        xrecs=tmp;
        clear tmp;
    end
end
if(length(shots)~=length(xrecs))
    error('shots and xrecs must be the same length');
end

%determine the live shots and total traces
nshots=length(shots);
ilive=1:nshots;
ntraces=0;
for k=1:nshots
    if(isempty(shots{k}))
        ilive(k)=0;
    end
    ntraces=ntraces+size(shots{k},2);
end
ind=find(ilive==0);
if(~isempty(ind))
    ilive(ind)=[];
end
nliveshots=length(ilive);
%check sizes
for k=1:nshots
   if(size(shots{k})~=size(xrecs{k}))
       error(['size of shot record ' int2str(k) ' does not match xrecs']);
   end
end

%determine unique receivers
%the following only works if the receiver spacing is regular
xrmin=inf;
xrmax=-inf;
for k=1:length(xrecs)
    if(min(xrecs{k})<xrmin)
        xrmin=min(xrecs{k});
    end
    if(max(xrecs{k})>xrmax)
        xrmax=max(xrecs{k});
    end
end
dx=xrecs{1}(2)-xrecs{1}(1);
xr=xrmin:dx:xrmax;%this will be the output x coordinate
        

stack=zeros(length(t),length(xr));
fold=stack;

% gathers=cell(size(x));
% offsets=gathers;
%prepopulate gathers so we don't have growing arrays
% gath0=zeros(length(z),nshots);
% off0=zeros(1,nshots);
% for k=1:length(x)
%     gathers{k}=gath0;
%     offsets{k}=off0;
% end


offmute=mute(:,1);
tmute=mute(:,2);
if(offmute(1)~=0)
    offmute=[0;offmute];
    tmute=[0;tmute];
end
offmax=max(offmute);
for j=1:nliveshots
    shot=shots{ilive(j)};
    x=xrecs{ilive(j)};
    for k=1:length(x);
        off=abs(x(k)-xshots(ilive(j)));
        if(off<=offmax)
            trace=shot(:,k);
            count=ones(size(trace));
            tbegin=interp1(offmute,tmute,off);
            ind=find(t<tbegin);
            if(~isempty(ind))
                trace(ind)=0;
                count(ind)=0;
            end
            istk=near(xr,x(k));
            stack(:,istk)=stack(:,istk)+trace;
            fold(:,istk)=fold(:,istk)+count;
        end
    end
end

%normalize
ind=find(fold>0);
stack(ind)=stack(ind)./fold(ind);

%taper if asked
% if(taper>0)
%     mask=ones(size(stack));
%     ind= fold==0;
%     mask(ind)=0;
%     mask=gaussian_smoother(mask,x,z,taper);
%     stack=stack.*mask;
% end

%process the gathers. Delete zero traces, sort into increasing offset
% for k=1:length(x)
%     gathtmp=gathers{k};
%     offtmp=offsets{k};
%     test=sum(abs(gathtmp));
%     ind= test==0;
%     gathtmp(:,ind)=[];
%     offtmp(ind)=[];
%     
%     [tmp,isort]=sort(offtmp);
%     gathers{k}=gathtmp(:,isort);
%     offsets{k}=tmp;
% end
    
        
    