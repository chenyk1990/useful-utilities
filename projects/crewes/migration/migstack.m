function [stack,gathers,offsets]=migstack(shots,x,z,xshots,mute,killrad,taper,gainopt)
% MIGSTACK: stack a 2D line of migrated shots whill applying a mute and gain
%
% stack=migstack(shots,x,z,xshots,mute,killrad,taper,gainopt)
%
% MIGSTACK simply stacks a line of migrated shots while applying a mute (to
% the shots not the CIG's). THe option exists to gain the shots as may be
% appropriate if a cross-correlation imaging condition was used. This
% function assumes the shots were all padded with zeros traces to the exact
% dimensions of the velocity model and then put through psdm (see
% pspi_shot).
% 
% shots ... cell array of shots to be stacked. Each shot is an (x,z) matrix
%   of identical size, however, null (completely empty) shots are allowed.
% x ... x coordinate vector of each shot and of the resulting stack
% z ... z coordinate vector of each shot and of the resulting stack
%NOTE: all shots must have the same x and z coordinates)
% xshots ... vector of x coordinates for each shot (one number per shot
% giving the shot location)
% mute ... (offset,z) values of a mute to be applied to each shot. Offsets
%   greater than the maximum offset in the mute will be zero'd. Specified
%   like mute=[off1 z1;off2 z2;....]. Mute must not be double valued. Here offset means
%   abs(x(k)-xshot(j)) where k is trace number and j is shot number. Thus offset horizontal
%   distance from the shotpoint to the image point.
% killrad ... within this radial distance of the shot, amplitudes will be
%   zero'd
% taper ... width of a taper around the nonzero part of the stack
% gainopt ... 1 means apply gain (appropriate if cc imaging condition)
%             0 means don't apply
%         see GAINCC to understand how gain is applied.
% ************* default = 0 *************
% stack ... stack of gathers normalized by live (non-zero) fold
% gathers ... cell array of common image gathers, one per each x coordinate
% offsets ... cell array of offsets for each gather. If a gather has 100 traces then there will
%   be 100 offsets in this array for that gather.
%
%G.F. Margrave 2011
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


if(nargin<8)
    gainopt=0;
end
if(length(shots)~=length(xshots))
    error('shots and xshots must be the same length');
end

%determine the live shots
nshots=length(shots);
ilive=1:nshots;
for k=1:nshots
    if(isempty(shots{k}))
        ilive(k)=0;
    end
end
ind=find(ilive==0);
if(~isempty(ind))
    ilive(ind)=[];
end
nshots=length(ilive);
        

stack=zeros(size(shots{ilive(1)}));
fold=stack;

gathers=cell(size(x));
offsets=gathers;
%prepopulate gathers so we don't have growing arrays
gath0=zeros(length(z),nshots);
off0=zeros(1,nshots);
for k=1:length(x)
    gathers{k}=gath0;
    offsets{k}=off0;
end


offmute=mute(:,1);
zmute=mute(:,2);
if(offmute(1)~=0)
    offmute=[0;offmute];
    zmute=[0;zmute];
end
offmax=max(offmute);
for j=1:nshots
    shot=shots{ilive(j)};
    if(gainopt)
        shot=gaincc(shot,x,z,xshots(ilive(j)),0);
    end
    for k=1:length(x);
        off=abs(x(k)-xshots(ilive(j)));
        r=sqrt(z.^2+off^2);
        if(off<=offmax)
            trace=shot(:,k);
            count=ones(size(trace));
            zbegin=interp1(offmute,zmute,off);
            ind=find(z<zbegin);
            if(~isempty(ind))
                trace(ind)=0;
                count(ind)=0;
            end
            ind=find(r<killrad);
            if(~isempty(ind))
                trace(ind)=0;
                count(ind)=0;
            end
            stack(:,k)=stack(:,k)+trace;
            fold(:,k)=fold(:,k)+count;
            gathers{k}(:,j)=trace;
            offsets{k}(j)=off;
        end
    end
end

%normalize
ind=find(fold>0);
stack(ind)=stack(ind)./fold(ind);

%taper if asked
if(taper>0)
    mask=ones(size(stack));
    ind= fold==0;
    mask(ind)=0;
    mask=gaussian_smoother(mask,x,z,taper);
    stack=stack.*mask;
end

%process the gathers. Delete zero traces, sort into increasing offset
for k=1:length(x)
    gathtmp=gathers{k};
    offtmp=offsets{k};
    test=sum(abs(gathtmp));
    ind= test==0;
    gathtmp(:,ind)=[];
    offtmp(ind)=[];
    
    [tmp,isort]=sort(offtmp);
    gathers{k}=gathtmp(:,isort);
    offsets{k}=tmp;
end
    
        
    