function [vint,vint2,zint]=fbpicks2vint(z,tfb,zint,nave,zref)
%
% z ... vector of depths at which fb picks are known
% tfb ... vector of fb picks (same size as z)
% zint ... vector of depths at which interval velocities are desired. If zint==z(1) is not included,
%   then it will be, Should be in ascending order.
% nave ... size of a moving average filter applied to the fb picks. Specified in samples. Each input
%   pick will be replaced by the averaage -nave:nave picks
% ************* default = 0 (no averaging) ***********
% zref ... depth at which time is to be preserved. If provided, then the interval velocity in the
%   first layer will be whatever is required to make the integrated total time match the fb time at
%   this depth. This means the actual pick at the first layer will be ignored.
% **** default = nan  which means no action and the first layer is determined by its pick ****
% 
% vint ... vector same size as z with interval velocities
% vint2 ... vector same size as zint with interval velocities
% zint ... same as zint on input unless z(1) needed to be inserted
%

z=z(:);
tfb=tfb(:);
zint=zint(:);

if(length(z)~=length(tfb))
    error('z and tfb must be the same length')
end

if(nargin<5)
    zref=nan;
end

if(nargin<4)
    nave=0;
end

if(nave>0)
    tfb=convz(tfb,ones(2*nave+1,1))/(2*nave+1);
end

ind= zint>=z(1) & zint<=z(end);
zint=zint(ind);

if(zint(1)~=z(1))
    zint=[z(1);zint];
end

vint2=zeros(length(zint),1);
vint=zeros(length(z),1);
for k=2:length(zint)
   i1=near(z,zint(k-1));
   i2=near(z,zint(k));
   z1=z(i1);
   t1=tfb(i1);
   z2=z(i2);
   t2=tfb(i2);
   %interval velocity is the average velocity of the interval
   if(~isnan(t2-t1)&&t2>t1)
       vint2(k)=(z2-z1)/(t2-t1);
   else
       vint2(k)=vprev;
   end
   vint(i1:i2)=vint2(k);
   vprev=vint2(k);
end
vint2(1)=vint2(2);
ind=find(vint==0);
if(~isempty(ind))
    vint(ind)=vint(ind(1)-1);
end

%adjust for reference depth
if(~isnan(zref))
    iref=near(z,zref);
    tref=tfb(iref(1));%actual first break time at reference depth
    if(isnan(tref))
       error('no first break time at reference depth'); 
    end
    t1=vint2t(vint,z);
    t1ref=t1(iref(1));%integrated time at reference depth
    dt=tref-t1ref;%if this is positive then vint(2) is too fast
    dz1=zint(2)-zint(1);%first layer extends from zint(1) to zint(2)
    v1=vint(1);%first layer velocity
    if((dz1+v1*dt)>0)
        v1_new=dz1*v1/(dz1+v1*dt);
        ind=near(z,zint(1),zint(2));
        vint(ind)=v1_new;
        vint2(1:2)=v1_new;
    else
        error('unable to adjust to reference depth');
    end
end