function [seis,t,x,v]=diffraction_section(dx,dt,xmax,zmax,v0,c,fmin,fmax,ieveryt,ieveryx)
%DIFFRACTION_SECTION: make a section full of diffractions on a regular grid
%
% [seis,t,x,seisz,z]=diffraction_section(dx,dt,xmax,tmax,v0,c,fmin,fmax,ieveryz,ieveryx)
%
% dx ... spatial grid size (in both x and z)
% ********** default=10 ************
% dt ... time sample size
% ********** default =.002 ***********
% xmax ... maximum lateral spatial coordinate
% ********** default = 2000 *************
% zmax ... maximum depth 
% ********** default = 2000 ************
% v0 ... velocity at z=0
% ********* default = 1800 **********
% c ... vertical velocity gradient
% ********* default = 0.6 ***********
% Note: velocity function is v=v0+c*z. Set c=0 for constant velocity
% fmin ... low frequency of pass band
% ********** default = 10 ***********
% fmax ... high frequency of passband
% ********** default = .4*fnyq ******** 
% ieveryt ... place a diffractor every this many time samples
% ********** default = 50 *********
% ieveryx .... place a diffractor every this many x samples
% ********** default = 20 ********
% 
% seis ... modelled diffraction section
% t ... time coordinate for seis
% x ... x coordinate for seis
% v ... instantaeous velocity as a function of time. Same size as t
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

if(nargin<1)
    dx=10;
end
if(nargin<2)
    dt=.002;
end
fnyq=.5/dt;
if(nargin<3)
    xmax=2000;
end
if(nargin<4)
    zmax=2000;
end
if(nargin<5)
    v0=1800;
end
if(nargin<6)
    c=.6;
end
if(nargin<7)
    fmin=10;
end
if(nargin<8)
   fmax=.4*fnyq;
end
if(nargin<9)
    ieveryt=50;
end
if(nargin<10)
    ieveryx=20;
end
x=0:dx:xmax;
nx=length(x);
%linear v(z) function
z=0:dx:zmax;
vz=v0+c*z;
tv=2*vint2t(vz,z);%2-way vertical traveltime
%vave=z./(.5*tv);%average velocity versus z or tv
tmax=dt*(ceil(tv(end)/dt)+1);
t=(0:dt:tmax)';
nt=length(t);
v=interp1(tv,vz,t);%instantaneous v(t)
ind=find(isnan(v));%look for nans at the end
if(~isempty(ind))
    v(ind)=v(ind(1)-1);
end

fmin2=[fmin .25*fmin];fmax2=[fmax min([20 .25*(fnyq-fmax)])];
seis=zeros(nt,nx);
%determine diffraction locations
itdiff=1:ieveryt:nt;
if(itdiff(end)~=nt)
    itdiff=[itdiff,nt];
end
ixdiff=1:ieveryx:nx;    
tmp=zeros(size(t));
tmp(itdiff)=1;
%bandlimit
tmp=filtf(tmp,t,fmin2,fmax2);
%install
for k=ixdiff
    seis(:,k)=tmp;
end

%set up for modelling
params=nan*(1:13);
params(1)=min([fmax2(1)+fmax2(2), fnyq]);
params(2)=0;
params(13)=50;
[seis,t,x]=vz_fkmod(seis,v,t,x,params);
