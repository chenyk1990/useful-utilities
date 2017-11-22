function [vel,x,z]=twodipmodel(dx,xmax,zmax,vone,vtwo,vthree,dips,znots)
% dipmodel: build a simple 2D model of a dipping reflector
% 
% [vel,x,z]=twodipmodel(dx,xmax,zmax,vone,vtwo,vthree,dips,znots)
%
% dx ... grid interval (distance between grid points in x and z)
% xmax ... maximum x coordinate (minimum is zero)
%  *********** default 2500 **********
% zmax ... maximum z coordinate (minimum is zero)
%  *********** default 1000 ************
% vlow ... velocity above the first dipping reflector
%  *********** default 2000 ************
% vtwo ... velocity below the first dipping reflector and above the second
%  *********** default 3000 ************
% vthree ... velocity below the second reflector
% ************ default =4000 ***********
% dips ... dips in degrees of both reflectors, positive is down to the right
%  *********** default = [10 -10] degrees ************
% znots ... depth to each dipping reflector at xmax/2
%  *********** default = [zmax/4 zmax/2] **************
%
% vel ... velocity model matrix
% x ... x coordinate vector for vel
% z ... z coordinate vector for vel
%
% NOTE: the simplest way to plot vel is: plotimage(vel-mean(vel(:)),z,x)
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

%

if(nargin==0)
    [vel,x,z]=twodipmodel(10);
    figure
    imagesc(x,z,vel);colorbar
    title('Two Dip model, colors indicate velocity')
    xlabel('distance (m)')
    zlabel('depth (m)')
    clear
    return;
end
% [vel,x,z]=twodipmodel(dx,xmax,zmax,vone,vtwo,vthree,dips,znots)
if(nargin<4)
    vone=2000;
end
if(nargin<5)
    vtwo=3000;
end
if(nargin<6)
    vthree=4000;
end
if(nargin<3)
    zmax=1000;
end
if(nargin<2)
    xmax=2500;
end
if(nargin<7)
    dips=[5 -5];
end
if(nargin<8)
    znots=[zmax/4 zmax/2];
end

x=0:dx:xmax; % x coordinate vector
z=0:dx:zmax; % z coordinate vector

%initialize velocity matrix as a constant matrix full of vlow
vel=vone*ones(length(z),length(x));

% define the first dipping horizon as a four point polygon
dx2=dx/2;
xnot=xmax/2;
x1=x(1)-dx2;x2=x(end)+dx2;
z1=znots(1)+(x1-xnot)*tand(dips(1));
z2=znots(1)+(x2-xnot)*tand(dips(1));

xpoly=[x1 xnot  x2 x2 x1];zpoly=[z1 znots(1) z2 zmax+dx2 zmax+dx2];

% install the basement in the velocity matrix
vel=afd_vmodel(dx,vel,vtwo,xpoly,zpoly);

% define the second dipping horizon as a four point polygon
dx2=dx/2;
xnot=xmax/2;
x1=x(1)-dx2;x2=x(end)+dx2;
z1=znots(2)+(x1-xnot)*tand(dips(2));
z2=znots(2)+(x2-xnot)*tand(dips(2));

xpoly=[x1 xnot  x2 x2 x1];zpoly=[z1 znots(2) z2 zmax+dx2 zmax+dx2];

% install the basement in the velocity matrix
vel=afd_vmodel(dx,vel,vthree,xpoly,zpoly);