function [g,x,y]=gaus_radial(dx,sigmax,dy,sigmay,nhw)
% GAUS_RADIAL: makes a 2D radial Gaussian suitable for convolutional smoothing
% 
% [g,x,y]=gaus_radial(dx,sigmax,dy,sigmay,nhw)
% 
% dx ... grid size in x (physical units) (distance between columns)
% sigmax ... Gaussian half-width (standard deviation) in x 
% dy ... grid size in y (physical units) (distance between rows)
% ********** default dy=dx **********
% sigmay ... Gaussian half-width (standard deviation) in y
% ********** default sigmy=sigmax ***********
% nhw ... number of half widths that the gaussian spans. The resulting radial gaussian will
% extend from -nhw*sigmax to +nhw*sigmx in x and similarly in y. The number of samples in both
% dimensions will be odd.
% ********** default nhw=3 **********
% 
% g ... the 2D radial Gaussian
% x ... x coordinate for g
% y ... y coordinate for g
%
% Example:
% [g,x,y]=gaus_radial(1,50,1,50,3);
% figure;surf(x,y,g);shading flat
%
%
% G.F. Margrave, Devon Energy, 2017
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

if(nargin<3)
    dy=dx;
end
if(nargin<4)
    sigmay=sigmax;
end
if(nargin<5)
    nhw=3;
end
if(sigmax<=dx)
    error('sigmax should be greater than dx')
end
if(sigmay<=dy)
    error('sigmay should be greater than dy')
end
if(nhw<1)
    error('nhw should be 1 or greater');
end

x=-sigmax*nhw:dx:sigmax*nhw;
y=(-sigmay*nhw:dy:sigmay*nhw)';

ixnot=near(x,0);
xnot=x(ixnot(1));
iynot=near(y,0);
ynot=y(iynot(1));

xx=x(ones(size(y)),:);
yy=y(:,ones(size(x)));

g=exp(-(xx-xnot).^2/sigmax^2-(yy-ynot).^2/sigmay^2);
