function seiszf=wavenumber_highpass(seisz,x,y,sigmax,sigmay)
% WAVENUMBER_HIGHPASS: apply a 2D high-pass wavenumber filter to a depth section or slice.
%
% seiszf=wavenumber_highpass(seisz,x,y,sigmax,sigmay)
%
% It is often desirable to apply a filter that suppresses the lower wavenumbers on a seismic time or
% depth slice or on a depth section.  In the latter case, this might be useful after RTM.  This
% function does that by performing a 2D FFT on the slice and then pointwise multiplying the
% wavenumber spectrum by a Gaussian mask. The Gaussian mask has the same dimensions as the spectrum
% and is just 1-g(kx,ky) where g(kx,ky) is a 2D Gaussian centered at the origin with the specified
% standard deviations. In the case of a 2D depth section, interprete ky as the vertical wavenumber.
%
% seisz ... input time slice or detph section. If a time slice, the column coordinate is x and the
%       row coordinate is y. If a depth section, the column coordinate is x and the depth (row)
%       coordinate is y.
% x ... coordinate for the columns of seisz. Can be either a scalar or a vector. If a scalar it is
%       assumed to be the x coordinate spacing. If a vector, the length(x) must equal size(slice,2).
%       x should be in physical units like feet or meters (i.e. not linenumber) so that wavenumbers
%       are calculated correctly. Must be regularly spaced for the FFT.
% y ... similar to x but for the rows of seisz.
% sigmax ... stdev of mask in kx expressed as a fraction of kxnyq. That is the actual stddev is
%               sigmax*kxnyq. Should be a positive number.
% sigmay ... similar for sigmax but in the ky direction
%   ********** default is sigmay=sigmax **************
%
% seiszf ... the wavenumber filtered slice or depth section
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

if(nargin<5)
    sigmay=sigmax;
end

if(length(x)==1)
    dx=x;
else
    if(size(seisz,2)~=length(x))
        error('x coordinate is the wrong size');
    end
    dx=abs(x(2)-x(1));
end
if(length(y)==1)
    dy=y;
else
    if(size(seisz,1)~=length(y))
        error('y coordinate is the wrong size');
    end
    dy=abs(y(2)-y(1));
end

kny=.5/dy;
knx=.5/dx;

%forward transform
[ny,nx]=size(seisz);
ny2=2^nextpow2(2*ny);
nx2=2^nextpow2(2*nx);
Seisz=fftshift(fft2(seisz,ny2,nx2));
dkx=1/(nx2*dx);
kx=-knx:dkx:knx-dkx;
dky=1/(ny2*dy);
ky=(-kny:dky:kny-dkx)';

%define Gaussian mask

kkx2=kx(ones(size(ky)),:).^2;
kky2=ky(:,ones(size(kx))).^2;
sigmay2=(kny*sigmay)^2;
sigmax2=(knx*sigmax)^2;
gy=exp(-kky2/sigmay2);
gx=exp(-kkx2/sigmax2);
g=1-gx.*gy;
%apply the mask
Seiszf=Seisz.*g;

%Inverse transform
tmp=ifft2(fftshift(Seiszf));
seiszf=real(tmp(1:ny,1:nx));