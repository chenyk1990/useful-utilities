function slicem=wavenumber_gaussmask(slice,x,y,sigmax,sigmay)
% WAVENUMBER_GAUSSMASK: apply a 2D Gaussian mask to the 2DFFT of a seismic timeslice or depthslice
%
% slicem=wavenumber_gaussmask(slice,x,y,sigmax,sigmay)
%
% It is often desirable to apply a filter that suppresses the higher wavenumbers on a seismic time
% or depth slice. This function does that by performing a 2D FFT on the slice and then pointwise
% multiplying the wavenumber spectrum by a Gaussian mask. The Gaussian mask has the same dimensions
% as the spectrum and is a either a decaying or growing Guassian in both kx and ky. Both masks are
% unity at the origin and decay (or grow) away from the origin. The total mask is the poitwise
% product of the kx and ky masks.
%
% slice ... input time slice, the column coordinate is x and the row coordinate is y.
% x ... coordinate for the columns of slice. The length(x) must equal size(slice,2). x should be
%       in physical units like feet or meters (i.e. not linenumber) so that wavenumbers are
%       calculated correctly. Must be regularly spaced for the FFT.
% y ... coordinate for the rows of slice. The length(y) must equal size(slice,1). y should be
%       in physical units like feet or meters (i.e. not linenumber) so that wavenumbers are
%       calculated correctly. Must be regularly spaced for the FFT.
% sigmax ... stdev of mask in kx expressed as a fraction of kxnyq. That is the actual stddev is
%               sigmax*kxnyq. Positive values give decay while negative values mean growth.
% sigmay ... similar for sigmax but in the ky direction
%
% slicem ... the wavenumber filtered slice
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

%forward transform
[ny,nx]=size(slice);
ny2=2^nextpow2(2*ny);
nx2=2^nextpow2(2*nx);
[Slice,ky,kx]=fktran(slice,y,x,ny2,nx2);

%define Gaussian mask
kny=max(ky);
knx=max(kx);
kkx2=kx(ones(size(ky)),:).^2;
kky2=ky(:,ones(size(kx))).^2;
sigmay2=(kny*sigmay)^2;
sigmax2=(knx*sigmax)^2;
gy=exp(-kky2/sigmay2);
gx=exp(-kkx2/sigmax2);

%apply the mask
Slicem=Slice.*gx.*gy;

%Inverse transform
tmp=ifktran(Slicem,ky,kx);
slicem=tmp(1:ny,1:nx);