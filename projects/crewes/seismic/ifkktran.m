function [seis,t,x,y]=ifkktran(spec,f,kx,ky,nt,nx,ny)
% IFKKTRAN inverse f-kx-ky transform (3D)
% [seis,t,x,y]=ifkktran(spec,f,kx,ky,nt,nx,ny)
%
% IFKKTRAN uses Matlab's built in fft to perform a 3D inverse f-kx-ky
% transform on a complex valued (presumably seismic f-kx-ky) matrix.  It is
% assumed that the forward transform was performed by fkktran which means
% that only non-negative f's will be present.
%
% spec ... complex valued f-kx-ky transform as a 3D matrix
% f ... vector of frequency coordinates for the first dimension of spec.  
%	length(f) must be the same as size(spec,1).
% kx ... vector of wavenumber coordinates for the second dimension of spec
%	length(kx) must equal size(spec,2).
% ky ... vector of wavenumber coordinates for the third dimension of spec
%	length(ky) must equal size(spec,3).
% nt ... size of the output seismic matrix in the first dimension. This can
%   be specified to remove any padding that was applied in the forward
%   transform. Must be no larger than 2*length(nf)-1
% ***************** default =2*length(f)-1 (no padding removed) ***********
% nx ... size of the output seismic matrix in the second dimension. This
%   can be specified to remove any padding that was applied in the forward
%   transform. Must be no larger than length(kx).
% ***************** default = length(kx) (no padding removed) ************
% ny ... size of the output seismic matrix in the third dimension. This
%   can be specified to remove any padding that was applied in the forward
%   transform. Must be no larger than length(ky).
% ***************** default = length(ky) (no padding removed) ************
% seis ... output 3D seismic matrix. Time is dimension 1, x is dimension 2
%   and y is dimension 3.
% t ... vector of time coordinates for seis. 
% x ... vector of first space coordinates for seis. 
% y ... vector of second space coordinates for seis
% 
% G.F. Margrave, Devon Energy 2017
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any
% purpose with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this
% software's Matlab source file.

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

[nf,nkx,nky]=size(spec);

if(nargin<5)
    nt=2*nf-1;
end
if(nargin<6)
    nx=nkx;
end
if(nargin<7)
    ny=nky;
end


if(length(f)~=nf)
	error(' Frequency coordinate vector is incorrect');
end
if(length(kx)~=nkx)
	error(' x wavenumber coordinate vector is incorrect');
end
if(length(kx)~=nkx)
	error(' y wavenumber coordinate vector is incorrect');
end


if(nx>nkx)
    error('nx cannot be greater than length(nkx)');
end
if(ny>nky)
    error('ny cannot be greater than length(nky)');
end

%determine if we need fftshift or not
ishift=0;
if(kx(1)~=0)
    ishift=1;
end


%ok 2D wavenumber transform on frequency samples
%disp('kx-x')
for k=1:nf
    if(ishift)
        tmp = fft2(fftshift(squeeze(spec(k,:,:))));
    else
        tmp = fft2(squeeze(spec(k,:,:)));
    end
    spec(k,:,:)=shiftdim(tmp,-1);
end

if(nx~=nkx || ny~=nky)
   %spec2=zeros(nf,nx,ny);
   spec2=spec(:,1:nx,1:ny);
else
    spec2=spec;
end

clear spec;

%use ifftrl for the f-t transform
%disp('f-t');
[seis,t]=ifftrl(spec2,f);
if(nt<size(seis,1))
    dt=t(2)-t(1);
    nt2=length(t);
    seis(nt+1:nt2,:,:)=[];
    t=dt*(0:nt-1)';
end
clear spec2;

% compute x and y
dx=.5/max(kx);
x=dx*(0:nx-1);

dy=.5/max(ky);
y=dy*(0:ny-1);