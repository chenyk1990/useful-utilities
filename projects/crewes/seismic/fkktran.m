function [spec,f,kx,ky]=fkktran(seis,t,x,y,ntpad,nxpad,nypad,percent,ishift)
% FKKTRAN: Forward f-kx-ky transform (3D)
% [spec,f,kx,ky]=fkktran(seis,t,x,y,ntpad,nxpad,nypad,percent,ishift)
%
% FKKTRAN uses Matlab's built in fft to perform a 3-D f-kx-ky transform
% on a real valued (presumably seismic time-space) matrix. Only the
% positive f's are calculated while all kxs and kys are. Matlab's 'fft' is used
% on each column for the t-f transform while 'ifft' is used on each time slice
% for the x-kx,y-ky  transform. The inverse transform is performed by ifkktran.
%
% seis ... input 2-d seismic matrix. One trace per column.
% t ... vector of time coordinates for seis. length(t) must be the
%	same as number of rows in seis.
% x ... vector of the first space coordinates for seis. length(x) must be the same
% 	as the size(seis,2).
% y ... vector of the second space coordinates for seis. length(y) must be the same
% 	as the size(seis,3).
% ntpad ... pad seis with zero filled rows until it is this size.
%	******* default = next power of 2 ******
% NOTE ... a value of 0 for ntpad is taken to mean no pad is desired
% nxpad ... pad seis with zeros in the second dimension until it is this size.
%   ******* default = next power of 2 ******
% NOTE ... a value of 0 for nxpad is taken to mean no pad is desired
% nypad ... pad seis with zeros in the third dimension until it is this size.
%   ******* default = next power of 2 ******
% NOTE ... a value of 0 for nypad is taken to mean no pad is desired
% percent ... apply a raised cosine taper to t, x, and y prior to zero pad.
%	length of taper is theis percentage of the length of the x and t axes
%   ******* default = 0 *********
% ishift ... if 1, then the k axis of the transform is unwrapped to put
%	kx=0 in the middle.
%   ******* default = 1 *******
% spec ... complex valued f-kx-ky transform of seis (3D matrix)
% f ... vector of frequency coordinates for the first dimension of spec
% kx ... vector of wavenumber coordinates for the second dimension of spec
% ky ... vector of wavenumber coordinates for the third dimension of spec
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

[nt,nx,ny]=size(seis);

if(ny==1)
    error('Data matrix is not 3D');
end

if(length(t)~=nt)
	error(' t coordinate vector is incorrect');
end
if(length(x)~=nx)
	error(' x coordinate vector is incorrect');
end
dx=abs(x(2)-x(1));
if(sum(abs(diff(diff(x))))>.00001*dx)
    error(' x coordinate must be regularly sampled');
end
if(dx==0)
    error('x sample size cannot be zero');
end
if(length(y)~=ny)
	error(' y coordinate vector is incorrect');
end
dy=abs(y(2)-y(1));
if(sum(abs(diff(diff(y))))>.00001*dy)
    error(' y coordinate must be regularly sampled');
end
if(dy==0)
    error('y sample size cannot be zero');
end

if(nargin<9); ishift=1; end
if(nargin<8); percent=0.; end
if(nargin<7); nypad=2^nextpow2(length(y)); end
if(nargin<6); nxpad=2^nextpow2(length(x)); end
if(nargin<5); ntpad=2^nextpow2(length(t)); end

%use fftrl for the t-f transform
[specfxy,f]=fftrl(seis,t,percent,ntpad);
nf=length(f);
clear seis;

% ok taper and pad in x and y
if(percent>0)
	mwx=mwindow(nx,percent);%column vector
    mwy=mwindow(ny,percent)';%row vector
	mw=shiftdim(mwy(ones(1,nx),:).*mwx(:,ones(ny,1)),-1);
    for k=1:nt
        specfxy(k,:,:)=mw.*specfxy(k,:,:);
    end
    clear mw;
end

if(nx<nxpad)
	nx=nxpad;%this causes ifft2 to apply the x pad
end
if(ny<nypad)
	ny=nypad;%this causes ifft2 to apply the y pad
end

%ifft on frequency samples
for k=1:nf
    tmp = ifft2(squeeze(specfxy(k,:,:)),nx,ny);
    if(k==1)
        if(isa(tmp,'single'))
            spec=single(zeros(nf,size(tmp,1),size(tmp,2)));
        else
            spec=zeros(nf,size(tmp,1),size(tmp,2));
        end
    end
    if(ishift==1)
        spec(k,:,:)=fftshift(shiftdim(tmp,-1));
    else
        spec(k,:,:)=shiftdim(tmp,-1);
    end
end

% compute kx and ky
kxnyq = 1/(2.*dx);
dkx = 2.*kxnyq/nx;
kx=[0:dkx:kxnyq-dkx -kxnyq:dkx:-dkx];

kynyq = 1/(2.*dy);
dky = 2.*kynyq/ny;
ky=[0:dky:kynyq-dky -kynyq:dky:-dky];

if(ishift==1)
	kx=sort(kx);
	ky=sort(ky);
end	