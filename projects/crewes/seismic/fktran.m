function [spec,f,kx]=fktran(seis,t,x,ntpad,nxpad,percent,ishift)
% FKTRAN: Forward fk transform
% [spec,f,kx]=fktran(seis,t,x,ntpad,nxpad,percent,ishift)
%
% FKTRAN uses Matlab's built in fft to perform a 2-d f-k transform
% on a real valued (presumably seismic time-space) matrix. Only the
% positive f's are calculated while all kxs are. Matlab's 'fft' is used
% on each column for the t-f transform while 'ifft' is used on each row
% for the x-kx transform. The inverse transform is performed by ifktran.
%
% seis ... input 2-d seismic matrix. One trace per column.
% t ... vector of time coordinates for seis. length(t) must be the
%	same as number of rows in seis.
% x ... vector of space coordinates for seis. length(x) must be the same
% 	as the number of columns in seis.
% ntpad ... pad seis with zero filled rows until it is this size. Entering nan will also get the
%       default.
%	******* default = next power of 2 ******
% NOTE ... a value of 0 for ntpad is taken to mean no pad is desired
% nxpad ... pad seis with zero filled columns until it is this size. Entering nan will also get the
%           default.
%   ******* default = next power of 2 ******
% NOTE ... a value of 0 for nxpad is taken to mean no pad is desired
% percent ... apply a raised cosine taper to both t and x prior to zero pad.
%	length of taper is theis percentage of the length of the x and t axes
%   ******* default = 0 *********
% ishift ... if 1, then the k axis of the transform is unwrapped to put
%	kx=0 in the middle.
%   ******* default = 1 *******
% spec ... complex valued f-k transform of seis
% f ... vector of frequency coordinates for the rows of spec
% kx ... vector of wavenumber coordinates for the columns of spec
% 
% G.F. Margrave, CREWES Project, U of Calgary, 1996
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

[nsamp,ntr]=size(seis);

if(length(t)~=nsamp)
	error(' Time coordinate vector is incorrect');
end
if(length(x)~=ntr)
	error(' Space coordinate vector is incorrect');
end

if(nargin<7); ishift=1; end
if(nargin<6); percent=0.; end
if(nargin<5); nxpad=nan; end
if(nargin<4); ntpad=2^nextpow2(length(t)); end

if(isnan(nxpad)); nxpad=2^nextpow2(length(x)); end
if(isnan(ntpad)); ntpad=2^nextpow2(length(t)); end

%use fftrl for the t-f transform
[specfx,f]=fftrl(seis,t,percent,ntpad);
clear seis;

% ok taper and pad in x
if(percent>0)
	mw=mwindow(ntr,percent)';
	mw=mw(ones(size(f)),:);
	specfx=specfx.*mw;
	clear mw;
end
if(ntr<nxpad)
	ntr=nxpad;%this causes ifft to apply the x pad
end

%fft on rows
spec = ifft(specfx.',ntr).';

% compute kx
if(x(2)==x(1))
    error('spatial sample rate appears to be zero')
end
kxnyq = 1/(2.*(x(2)-x(1)));
dkx = 2.*kxnyq/ntr;
kx=[0:dkx:kxnyq-dkx -kxnyq:dkx:-dkx];

if(ishift==1)
	[kx,ikx]=sort(kx);
	spec=spec(:,ikx);
end	