function [trout,specinv,smoo]= deconf(trin,trdsign,n,stab,phase,varargin)
% DECONF: frequency domain spiking deconvolution
%
% [trout,specinv]= deconf(trin,trdsign,n,stab,phase,'smoothertype',gaussian,'staboption','mean')
% [trout,specinv]= deconf(trin,trdsign,n,stab,phase,'smoothertype','gaussian')
% [trout,specinv]=deconf(trin,trdsign,n,stab,phase)
% [trout,specinv]=deconf(trin,trdsign,n,stab)
% [trout,specinv]=deconf(trin,trdsign,n)
%
% DECONF performs a frequency domain deconvolution of the
% input trace.  Running the function without any output arguments produces
% a diagnostic plot useful in choosing stab.
%
% trin= input trace to be deconvolved
% trdsign= input trace to be used for operator design
% n= number of points in frequency domain boxcar smoother
% stab= stabilization factor expressed as a fraction of the
%       zero lag of the autocorrelation. This is equivalent to being
%		a fraction of the mean power.
%      ********* default= .0001 **********
% phase= 0 ... zero phase whitening is performed
%        1 ... minimum phase deconvolution is performed
% ************** default= 1 ***************
% Extra arguments. Can be in any order but must be namve, value pairs. That is
% 'smoothertype','gaussian' or 'staboption','mean' or both.
% 'smoothertype' ... smoother. This is a pair of inputs prescribing the
%        type of spectral smoother. 'smoothertype','boxcar' is the default
%        and prescrives boxcar (convolutional) smoothing. Other
%        possibilities are 'smoothertype','triangle' and 'smoothertype','gaussian'
%        which have obvious interpretations.
% 'staboption' ... determines how the stab constant is used. Two possible choices either 'max'
%       or 'mean'. 'max' is the default because this allows stab*pmax to be easily interpreted
%       as a decibel level. However, such stab values are then not equivalent to those used in
%       deconw. For direct equivalence, use 'mean' and then the stab power is stab*pmean.
%
% trout= output trace which is the deconvolution of trin
% specinv= output inverse operator spectrum. The time domain operator can be
%          recovered by real(ifft(fftshift(specinv)))
%
% by: G.F. Margrave, May 1991
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

% set defaults
nargs=nargin-5;
if(nargin==6)
    nargs=length(varargin);
    if(2*floor(nargs/2)~=nargs)
        error('extra arguments must be name-value pairs')
    end
end
if nargin < 5
    phase=1;
end
if nargin<4
    stab=.0001;
end
smoother='boxcar';
staboption='max';
if(nargs>0)
    for k=1:2:nargs
        recognized=0;
        if(strcmp(varargin{k},'smoothertype'))
            smoother=varargin{k+1};
            
            if(strcmp(smoother,'boxcar'))
                recognized=1;
            elseif(strcmp(smoother,'triangle'))
                recognized=1;
            elseif(strcmp(smoother,'gaussian'))
                recognized=1;
            else
                error('unknown smoother type');
            end
        elseif(strcmp(varargin{k},'staboption'))
            staboption=varargin{k+1};
            if(strcmp(staboption,'max'))
                recognized=1;
            elseif(strcmp(staboption,'mean'))
                recognized=1;
            else
                error('unknown stab option');
            end
                
        end
        if(~recognized)
            error('unrecognized argument name')
        end
    end
end

% switch to column vectors
[k,j]=size(trin);
if((k-1)*(j-1)~=0)
    error('deconf takes one trace at a time only')
end

if(k==1); trin=trin.'; trdsign=trdsign.'; end

% pad input traces to power of 2
N=length(trin);
trin=padpow2(trin);
trdsign=pad_trace(trdsign,trin);
Npad=length(trin);
%increase n to account for padding
n = Npad*n/N;
% generate the power spectrum (pad to length of trin)
%m=length(trin)-length(trdsign);
spec= fftshift(fft(trdsign));
power= real(spec).^2 + imag(spec).^2;
% stabilize the power spectrum
if(strcmp(staboption,'max'))
    P0=max(power);
else
    P0=sum(power)/length(power);
end
delta_p= stab*P0;
if(nargout==0)
    figure
    a=sqrt(power(length(power)/2:end));
    ma=max(a);
    %fake frequency axis
    nf=length(a);
    f=linspace(0,pi,nf);
    s1=sqrt(stab)*ma*ones(size(a));
    s2=10*s1;
    if(max(s2)>ma); s2=.9*ma; end
    s3=s1/10;
    plot(f,todb(a,ma),f,todb(s1,ma),f,todb(s2,ma),f,todb(s3,ma));
    legend('A(f) of trace','stab*max(A)','10*stab*max(A)','stab*max(A)/10')
    xlabel('frequency (pi=Nyquist)')
    return
end
power=power+delta_p;
% create frequency smoother (odd number of points)
nn=2*fix(n/2)+1;
if(strcmp(smoother,'boxcar'))
    smoo=ones(nn,1);%want an odd number of samples
elseif(strcmp(smoother,'triangle'))
    smoo=triangle(1:2*nn-1,nn,2*nn);
elseif(strcmp(smoother,'gaussian'))
    smoo=gaussian(1:2*nn-1,nn,2*nn);
end
%power=convz(power,smoo)/sum(abs(smoo));
power=specsmth(power,smoo);
% symmetrize the power spectrum
n2=length(power);
power(n2/2+2:end)=power(n2/2:-1:2);
% compute the minimum phase spectrum
if phase==1
    %logspec=hilbert(.5*log(power));
    logspec=herbert(.5*log(power));
    % compute the complex spectrum of the inverse operator
    specinv= exp(-conj(logspec));
else
    specinv=power.^(-.5);
end
% deconvolve the input trace
specin=fftshift(fft(trin));
specout=specin.*specinv;
trout=real(ifft(fftshift(specout)));
%   trout=balans(trout,trin);

%unpad to length of trin
trout = pad_trace(trout,1:N);

if( k==1 ); trout =trout.'; end

end

function pwr=specsmth(pwr,smo)
    %extend pwr on both ends
    np=length(pwr);
    ns=length(smo);
%     pwr2=zeros(np+ns,1);
%     ipw=floor(ns/2)+1:floor(ns/2)+np;
    pwr2=zeros(np+2*ns,1);
    ipw=ns+1:ns+np;
    pwr2(ipw)=pwr;
    pwr2(1:ipw(1)-1)=pwr(1);
    pwr2(ipw(end)+1:end)=pwr(end);
    tmp=convz(pwr2,smo,ceil(ns/2),length(pwr2),0)/sum(abs(smo));
    pwr=tmp(ipw);
end