function [seis,specinfo]=tvgaussianmask(seis,t,x,y,tspec,sigmax,sigmay)
%
% [seis,specinfo]=tvgaussianmask(seis,t,x,y,tspec,sigmax,sigmay)
%
% Apply a time-variant Gaussian mask in the wavenumber domain to a 3D dataset. This is equivalent to
% a time-variant spatial convolution on time slices. The normal gaussian is a decaying exponential
% but a growing one, a spatial deconvolution, is also accomodated.
%
% seis ... input 3D dataset. First dimension is time, second ix x (xline) and third is y (iline).
%       Must be resularly sampled in all three dimensions. 
% t ... time coordinate vector, length(t) must equal size(seis,1).
% x ... x coordinate vector, length(x) must equal size(seis,2).
% y ... y coordinate vector, length(y) must equal size(seis,3).
% tspec ... vector of times at which the spatial filter sizes are specified.
% sigmax ... vector of standard deviations of the Gaussian mask in the x direction expressed as a 
%       fraction of x spatial Nyquist, eg a value of .25 means .25*kxnyq. length(tspec) must equal
%       length(sigmax). Negative values for sigmax signify that the mask is division by a Gaissian
%       which causes spatial deconvolution (sharpening) while positive values cause spatial
%       convolution (blurring). Generally, sigmax should be a number between -1 and +1. However,
%       values greater than 1 are allowed and the larger it is the less action. Values of 10 or
%       greater are taken as a flag for no action and such time slices are simply passed to output.
% sigmay ... Just like sigmax except for the y direction. x and y filtering actions can be completley 
%       independent.
% *********** default sigmay=sigmay **************
% 
% seis (outout) ... 3D seismic matrix the same size as the input with the spatial filtering applied
% specinfo .... cell array of length 6 with information about the input and output wavenumber
%       spectra. The 6 entries are
%       specinfo{1} ... vector of length(tspec) giving the average kx spectrum at each time tspec
%               for the input data. Each spectrum is a vector to be plotted versus kx.
%       specinfo{2} ... vector of length(tspec) giving the average ky spectrum at each time tspec
%               for the input data. Each spectrum is a vector to be plotted versus ky.
%       specinfo{3} ... like specinfo{1} except for the output data
%       specinfo{4} ... like specinfo{2} except for the output data
%       specinfo{5} ... vector of kx values for plotting specinfo 1 and 3
%       specinfo{6} ... vector of ky values for plotting specinfo 2 and 4
%
%
% NOTE: Because we are effectively convolving each time slice with a different size spatial filter,
% there is a possible time-variant amplitude effect. At present, this is addressed by balancing
% input to output rms power on each time slice.
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

if(nargin<7)
    sigmay=sigmax;
end
disp('tvgaussianmask')
[nt,nx,ny]=size(seis);
if(length(t)~=nt)
    error('t coordinate is incompatible with seis')
end
if(length(x)~=nx)
    error('x coordinate is incompatible with seis')
end
if(length(y)~=ny)
    error('y coordinate is incompatible with seis')
end
if(length(tspec)~=length(sigmax))
    error('tspec and sigmax must have the same length')
end
if(length(tspec)~=length(sigmay))
    error('tspec and sigmay must have the same length')
end

tspec=tspec(:)';
sigmax=sigmax(:)';
sigmay=sigmay(:)';

%expand filter specification to cover the entire time range
ind=find(tspec==t(1), 1);
if(isempty(ind))
    tspec=[t(1) tspec];
    sigmax=[sigmax(1) sigmax];
    sigmay=[sigmay(1) sigmay];
end
ind=find(tspec==t(end), 1);
if(isempty(ind))
    tspec=[tspec t(end)];
    sigmax=[sigmax sigmax(end)];
    sigmay=[sigmay sigmay(end)];
end

sigmaxt=interp1(tspec,sigmax,t);
sigmayt=interp1(tspec,sigmay,t);
dt=t(2)-t(1);

if(ny>1)
    %3D case
    nx2=2^nextpow2(nx);
    ny2=2^nextpow2(ny);
    knx=.5/abs(x(2)-x(1));
    kny=.5/abs(y(2)-y(1));
    ievery=100;
    t0=clock;
    specavekxin=cell(size(tspec));
    specavekyin=specavekxin;
    specavekxout=specavekxin;
    specavekyout=specavekxin;
    kspec=round(tspec/dt)+1;
    iout=1;
    for k=1:length(t)
        slice=squeeze(seis(k,:,:));
        [Slice,kx,ky]=fktran(slice,x,y,nx2,ny2);
        kky2=ky(ones(size(kx)),:).^2;
        kkx2=kx(:,ones(size(ky))).^2;
        sigmax=sigmaxt(k)*knx;
        sigmay=sigmayt(k)*kny;
        if(sigmax>=10)
            gx=ones(size(Slice));
        else
            gx=exp(-kkx2/sigmax^2);
        end
        if(sigmay>=10)
            gy=ones(size(Slice));
        else
            gy=exp(-kky2/sigmay^2);
        end
        if(sigmax<=10 || sigmay<=10)
            Slice2=Slice.*gx.*gy;
            tmp=ifktran(Slice2,kx,ky);
            slice2=tmp(1:nx,1:ny);
            slice2=slice2*norm(slice)/norm(slice2);
            seis(k,:,:)=shiftdim(slice2,-1);
        end
        
        if(rem(k,ievery)==0)
            tnow=clock;
            timeused=etime(tnow,t0);
            timeperslice=timeused/k;
            timeleft=timeperslice*(nt-k);
            disp(['Finished time ' time2str(t(k)) ' after ' int2str(timeused) ' seconds'])
            disp(['Time remaining '  int2str(timeleft) ' seconds'])
        end
        if(k==kspec(iout))
            specavekxin{iout}=mean(abs(Slice),2);
            specavekyin{iout}=mean(abs(Slice),1);
            specavekxout{iout}=mean(abs(Slice2),2);
            specavekyout{iout}=mean(abs(Slice2),1);
            if(iout<length(kspec))
                iout=iout+1;
            end
        end
    end
else
    %2D case
    nx2=2^nextpow2(nx);
    dkx=1/nx2;
    kx=-.5:dkx:.5-dkx;
    knx=.5;
    t0=clock;
    for k=1:length(t)
        tmp=fftshift(fft(seis(k,:),nx2));

        sigmax=sigmaxt(k)*knx;
        if(sigmax>=10)
            gx=ones(size(tmp));
        else
            gx=exp(-kx.^2/sigmax^2);
        end
        if(sigmax<=10)
            tmp2=tmp.*gx;
            tmp=ifft(fftshift(tmp2));
            seis(k,:)=tmp(1:nx);
        end
    end
end

totaltime=etime(clock,t0);
disp(['Finished tvgaussianmask after ' int2str(totaltime) ' seconds or ' num2str(totaltime/60) ' minutes'])
if(ny>1)
    specinfo={specavekxin,specavekyin,specavekxout,specavekyout,kx,ky};
else
    specinfo=[];
end

    



