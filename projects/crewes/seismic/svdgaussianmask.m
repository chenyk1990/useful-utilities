function seis=svdgaussianmask(seis,t,x,y,tspec,sigmax,sigmay,sigmaxB,sigmayB,singcut,flag)
%
% seis=svdgaussianmask(seis,t,x,y,tspec,sigmax,sigmay,sigmaxB,sigmayB,singcut,flag)
%
% For each time-slice of a 3D dataset, apply SVD separation to separate the slice into Gross and
% Detail (see svd_sep) and then apply wavenumber filters to both Gross and Detail before recombining. 
% The wavenumber filters are generally different for Gross and Detail usually being more severe for
% the former. The filters are applied as Gaussian masks in the wavenumber domain with each mask
% being specified by its standard deviations in x and y. The five parameters of the process can all
% be time-variant. The normal Gaussian mask is a decaying exponential but a growing one, a spatial
% deconvolution, is also allowed.
%
% seis ... input 3D dataset. First dimension is time, second ix x (xline) and third is y (iline).
%       Must be resularly sampled in all three dimensions. 
% t ... time coordinate vector, length(t) must equal size(seis,1).
% x ... x coordinate vector, length(x) must equal size(seis,2).
% y ... y coordinate vector, length(y) must equal size(seis,3).
% tspec ... vector of times at which the spatial filter sizes are specified.
% sigmax ... vector of standard deviations of the Gaussian mask to be applied to Gross in the x
%       direction expressed as a fraction of x spatial Nyquist, eg a value of .25 means .25*kxnyq.
%       length(tspec) must equal length(sigmax). Negative values for sigmax signify that the mask is
%       division by a Gaussian which causes spatial deconvolution (sharpening) while positive values
%       cause spatial convolution (blurring). Generally, sigmax should be a number between -1 and
%       +1. However, values greater than 1 are allowed and the larger it is the less action. Values
%       of 10 or greater are taken as a flag for no action and such time slices are simply passed to
%       output.
% sigmay ... Just like sigmax except for the y direction. x and y filtering actions can be
%       completley independent.
% sigmaxB and sigmayB ... Just like sigmax and sigmay except they are for filtering Detail.
% singcut ... vector of singular value cuttoffs, the same length as tspec
% flag ... 0 means the singular-value cutoff function is a bimodal boxcar
%          1 means the singular-value cutoff function is a Gaussian
% *********** default flag=1 **************
% 
% seis (output) ... 3D seismic matrix the same size as the input with the spatial filtering applied
%
% NOTE: Because we are effectively convolving each time slice with a different size spatial filter,
% there is a possible time-variant amplitude effect. At present, this is addressed by balancing
% input to output rms power on each time slice.
%
% NOTE2: To determine the input parameters for this process, view your data in plotimage3D. View a
% time-slice of interest and right-click in to to launch the wavenumber filtering tool.
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

if(nargin<11)
    flag=1;
end

disp('svdgaussianmask')
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
if(length(tspec)~=length(sigmaxB))
    error('tspec and sigmax must have the same length')
end
if(length(tspec)~=length(sigmayB))
    error('tspec and sigmay must have the same length')
end
if(length(tspec)~=length(singcut))
    error('tspec and singcut must have the same length')
end

tspec=tspec(:)';
sigmax=sigmax(:)';
sigmay=sigmay(:)';
singcut=singcut(:)';

%expand filter specification to cover the entire time range
ind=find(tspec==t(1), 1);
if(isempty(ind))
    tspec=[t(1) tspec];
    sigmax=[sigmax(1) sigmax];
    sigmay=[sigmay(1) sigmay];
    sigmaxB=[sigmaxB(1) sigmaxB];
    sigmayB=[sigmayB(1) sigmayB];    
    singcut=[singcut(1) singcut];
end
ind=find(tspec==t(end), 1);
if(isempty(ind))
    tspec=[tspec t(end)];
    sigmax=[sigmax sigmax(end)];
    sigmay=[sigmay sigmay(end)];
    sigmaxB=[sigmaxB sigmaxB(end)];
    sigmayB=[sigmayB sigmayB(end)];    
    singcut=[singcut singcut(end)];
end

sigmaxt=interp1(tspec,sigmax,t);
sigmayt=interp1(tspec,sigmay,t);
sigmaxtB=interp1(tspec,sigmaxB,t);
sigmaytB=interp1(tspec,sigmayB,t);
singcutt=interp1(tspec,singcut,t);


ievery=100;
t0=clock;

for k=1:length(t)
    slice=squeeze(seis(k,:,:));
    slicesvd=svd_wavenumber_filt(slice,x,y,sigmaxt(k),sigmayt(k),sigmaxtB(k),sigmaytB(k),singcutt(k),flag);
    seis(k,:,:)=shiftdim(slicesvd,-1);
    
    if(rem(k,ievery)==0)
        tnow=clock;
        timeused=etime(tnow,t0);
        timeperslice=timeused/k;
        timeleft=timeperslice*(nt-k);
        disp(['Finished time ' time2str(t(k)) ' after ' int2str(timeused) ' seconds'])
        disp(['Time remaining '  int2str(timeleft) ' seconds'])
    end

end

totaltime=etime(clock,t0);
disp(['Finished svdgaussianmask after ' int2str(totaltime) ' seconds or ' num2str(totaltime/60) ' minutes'])

    



