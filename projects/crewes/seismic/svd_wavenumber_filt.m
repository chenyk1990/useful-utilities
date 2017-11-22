function [slicesvd,more]=svd_wavenumber_filt(slice,x,y,sigmax,sigmay,sigmaxB,sigmayB,singcut,flag)
% SVD_WAVENUMBER_FILT: wavenumber filtering on the largest singular values
%
% [slicesvd,more]=svd_wavenumber_filt(slice,x,y,sigmax,sigmay,singcut)
%
% Supppression of the acquisition footprint is often addressed by wavenumber filtering of time
% slices designed to suppress higher wavenumbers. However, this causes a loss of detail in the
% result. This function provides a means to apply a wavenumber filter to the gross structure
% while preserving the detail. It has been observed that the footprint is mostly in the larger
% singular values therefore filtering the larger ones while preserving the smaller ones if the main
% idea. To get a better understanding, select a time slice from a 3D survey and view it with
% seisplotsvd before attempting this filtering. This function is inherently 2D. The 3D version of
% the same thing is svdgaussianmask. 
%
% Given an input time slice from a 3D survey, an SVD (see svd) is performed on the slice to divide it
% into two components: slicea is dominated by the largest singular values while sliceb is dominanted
% by all of the others. Adding slicea and sliceb reproduces the input slice to machine accuracy. This
% decomposition is achieved by defining a Gaussian of the form g(j)=exp(-(j-1)^2/jcut^2) where j is
% the singular value number (j=1:nsing where nsing is the number of singular values and j=1 is the
% largest -most important- singular value). Let the singular value decomposition of the input slice
% be slice=U*S*V' where S contains the singular values on its diagonal (again see help on svd). Then
% let s be the diagonal entries of S (all other entries are zero) and slicea is definded by
% slicea=U*diag(s.*g)*V' where diag(s.*g) is the diagonal matrix like S except that the singular
% values are windowed by g. Defining h=1-g, sliceb is given by sliceb=U*diag(s.*h)*V' and it follows
% that slicea+sliceb=slice and slicea and sliceb have the mentioned properties. After this
% decomposition, both slices are passed through wavenumber filtering by Gaussian mask (see
% wavenumber_gaussmask) to give slicea2 and slice b2. In general slicea should be filtered much more
% than sliceb. (Consider using seisplotsvd_foot to determine the parameters for this process.)  As a
% final step, scalars a and b are determined such that ||slicea-a*slicea2|| and
% ||sliceb-b*sliceb2||are minimal (least-squares subtraction).  The final result is formed as
% slicesvd=a*slicea2+b*sliceb2.
%
% slice ... input time slice, the column coordinate is x and the row coordinate is y.
% x ... coordinate for the rows of slice. The length(x) must equal size(slice,1). x should be
%       in physical units like feet or meters (i.e. not linenumber) so that wavenumbers are
%       calculated correctly. Must be regularly spaced for the FFT.
% y ... coordinate for the columns of slice. The length(y) must equal size(slice,2). y should be
%       in physical units like feet or meters (i.e. not linenumber) so that wavenumbers are
%       calculated correctly. Must be regularly spaced for the FFT.
% sigmax,sigmay ... standard deviations of the Gaussian mask to be applied to slicea. Expressed as
%       fractions of the corresponding Nyquist wavenumbers. These can be positive or negative
%       numbers. Positive means suppression, negative means amplification. See wavenumber_gaussmask
%       for more detail.
% sigmaxB,sigmayB ... standard deviations of the Gaussian mask to be applied to sliceb. Expressed as
%       fractions of the corresponding Nyquist wavenumbers. These can be positive or negative
%       numbers. Positive means suppression, negative means amplification. See wavenumber_gaussmask
%       for more detail.
% singcut ... the number of the "cutoff" singular value. Must be an integer between 1 and the total
%       number of singular values nsing. If slice has m rows and n columns then nsing=min([m n]).
%       singcut is the same as jcut mentioned in the method description. Roughly speaking, slicea is
%       dominated by singular values from 1 to singcut while slice be is dominated by those from
%       singcut to nsing. Viewing a representative slice with seisplotsvd is recommended prior to
%       choosing this value. Expect it to be a a smallinteger perhaps around 10.
% flag ... 0 means the singular-value cutoff function is a bimodal boxcar
%          1 means the singular-value cutoff function is a Gaussian
% *********** default flag=1 **************
%
% slicesvd ... the output svd-wavenumber filtered slice
% more is a cell array with more information
% more{1} = slicea ... the portion of the input slice which will be wavenumber filtered
% more{2} = sliceb ... the portion of the input slice which will NOT be wavenumber filtered
% more{3} = slicea2 ... the wavenumber filtered version of slicea (scaled by "a")
% more{4} = sliceb2 ... the wavenumber filtered version of sliceb (scaled by "b")
% more{5} = [a b] ... values of the scalars a and b
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
sz=size(slice);
if(length(sz)>2)
    error('svd_wavenumber_filt works only in 2D, see svdgaussianmask for 3D')
end
m=sz(1);
n=sz(2);
nsing=min([m n]);
if(nargin<9)
    flag=1;
end
if(length(x)~=m)
    error('x has the wrong size');
end
if(length(y)~=n)
    error('y has the wrong size');
end

small=10^6*eps;
if(sum(abs(diff(diff(x))))>small)
    error('x must be resularly spaced');
end

if(sum(abs(diff(diff(y))))>small)
    error('y must be resularly spaced');
end

if(singcut<1 || singcut>nsing)
    error(['singcut must be between 1 and ' int2str(nsing) ' inclusive'])
end

[slicea,sliceb]=svd_sep(slice,singcut,1,flag);

slicea2=wavenumber_gaussmask(slicea,y,x,sigmay,sigmax);
sliceb2=wavenumber_gaussmask(sliceb,y,x,sigmayB,sigmaxB);

[~,a]=lsqsubtract(slicea(:),slicea2(:));
slicea2=a*slicea2;
[~,b]=lsqsubtract(sliceb(:),sliceb2(:));
sliceb2=b*sliceb2;
slicesvd=slicea2+sliceb2;
if(nargout>1)
    more=cell(1,5);
    more{1}=sliceaq;
    more{2}=sliceb;
    more{3}=slicea2;
    more{4}=sliceb2;
    more{5}=[a b];
end