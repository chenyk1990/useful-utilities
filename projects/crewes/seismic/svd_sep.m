function [seisg,seisd,singg,singd]=svd_sep(seis,singcut,dim,flag)
% SVD_SEP: singular value decomposition separation into Gross and Detail
%
% [seisg,seisd]=svd_sep(seis,singcut,dim,flag)
%
% Those unfamiliar with SVD (singular value decomposition) should examine the function svd first. A
% second step is to view your matrix with seisplotsvd_sep in order to determine the input parameter
% singcut. The basic idea here is to use SVD to separate a matrix into Gross and Detail parts. These
% are matrices the same size and class as the input where Gross contains the large-scale features
% and Detail contains the small-scale features. The input matrix is recreated by the simple
% summation of Gross and Detail. This is a 2D operation so if the input matrix is 3D then the
% separation is done on each "slice" where the dimension "dim" (usually the first) is held constant.
% The input can be single or double precision but the SVD is always done at double precision and the
% result is returned in the same precision as the input. For each slice, an SVD (see svd) is
% performed on the slice to divide it into two components: sliceg (Gross) is dominated by the
% largest singular values while sliced (Detail) is dominanted by all of the others. Adding sliceg
% and sliced reproduces the input slice to machine accuracy. This decomposition is achieved by
% defining a Gaussian of the form g(j)=exp(-(j-1)^2/jcut^2) where j is the singular value number
% (j=1:nsing where nsing is the number of singular values and j=1 is the largest -most important-
% singular value). Let the singular value decomposition of the input slice be slice=U*S*V' where S
% contains the singular values on its diagonal (again see help on svd). Then let s be the diagonal
% entries of S (all other entries are zero) and sliceg is defined by sliceg=U*diag(s.*g)*V' where
% diag(s.*g) is the diagonal matrix like S except that the singular values are windowed by g.
% Defining d=1-g, sliced is given by sliced=U*diag(s.*d)*V' and it follows that sliceg+sliced=slice
% and sliceg and sliced have the mentioned properties. Given this separation, many further
% possibilities exist. For example, see svd_wavenumber_filt and svdgaussianmask.
%
% seis ... input seismic matrix, either 2D or 3D
% singcut ... cutoff singular value defining the separation between Grodd and Detail. This is jcut
%           in the description above.
% dim ... dimension held constant for the separation. Has no effect if seis is 2D.
% *********** default = 1 *************
% flag ... 0 means the singular-value cutoff function is a bimodal boxcar
%          1 means the singular-value cutoff function is a Gaussian
%
% seisg ... matrix of the same size and class as seis containing the Gross structure of seis
% seisd ... matrix of the same size and class as seis containing the Detail of seis
% singg ... singular values for Gross. Each column of singg contains the singular values for the
%           correswponding slice of seisg. There is nothing useful that can be done with these
%           without the U-V matrices. They are just provided for interest.
% singd ... like singg except for Detail.
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
    dim=1;
end

sz=size(seis);
nt=sz(1);
nx=sz(2);

if(length(sz)==3)
    ny=sz(3);
    if(dim==1)
        tmp=double(squeeze(seis(1,:,:)));
        nslices=nt;
    elseif(dim==2)
        tmp=double(squeeze(seis(:,1,:)));
        nslices=nx;
    else
        tmp=double(squeeze(seis(:,:,1)));
        nslices=ny;
    end
   
elseif(length(sz)==2)
    nslices=1;
    tmp=double(seis);
else
    error('svd_sep works only on 2D or 3D matrices');
end

[m,n]=size(tmp);
nsing=min([m n]);

singcut=round(singcut);

if(singcut<1 || singcut>nsing)
    error(['singcut must be between 1 and ' int2str(nsing) ' inclusive'])
end

seisg=seis;
seisd=zeros(size(seis),'like',seis);
clear seis;

singg=zeros(nsing,nslices);
singd=singg;

ievery=round(nslices/100);

tnot=clock;
if(nslices>1)
    disp('svd_sep')
end
for k=1:nslices
    
    [U,S,V]=svd(tmp);
    singvals=diag(S);
    if(k==1)
        j=1:nsing;
        if(flag==1)
            g=exp(-(j-1).^2/singcut^2)';
        else
            g=zeros(size(j))';
            ind=j<=singcut;
            g(ind)=1;
        end
        d=1-g;
    end
    tmpsingg=singvals.*g;
    tmp2=diag(tmpsingg);
    if(m>n)
        Sg=[tmp2;zeros(m-n,n)];
    elseif(n>m)
        Sg=[tmp2 zeros(m,n-m)];
    else
        Sg=tmp2;
    end
    tmpsingd=singvals.*d;
    tmp=diag(tmpsingd);
    if(m>n)
        Sd=[tmp;zeros(m-n,n)];
    elseif(n>m)
        Sd=[tmp zeros(m,n-m)];
    else
        Sd=tmp;
    end
    sliceg=U*Sg*V';
    sliced=U*Sd*V';
    %sliced=tmp-sliceg;
    if(nslices==1)
        seisg=sliceg;
        seisd=sliced;
        singg=tmpsingg(:);
        singd=tmpsingd(:);
    else
        if(dim==1)
            seisg(k,:,:)=shiftdim(sliceg,-1);
            seisd(k,:,:)=shiftdim(sliced,-1);
        elseif(dim==2)
            seisg(:,k,:)=reshape(sliceg,m,1,n);
            seisd(k,:,:)=shiftdim(sliced,m,1,n);
        else
            seisg(:,k,:)=reshape(sliceg,m,n,1);
            seisd(k,:,:)=shiftdim(sliced,m,n,1);
        end
       singg(:,k)=tmpsingg(:);
       singd(:,k)=tmpsingd(:);
    end
    if(k<nslices)
        if(dim==1)
            tmp=double(squeeze(seis(k+1,:,:)));
        elseif(dim==2)
            tmp=double(squeeze(seis(:,k+1,:)));
        else
            tmp=double(squeeze(seis(:,:,k+1)));
        end
    end
    
    if(rem(k,ievery)==0)
        tnow=clock;
        timeused=etime(tnow,tnot);
        timeperslice=timeused/k;
        timeremaining=timeperslice*(nslices-k);
        disp(['finished slice ' int2str(k) ' of ' int2str(nslices)]);
        disp(['elapsed time ' int2str(timeused) ', remaining time ' int2str(timeremaining) ' s']);
    end
        
end


