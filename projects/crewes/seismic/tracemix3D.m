function seis=tracemix3D(seis,weights)
% TRACEMIX3D ... weighted trace mixing on a 3D seismic volume
%
% seis=tracemix3D(seis,weights)
%
% seis ... input seismic section. Should be a 3D matrix with one trace per
%           column, t is dimenstion 1, x is dimension 2 and y is dimension 3.
% weights ... 2D matrix of mixing weights. (Can be a 1D vector if the weights are the same in x and y).
%           Mixing is accomplished as a 2D convolution with the time slices. Therefore the major
%           computation loop is over dimension 1. Dimension 1 of weights corresponds to dimension 2
%           of seis and dimension 2 of weights is dimension 3 of seis. The input weights are
%           normalized by division by sum(abs(weights(:))).
%
% seis ... output seismic section
%
% by G.F. Margrave, Devon, 2017
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

[nwx,nwy]=size(weights);
if(nwx==1 || nwy==1)
    n=max([nwx nwy]);
    tmp=weights(:);
    weights=tmp(:,ones(1,n));
end

a=sum(abs(weights(:)));
ievery=250;
for k=1:nt
    tmp=double(squeeze(seis(k,:,:)));
    tmp2=conv2(tmp,weights/a,'same');
    if(isa(seis,'single'))
        seis(k,:,:)=shiftdim(single(tmp2),-1);
    else
        seis(k,:,:)=shiftdim(tmp2,-1);
    end
    if(rem(k,ievery)==0)
        disp(['completed sample ' int2str(k)])
    end
end
    