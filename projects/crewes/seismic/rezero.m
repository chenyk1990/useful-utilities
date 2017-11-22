function seis3D=rezero(seis3D,zeromap)
% REZERO: return empty traces in a 3D seismic matrix to zero
%
% seis3D=rezero(seis3D,zeromap)
%
% A 3D seismic survey (post-stack) will usually have an iregular footprint in the lateral
% coordinates. That is when viewed in map view, the survey boundary is not a perfect rectangle but
% rather is some irregular shape. In order to store it in a 3D matrix, it must be padded with zero
% traces out to the smallest rectangle that encloses the survey. seis3D should be a 3D seismic
% matrix whose first dimension is time, the second dimension is x, and the third dimension is y.
% zeromap is a 2D matrix whose row dimension corresponds to the seismic x and whose column dimension
% corresponds to the seismic y. Wherever zeromap is zero, the corresponding trace in seis will be
% set to zero. See the note below for an example of creating zeromap.
%
% seis3D (input) ... 3D survey with traces that need to be set to zero.
% zeromap ... 2D matrix inidcating which traces in seis3D should be set to zero. Dimension 1 of
%           zeromap must match dimension 2 of seis3D and dimension 2 of zeromap must match dimension
%           3 of seis3D. Wherever zeromap==0, the corresponding trace will be set to zero.
% seis3D (output) ... identical to the input seis3D except that the trace locations indicated in
%       zeromap are now also zero in seis3D.
%
% NOTE: The easiest way to create zeromap is to sum the last few samples of a same-size seismic
% volume in which the proper traces are already zero. For example, let seis3Dz be a a 3D seismic
% volume with the proper traces set to zero, and let seis3D be the result of processing seis3Dz in
% such a way that the zero traces become nonzero. For example, a spatial filter run on time slices
% (e.g. tvwavenfilt) will have this effect. Then, if nt is the number of time samples in seis3Dz
% (the size of dimension 1), the zero map can be created by
% zeromap=squeeze(sum(abs(seis3Dz(nt-10:nt,:,:))));
% The choice to sum the last 10 samples in each trace here is arbitrary but it does assume that if
% the final 10 samples of a trace are zero, then the entire trace is zero. This is usually the case
% but there are bound to be exceptions.
%
% G.F. Margrave, Devon Canada, 2017
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


sz1=size(seis3D);
sz2=size(zeromap);

if(length(sz1)~=3);
    error('seis3D must be a 3D matrix');
end
if(length(sz2)~=2)
    error('zeromap must be a 2D matrix');
end
if(sz1(2)~=sz2(1))
    error('dimension 2 of seis3D must match dimension 1 of zeromap');
end
if(sz1(3)~=sz2(2))
    error('dimension 3 of seis3D must match dimension 2 of zeromap');
end

nx=sz1(2);ny=sz1(3);
iz=0;
ntraces=nx*ny;
for k=1:nx
    for j=1:ny
        if(zeromap(k,j)==0)
            seis3D(:,k,j)=0;
            iz=iz+1;
        end
    end
end

disp(['REZERO: ' int2str(iz) ' traces set to zero from a total of ' int2str(ntraces) ' traces'])
