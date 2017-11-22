function d = uint24_2num(u,fe)
% function d = uint24_2num(u,fe)
% Where:
%     u = 3-byte signed int stored as 8-bit unsigned integers
%     d = Any number(s) stored in any Matlab datatype. vectors and arrays
%         are OK
%    fe = file endian 'l' (little-endian) or 'b' (big-endian)
%
% Example:
%
% fid = fopen('test.uint24,'r')
% u = fread(fid,nbytes,'uint8=>uint8')
% d = int242num(u,'b');
% fclose(fid)
%
% Authors: Kevin Hall, 2017
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.
%

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

narginchk(1,2)

[~,~,e]=computer;

if nargin<2 || isempty(fe)
    fe=e;
end

if ~isa(u,'uint8') || ~isvector(u)
    error('input vector must be data type ''uint8''')
end

inlen = length(u);
outlen = inlen+inlen/3;

c = mod(inlen,3);
if c
    error ('input array must have a length that is a multiple of three');
end

u8  = zeros(1,outlen,'uint8');
idx = ones(1,outlen,'logical');

switch lower(fe)
    case 'l'
        idx(4:4:end)=false;
    case 'b'
        idx(1:4:end)=false;
end

u8(idx)=u;
u32 = typecast(u8,'uint32');

if ~isequal(lower(e),lower(fe))
    u32 = swapbytes(u32);
end

d = double(u32);

end %end function