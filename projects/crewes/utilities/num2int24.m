function u8 = num2int24(d,fe)
% function u = num2int24(d,fe)
% Where:
%     u = 3-byte signed integers stored as 8-bit unsigned integers
%     d = Any number(s) stored in any Matlab datatype. vectors and arrays
%         are OK
%    fe = file endian 'l' (little-endian) or 'b' (big-endian). Default is
%         computer endian.
%
% Example:
%
% fid = fopen('test.int24,'w')
% u = num2int24([-7,-8,-9],'b');
% fwrite(fid,u,'uint8')
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

% convert values to 4-byte integers
u32 = int32(d);

% byte swapping
if ~isequal(lower(e),lower(fe));
    u32=swapbytes(u32);
    e = fe;    
end

u8 = typecast(u32,'uint8');
% remove the fourth byte from each value
switch lower(e)
    case 'l'
        u8(4:4:end)=[];
    case 'b'
        u8(1:4:end)=[];
end

end %end function