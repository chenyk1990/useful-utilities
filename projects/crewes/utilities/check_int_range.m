function tf = check_int_range(d,datafmt)
%function check_int_range(d,datafmt)
%
% Where:
%   tf = true or false
%   d = data to check
%   datafmt = int8   - 1-byte integer
%           = int16  - 2-byte integer
%           = int24  - 3-byte integer
%           = int32  - 4-byte integer
%           = int64  - 8-byte integer
%           = uint8  - 1-byte unsigned integer
%           = uint16 - 2-byte unsigned integer
%           = uint24 - 3-byte unsigned integer
%           = uint32 - 4-byte unsigned integer
%           = uint64 - 8-byte unsigned integer
%
% Usage:
%   Call this function before converting d to an integer
%   format
%
% Example:
%   if ~check_int_range(d,'uint16')
%      error ('data is out of range')
%   end
%
% If d contains NaN, Inf or any numbers outside of: 
%   d < intmin(datafmt) || d > intmax(datafmt)
% this function will exit with an error
%
% Authors: Kevin Hall, 2017
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

tf = true; %assume true
if ~isnumeric(d)
    tf = false; return;
end

if isnan(d)
    tf = false; return;
end

if isinf(d)
    tf = false; return;
end

switch datafmt
    case 'int24';        
       imin = -8388608;
       imax = 8388607;
    case 'uint24';        
       imin = 0;
       imax = 16777215;
    otherwise
       imin = intmin(datafmt);
       imax = intmax(datafmt);
end

if sum(sum(d<imin) + sum(d>imax))
    tf = false;
end

end %end function check_int_range