function tf=check_fp_range(d,datafmt)
%function tf=check_fp_range(d,datafmt)
%
% Where:
%   tf = true or false
%   d = data to check
%   datafmt = ibm32  - 4-byte IBM float
%           = single - 4-byte float (assumed to be IEEE)
%           = ieee32 - 4-byte IEEE float (see single)
%           = double - 8-byte float (assumed to be IEEE
%           = ieee64 - 8-byte IEEE float (see double)
%
% Usage:
%   Call this function before converting d to an single or a double
%
% Example:
%   if ~check_fp_range(d,'single')
%      error ('data is out of range')
%   end
%
% If d contains NaN, Inf or any numbers outside of:
% fmin = realmin(datafmt)
% fmax = realmax(datafmt)
%      -fmax>d || -fmin<d<0.0 || 0.0<d<fmin || fmax<d
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
    case 'ieee32'
        fmin = realmin('single');
        fmax = realmax('single');
    case 'ieee64'
        fmin = realmin('double');
        fmax = realmax('double');
    case 'ibm32'
        fmin = 16.0^-65.0;                % =~ 5.3976e-79 (+) 4-byte IBM min
        fmax = (1.0-16.0^-6.0)*16.0^63.0; % =~ 7.2370e+75 (+) 4-byte IBM max
    otherwise
        fmin = realmin(datafmt);
        fmax = realmax(datafmt);
end

if sum(sum(d<-fmax) +sum(d<0.0&d>-fmin) +sum(d>0.0&d<fmin) +sum(d>fmax))
    tf = false;
end



end %end function check_int_range