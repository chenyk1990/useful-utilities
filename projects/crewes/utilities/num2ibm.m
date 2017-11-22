function u = num2ibm(d,lims,warn)
% function u = num2ibm(d,lims)
% Where:
%     u = IBM 4-byte floats stored as a 32-bit unsigned integers.
%         Bits can be examined using dec2hex(u)
%     d = Any number(s) stored in any Matlab datatype. vectors and arrays
%         are OK.
%  lims = 'ibm' or [], u is forced into the range -IBM_MAX:IBM_MAX
%       = 'ieee' u is forced into the range -IEEE_MAX:IEEE_MAX
%  warn = 0 or [] no warnings (default)
%       = 1 text warnings
%       = 2 fatal error
%
% NOTE: IBM 4-byte floats are encoded using the formula
%          ( -1)^sign * 0.fraction * 16^(ibm_exponent +ibm_bias),
%       where ibm_bias is defined to be +64.
%
% The sign is stored in bit 32 (the most significant bit of the most
%     significant byte)
% The exponent+bias is stored in bits 25-31 (the most significant byte)
% The fraction is stored in bits 1-24 (three least significant bytes)
%
% WARNING! 
%   IBM_MAX = (1.0-16.0^-6.0)*16.0^63.0 =~ 7.2370e+75 (+) max
%
%   **NOTE*: NaN and Inf are not defined for IBM floats
%   NaN is output as IBM_MAX with a warning
%   Inf is output as IBM_MAX with a warning
%   d > IBM_MAX is output as IBM_MAX with a warning
%
% Example:
%
% fid = fopen('test_ibm.sgy','w','ieee-be')
% u = num2ibm([realmin('single') realmax('single')])
% fwrite(fid,u,'uint32')
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

if nargin<2 || isempty(lims)
    lims='ibm';
end
if nargin<3 || isempty(warn)
    warn=0;
end

d = double(d);

switch lims    
    case 'ieee'
        FP_MAX     = double(realmax('single')); % =~ 3.4028e+38
        FP_MIN     = double(realmin('single')); % =~ 1.1755e-38
    otherwise
        FP_MAX     = (1-16^-6)*16^63;           % =~ 7.2370e+75
        FP_MIN     = 16^-65;                    % =~ 5.3976e-79
end
IBM_BIAS = 64.0; %defined

% get sign
s = sign(d);  %range -1:1
s(s>0.0)=0.0; %set +1 to 0
s=abs(s);     %set -1 to +1; sign range is now 0:1

% Deal with numbers greater than FP_MAX (includes +Inf)
idx = d > FP_MAX;
if sum(sum(idx))
    got_a_problem('Input data contains numbers greater than FP_MAX, setting to FP_MAX',warn);
    d(idx) = FP_MAX; %set to +FP_MAX
end

% Deal with numbers less than -FP_MAX (includes -Inf)
idx = d < -FP_MAX;
if sum(sum(idx))
    got_a_problem('Input data contains numbers less than -FP_MAX, setting to -FP_MAX',warn);
    d(idx) = FP_MAX; %set to -FP_MAX
    s(idx) = 1;
end

% Deal with numbers between 0.0 and FP_MIN
idx = abs(d)<FP_MIN;
if sum(sum(idx))
    got_a_problem('Input data contains numbers less than FP_MIN, setting to 0.0',warn);
    d(idx) = 0.0; %set to 0.0
end

% Deal with IEEE NaN
idx = isnan(d);
if sum(sum(idx))
    got_a_problem('Input data contains NaN, setting to FP_MAX',warn);
    d(idx) = FP_MAX; %set to +FP_MAX
    s(idx) = 0.0;
end    

% get IBM fraction and exponent
[f,e] = log16(abs(d));

%encode IBM sign, fraction and exponent
% s should be an integer stored as a double
% e +IBM_BIAS should be an integer stored as a double
% f is rounded by uint32(), not truncated (this is desirable)
u = uint32(s*2.0^31.0) + uint32((e +IBM_BIAS)*2.0^24.0) +uint32(f*2.0^24.0);

% Deal with hard zeros
idx = d==0.0; %do NOT use isequal(); true for +0 AND -0
if sum(sum(idx))
    u(idx)=0; %set to uint32(0)
end

end %end function num2ibm

function got_a_problem(txt,warn)
    switch warn
        case 1
            warning(txt)
        case 2
            error(txt)
        otherwise
    end
end %end function got_a_problem

