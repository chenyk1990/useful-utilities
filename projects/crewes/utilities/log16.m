function [f,e] = log16(x)

%function [f,e] = log16(x)
% Emulates the behaviour of the Matlab built-in function log2, but using 
% base 16.
%
%   x = log16(d) returns the base 16 logarithm(s) of d, where d can be a real 
%       array of any datatype (integer, single, double).
%
%   [f,e] = log16(d) returns the (f)raction and (e)xponent that satisfy
%       d = f.*16.^e such as required to encode IBM floating point numbers.
%       Note that [f,e]=log16(d) does NOT check if d is suitable for
%       encoding as IBM floating point (ie, d>IBM_MAX, d==NaN, d==Inf)
%
%   See also 'help log2'
%
% Kevin Hall, CREWES, 2017
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

LOGCONST = 1.0/log10(16.0); %divide once
x=double(x); %make sure we are playing with doubles

if nargout<2
    %return log16(x)
    f=log10(x)*LOGCONST; %multiply many times (may be faster than division)
    return
end

s = sign(x);
x = abs(x);

e = ceil(log10(x)*LOGCONST); %multiply many times (may be faster than division)
% f1 = s.*x./(16.^e) %divide
f = s.*x.*16.0.^-e; %or multiply (may be faster than division)

% x==power of 16 gives a fraction of exactly 1, but for IBM floats the
% fraction must be less than 1.
idx = (f==1);
e(idx) = e(idx)+1; %increment exponent
f(idx) = s(idx).*x(idx)./(16.0.^e(idx)); %recalculate fraction

% x==0.0
idx = x==0.0; %do NOT use isequal(x,0) !!!
e(idx) = 0.0;
f(idx) = 0.0;

% x==Inf
idx = isinf(x);
e(idx) = 0.0;
f(idx) = s(idx).*Inf; %need s, as this could be +Inf or -Inf

% x ==NaN
idx = isnan(x);
e(idx) = 0.0;
f(idx) = NaN;




