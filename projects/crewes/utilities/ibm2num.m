function d = ibm2num(u)
%
% function d = ibm2num(u)
% Where:
%     u = IBM 4-byte float stored as a 32-bit unsigned integer
%     d = IBM stored value converted to double
%
% Example
%
% fid = fopen('test.ibm','r','ieee-be')
% u = fread(fid,'uint32=>uint32')
% d = ibm2num(u)
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

%https://en.wikipedia.org/wiki/...
%   IBM_Floating_Point_Architecture#Single-precision_32-bit
%   Accessed Feb, 2017
%
%IBM 4-byte floats are encoded using the formula:
% (-1)^sign * 0.fraction * 16^(ibm_exponent +ibm_bias)
% where ibm_bias is defined to be +64.
%
% The sign is stored in bit 32 (the most significant bit of the most
%     significant byte)
% The exponent+bias is stored in bits 25-31 (the most significant byte)
% The fraction is stored in bits 1-24 (three least significant bytes)
%

if ~isa(u,'uint32')
    error('ibm2num: Input must be a 4-byte unsigned integer');
end

%get sign bit 0=(+), 1=(-), stored in bit 32
sbit = double(bitget(u,32));

%get exponent
u = bitset(u,32,0); %remove sign
%hex2dec('ff000000') = 4278190080; hex2dec is slow!
eibm = bitand(u,uint32(4278190080)); %mask significand
eibm = double(bitshift(eibm,-24))-70.0; %bitshift exponent and unbias

%get fraction
%hex2dec('00ffffff') = 16777215; hex2dec is slow!
fibm = double(bitand(u,uint32(16777215))); %mask exponent

%Calculate answer
%d = (-1.0)^signbit * significand * 16.0^eibm;
d = fibm.*(-1.0).^sbit.*(16.0).^eibm;

end %end function d = ibm2num(u)