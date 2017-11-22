function [fn,bl] = byte2word(obj, byteloc)
%
% function [fn,bl] = byte2word(obj, byteloc)
%
% BYTE2WORD returns the header word nearest a byte location in traceheaders
%
% fn=byte2word(byteloc)
%
% byteloc ... a SEGY byte location. Must be between 1 and 240. 
% fn      ... the name of the SEGY header word in obj.HdrDef whose starting 
%             byte location is nearest to byteloc
% bl      ... actual byte location for fn
%
% Example, retrieve the SEGY header value stored at byte location 17
% [fn,bl] = byte2word(17)
%
% Authors: Gary Margrave,
%          Kevin Hall, 2017
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

if byteloc < 0  || byteloc > obj.SIZE
    error(['@BinaryHeader/byte2word: Byteloc must be between 0 and ' num2str(obj.SIZE)])
end

v = [obj.HdrDef{:,3}];   %convert header byte locations to a double vector
ind = near (v, byteloc); %find index containing closest match to byteloc
fn = obj.HdrDef{ind,1};  %return fieldname for ind
bl = obj.HdrDef{ind,3};  %return closest byte in header definition to byteloc

end %end function
