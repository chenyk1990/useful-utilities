function write ( obj, thead )
%
%function write ( obj, thead )
% Warning: Input textual header is assumed to be a 2D ASCII matrix
%   Conversions are performed, this function writes:
%     EBCDIC if obj.TxtFormat = 'ebcdic'
%     ASCII  if obj.TxtFormat = 'ascii'
%
% Authors: Chad Hogan, 2004
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

if ~isempty(thead)
    warning ('@ExtendedTextHeader/write: Extended Text Headers write not yet implemented')
end

% %make sure we're at the start of the header in the file
% fseek(obj.FileID,obj.HdrOffset,'bof');
% 
% %convert ascii to ebcdic if needed
% switch(obj.TxtFormat)
%     case('ebcdic')
%         thead = obj.ascii2ebcdic(thead);
% end   
%     
% %write text header to file
% c=fwrite(obj.FileID, thead', 'uchar');
% 
% %check to make sure the fwrite was successful
% if ~isequal(c,obj.HDRSIZE)
%     error('@ExtendedTextHeader/write: Failed to write textual header to file');
% end

end