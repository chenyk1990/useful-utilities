function thead = read ( obj )
%
%function thead = read ( obj )
%
% Read the textual file header from a SEG-Y file
% Returns:
%   SEG-Y textual file header as a 2D char matrix
%
% Warning:
%   This function assumes the text header on disk is obj.TxtFormat ('ascii'
%   or 'ebcdic') and converts as appropriate
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

if isequal(obj.NumHdrs,0)
    thead = [];
else
    error('@ExtendedTextHeader/read: Extended Text Headers read not yet implemented')
end
% %make sure we are positioned at the start of the header
% fseek(obj.FileID,obj.HdrOffset,'bof');
% %read the header
% thead = fread(obj.FileID, obj.HDRSIZE, 'char=>char');
% 
% if ~isequal(length(thead),obj.HDRSIZE)
%     error('@ExtendedTextHeader/read: Failed to read textual header');
% end
% 
% %convert header vector to 2D matrix
% thead = reshape(thead,80,40)';
% %convert ebcdic to ascii if needed
% switch(obj.TxtFormat)
%     case 'ebcdic'
%         thead = obj.ebcdic2ascii(thead);
% end

end %end function