function write ( obj, bh )
%
% function write ( obj, bh )
%
% Write a binary file header structure to a SEG-Y file
%
% WARNING: Error checking is limited to cross-referencing the struct to 
% obj.HdrDef. It is the users responsibility to verify values contained in
% the struct before writing to disk
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

%expecting exactly two arguments, obj and bh
narginchk(2,2)

%make sure bh is a struct
if ~isstruct(bh)
    error('@BinaryHeader/write: Header to write to disk must be a struct')
end

%check bh struct
obj.check(bh);

%make sure text header has already been written
fs = obj.fsize();
if isequal(fs,0)
    error('@BinaryHeader/write: File on disk does not contain a text file header');
elseif fs >0 && fs <obj.OFFSET
    error(['@BinaryHeader/write: File on disk text header is not ' num2str(obj.OFFSET) ' bytes']);
elseif fs >obj.OFFSET
    error('@BinaryHeader/write: Refusing to overwrite existing binary file header');
end

%position pointer at start of binary header
obj.fseek(obj.OFFSET,'bof');

%Write header value
for ii = 1:length(obj.HdrDef)
    v = bh.(obj.HdrDef{ii,1});    
    obj.fwrite(v,obj.HdrDef{ii,2});
end
        
end