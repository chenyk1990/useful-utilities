classdef TextHeader < File
%
%classdef TextHeader
%
%Class to deal with SEG-Y textual file headers
%
% Usage: th = TextHeader(fid)
%   where fid is a valid file identifier provided by fopen
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.
%
% Authors: Kevin Hall 2009, 2017
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

properties (Hidden = false)
    OFFSET = 0; %bytes from BOF
    SIZE   = 3200.0; %size in bytes
end

properties
    TxtFormat    = 'ascii';
    SegyRevision = 1;  
end

methods %public methods
    
    %Constructor
    function obj = TextHeader(filename,permission,byteorder,segyrevision,gui)
        if nargin <1 || isempty(filename)
            filename = -1;
        end
        if nargin <2 || isempty(permission)
            permission = [];
        end
        if nargin <3 || isempty(byteorder)
            byteorder = [];
        end
        if nargin <4 || isempty(segyrevision)
            segyrevision = 1;
        end
        if nargin <5
            gui = 1;
        end
        
        %Call superclass constructor
        obj = obj@File(filename,permission,byteorder,gui);
        if obj.FileID > 0 %File.openFile failed
%             error(['@TextHeader: Unable to open file ' filename ' for ' permission]);
%         else
            obj = obj.guessTextFormat;
        end
    
        obj.SegyRevision = segyrevision;
    end
    
    %Set methods

    function set.SegyRevision(obj,v)
        if isnumeric(v)
            obj.SegyRevision=v;
        else
            error('@TextHeader: SegyRevision must be numeric')
        end
    end
    
    function set.TxtFormat(obj,v)
        if ischar(v)
            v=lower(v);
            switch(v)
                case('ascii')
                    obj.TxtFormat = v;
                case('ebcdic')
                    obj.TxtFormat = v;
                otherwise
                    error('@TextHeader: TxtFormat must be ''ascii'' or ''ebcdic''');
            end
        else
            error('@TextHeader: TxtFormat must be char')
        end
    end    
    
end % end methods
    
end