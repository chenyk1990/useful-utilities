function varargout = read(obj,trcrange,whattoread)
%
% function varargout = read(obj,trcrange,whattoread)
%
% Where:
% trcrange (optional)   = 1:obj.FileInfo.TracesInFile (default)
% whattoread (optional) = 'both' (default)
%                           Returns [trcdata,trchdr] where trcdata is a
%                           matrix containing the seismic data stored as
%                           per Trace.FormatCode (1 => IBM converted to
%                           'single', 2 = 'int32', 3 = 'int16', 4 is
%                           unsupported, 5 is 'single', 8 is 'int8'), and
%                           trchdr is a structure of arrays that is 
%                           constructed using Trace.HdrDef
%                       = 'headers'
%                           Returns trchdr (see 'both')
%                       = 'data'
%                           Returns trcdata (see 'both')
%                       = A fieldname from Trace.HdrDef column 1
%                           Returns a vector of the datatype stored in
%                           Trace.HdrDef column 2
%
% Examples:
%   >> sf = SegyFile('test.sgy','r')
%
% Read all trace headers and trace data in SEGY file (Next four examples 
% have exactly the same result):
%   >> [trcdata,trchdr] = sf.Trc.read();
%   >> [trcdata,trchdr] = sf.Trc.read(1:sf.Trc.FileInfo.TracesInFile);
%   >> [trcdata,trchdr] = sf.Trc.read(1:sf.Trc.FileInfo.TracesInFile,'both');
%   >> [trcdata,trchdr] = sf.Trc.read([],'both',true);
%
% Read all trace headers (Next two examples do exactly the same thing):
%   >> trchdr  = sf.Trc.read(1:sf.Trc.FileInfo.TracesInFile,'header');
%   >> trchdr  = sf.Trc.read([],'header');
%
% Read all trace data (Next two examples do exactly the same thing):
%   >> trcdata = sf.Trc.read(1:sf.Trc.FileInfo.TracesInFile,'data');
%   >> trcdata = sf.Trc.read([],'data');
%
% Read all values for a given header word (Next two examples do exactly 
% the same thing):
%   spt = sf.Trc.read(1:sf.Trc.FileInfo.TracesInFile,'SourcePoint');
%   spt = sf.Trc.read([],'SourcePoint',);
%
% Read some trace headers and trace data
%   trchdr  = sf.Trc.read(100:2:200,'header');
%   trcdata = sf.Trc.read(100:2:200,'data');
%   spt     = sf.Trc.read(100:2:200,'SourcePoint');
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

%check inputs
narginchk(1,4)

if nargin == 1
    trcrange = 1:obj.FileInfo.TracesInFile;
    whattoread = 'both';
    gui = 1;
end

if nargin == 2
    trcrange = checkTrcRange(obj,trcrange);    
    whattoread = 'both';
    gui = 1;
end

if nargin == 3    
    trcrange = checkTrcRange(obj,trcrange);
    if isempty(whattoread)
        whattoread = 'both';    
    elseif ~ischar(whattoread)
        error('@Trace/read: whattoread must be empty or char');
    end
    gui = 1;
end

if nargin == 4
    trcrange = checkTrcRange(obj,trcrange);
    if isempty(whattoread)
        whattoread = 'both';    
    elseif ~ischar(whattoread)
        error('@Trace/read: whattoread must be empty or char');
    end
    if ~isnumeric(gui)
        error('@Trace/read: gui must be numeric');
    end
end

%convert trace numbers to byte locations from start of file
trcrange=obj.FileInfo.TraceOneOffset+(trcrange-1)*obj.FileInfo.TraceSize;

switch(whattoread)
    case 'both'
        [varargout{1},varargout{2}] = obj.readBoth(trcrange);
    case 'headers'
        varargout{1} = obj.readHeader(trcrange);
        varargout{2} = [];
    case 'data'
        varargout{1} = obj.readData(trcrange);
        varargout{2} = [];
    otherwise
        varargout{1} = obj.readHeaderWord(trcrange,whattoread);
        varargout{2} = [];
end %end switch(whattoread)

end %end function read

function trcrange = checkTrcRange(obj,trcrange)

maxtrace = obj.FileInfo.TracesInFile;

if isempty(trcrange)
    trcrange = 1:maxtrace;
else
    if ~isnumeric(trcrange)
        error('@Trace/read: trcrange must be numeric');
    end
    %check trace range to make sure it is in range
    if(min(trcrange)<1)
        error('@Trace/read: Minimum trace number in file is 1');
    end
    
    if(max(trcrange)>maxtrace)
        error(['@Trace/read: Maximum trace number in file is ',...
            num2str(maxtrace)]);
    end
end

end % end function checkTrcRange