classdef BinaryHeader <File
    %
    %classdef BinaryHeader
    %
    % Class to deal with SEG-Y binary file headers
    %
    % Usage: bh = BinaryHeader(filename,permission,byteorder,gui,varargin)
    %
    %   Where filename can be char or fid from fopen
    %
    % Authors: Kevin Hall, 2009, 2017
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
        SIZE   = 400.0;  %bytes
        OFFSET = 3200.0; %bytes
    end
    
    properties
        FileInfo;
        HdrDef = {}; %cell array containing the binary file header definition
       % SegyRevision = 1;
    end
    
    properties (Dependent)
        HdrFieldNames; %HdrDef col 1
        HdrDataTypes;  %HdrDef col 2
        HdrStartBytes; %HdrDef col 3
        HdrLongNames;  %HdrDef col 4
    end
    
    events
        SegyRevisionChanged;
    end
    
    methods
        function obj = BinaryHeader(filename,permission,byteorder,...
                fileinfo,gui,varargin)
            
            if nargin <1 || isempty (filename)
                filename = -1;
            end
            if nargin <2
                permission = [];
            end
            if nargin <3
                byteorder = [];
            end
            if nargin <4
                fileinfo = [];
            end            
            if nargin <5
                gui = 1;
            end
            
            %Call superclass constructor
            obj = obj@File(filename,permission,byteorder,gui,varargin{:});
                       
            if isempty(fileinfo)
                obj.FileInfo = obj.readFileInfo();
            else
                obj.FileInfo = fileinfo;
            end

            obj.HdrDef = obj.newDefinition(obj.FileInfo.SegyRevision);
            
        end %end constructor        
        
        function set.FileInfo(obj,v)
            if isstruct(v)
                obj.FileInfo=v;
            else
                error('@BinaryHeader: FileInfo must be a struct')
            end
        end
        
        function set.HdrDef(obj,v)
            if iscell(v)
                obj.check(v)
                obj.HdrDef=v;
            else
                error('@BinaryHeader: HdrDef must be a valid cell array; see BinaryHeader/checkDefinition')
            end
        end
        
        function fn = get.HdrFieldNames(obj)
            %returns contents of col 1 of obj.HdrDef
            fn = obj.HdrDef(:,1);
        end %end get.HdrFieldNames
        
        function st = get.HdrDataTypes(obj)
            %converts cols 1 and 2 of obj.HdrDef into a struct
            st = cell2struct(obj.HdrDef(:,2),obj.HdrDef(:,1));
        end %end get.HdrDataTypes
        
        function st = get.HdrStartBytes(obj)
            %converts cols 1 and 3 of obj.HdrDef into a struct
            st = cell2struct(obj.HdrDef(:,3),obj.HdrDef(:,1));
        end %end get.HdrStartBytes
        
        function st = get.HdrLongNames(obj)
            %converts cols 1 and 4 of obj.HdrDef into a struct
            idx = cellfun(@(X) ~isempty(X),obj.HdrDef(:,4));
            st = cell2struct(obj.HdrDef(idx,4),obj.HdrDef(idx,1));
        end %end get.HdrLongName
        
        %Functions to read specific header words independent of obj.HdrDef
        function hw = readFormatCode(obj)
            hw = readHeaderWord(obj,3224,'uint16');
        end
        
        function hw = readSamplesPerTrace(obj)
            hw = readHeaderWord(obj,3220,'uint16');
        end
        
        function hw = readTracesPerRec(obj)
            ndathw = readHeaderWord(obj,3212,'uint16');
            nauxhw = readHeaderWord(obj,3214,'uint16');
            if isequal(ndathw,nauxhw)
                hw=ndathw;
            else
                hw = ndathw+nauxhw;
            end
        end
        
        function hw = readSampleInterval(obj)
            hw = readHeaderWord(obj,3216,'uint16');
        end
 
        function hw = readExtTracesPerRec(obj)
            ndathw = readHeaderWord(obj,3260,'uint32');
            nauxhw = readHeaderWord(obj,3264,'uint32');
            hw = ndathw+nauxhw;
        end
 
        function hw = readExtSamplesPerTrace(obj)
            hw = readHeaderWord(obj,3268,'uint32');
        end
        
        function hw = readExtSampleInterval(obj)
            hw = readHeaderWord(obj,3272,'ieee64');
        end
        
        function hw = readIntegerConstant(obj)
            hw = readHeaderWord(obj,3296,'uint32');
        end          
        
        function hw = readSegyRevision(obj)
            segmaj = double(readHeaderWord(obj,3500,'uint8'));
            segmin = double(readHeaderWord(obj,3501,'uint8'));
            hw = segmaj+segmin/10;
            
            %hw _should be_ 0, 1, 2
            if hw < 0 || hw > 2
                warning(['Unknown SEG-Y revision ' num2str(hw) ...
                    ': Attempting to continue using SEG-Y revision 0'])
                hw=0;
            end
        end
 
        function hw = readFixedTrcLength(obj)
            hw = readHeaderWord(obj,3502,'uint16');
        end
        
        function hw = readNumExtHeaders(obj)
            hw = readHeaderWord(obj,3504,'int16');
        end
        
        function hw = readNumExtTrcHeaders(obj)
            hw = readHeaderWord(obj,3506,'uint32');
        end                

        function hw = readTraceOneOffset(obj)
            hw = readHeaderWord(obj,3520,'uint64');
        end        
 
        function hw = readNumDataTrailers(obj)
            hw = readHeaderWord(obj,3528,'int32');
        end        
        
        function hw = readHeaderWord(obj,startbyte,datatype)
            %function readHeaderWord(obj,startbyte,datatype)
            obj.fseek(startbyte,'bof');
            hw = obj.fread(1, datatype);
        end
        
        %Listeners
        function obj = listenSegyRevision(obj, varargin)
            %disp('in listenSegyRevision')
            obj.HdrDef = obj.newDefinition(obj.SegyRevision);
        end
    end %end methods
    
    methods (Static)
        bd = newDefinition(segyrev);
    end
    
end %end classdef