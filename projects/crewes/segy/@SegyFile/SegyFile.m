classdef SegyFile < File
%function sf = SegyFile(filename,permission,segyrev,sampint,nsamps,...
%                       segfmt,txtfmt,bytord,bindef,trcdef,gui)
%
% Required Inputs:
% filename   - SEG-Y disk file name
%
% Optional Inputs ([] is accepted):
% permission - 'r' (default), 'w', or 'a' (see help fopen)
% trcrange   - Traces to read from SEG-Y file. eg. 1:2:10.
%              Default ([]) attempts to read all traces in the file
% segyrev    - segy revision (0,1,2); Overrides SEG-Y Binary File Header on
%              disk.
% segfmt     - SEG-Y trace data format code. eg. 1 = IBM float, 5 = IEEE
%              float
% nsamps     - Samples per trace; Overrides SEG-Y Binary File Header on
%              disk.
% sampint    - Sample interval in (s); Overrides SEG-Y Binary File Header on
%              disk.
% bytord     - byte order of disk file 'l'=little-endian, 'b'=big-endian
% bindef     - 4 column binary header definition cell array such as provided by
%              @BinaryHeader/new; See uiSegyDefinition().
%              NOTE! writesegy will require the same bindef unless you
%              modify binhdr!
% trcdef     - 5 column trace header definition cell array such as provided by
%              @BinaryHeader/new; See uiSegyDefinition()
%              NOTE! writesegy will require the same trcdef unless you
%              modify trchdr!
% gui        - 0 (no progress bar), 1 (text progress bar), 
%              [] (default; gui progress bar and warnings), figure handle 
%              (same as [], but an attempt is made to center GUI popups on 
%              the figure represented by the figure handle)
% Outputs:
% sf        - A SEG-Y file object
%
% Example:
%   s=SegyFile('file.sgy','r')
%   s.TxtHdr.read
%   s.BinHdr.read
%   s.Trc.read(1:2:100)
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

properties
    FileInfo;  %contains struct
    TxtHdr;    %contains TextHeader object
    BinHdr;    %contains BinaryHeader object
    ExtTxtHdr; %contains ExtendedTextHeader object
    Trc;       %contains Trace object
end % end properties

methods
    %Constructor
    function obj = SegyFile(varargin)
        narginchk(1,11);
        
        filename = varargin{1};

        if nargin < 2 || isempty(varargin{2})
            permission='r';
        else
            permission = varargin{2};
        end

        %Doctor file permission so we can read info as well as write
        switch permission
            case 'w'
                permission = 'w+';
            case 'a'
                permission = 'a+';
        end
        
        if nargin < 3 %segyrevision
            segyrev=[];
        else
            segyrev = varargin{3};
        end
        if nargin < 4 % sampint
            sampint = [];
        else
            sampint = varargin{4};
        end 
        if nargin < 5 % nsamps
            nsamps = [];
        else
            nsamps=varargin{5};
        end    
        
        if nargin < 6 %fmtcode
            fmtcode = [];
        else
            fmtcode=varargin{6};
        end
        if nargin < 7 %fmtcode
            txtfmt = [];
        else
            txtfmt=varargin{7};
        end

        if nargin < 8 %bytord
            bytord = [];
        else
            bytord=varargin{8};
        end
        if nargin < 9 %bindef
            bindef = [];
        else
            bindef=varargin{9};
        end
        if nargin < 10 %trcdef
            trcdef = [];
        else
            trcdef=varargin{10};
        end
        if nargin < 11 %gui
            gui = 1;
        else
            gui = varargin{11};
        end
        
        %add .sgy filename extension if it does not already exist
        f = fliplr(filename);
        if strcmp(permission,'w+') && ~(strncmpi (f,'yges.',5) || strncmpi (f,'ygs.',4))
            filename = [filename '.sgy'];            
        end
        
        %Call @File constructor
        obj = obj@File(filename,permission,bytord,gui);
        if obj.FileID < 1 %File.openFile failed
            error(['@SegyFile: Unable to open file ' filename ' for ' permission]);
        end

        if isempty(bytord)
            obj = obj.whatEndian();
        end
        
        if(obj.Debug),disp('SegyFile: start of constructor');end

        %Create binary file header object based on disk file
        if(obj.Debug),disp('Creating BinHdr');end
        obj.BinHdr = BinaryHeader(obj.FileID,obj.Permission,...
            obj.ByteOrder,[],obj.GUI);
        
        %Override Binary File Header on disk before propagating to other
        %objects
        fileinfo = obj.BinHdr.FileInfo;
        if ~isempty(segyrev)
            fileinfo.SegyRevision=segyrev;
        end
        if ~isempty(fmtcode)
            fileinfo.FormatCode=fmtcode;
        end
        if ~isempty(nsamps)
            fileinfo.SamplesPerTrace=nsamps;
            
        end
        if ~isempty(sampint)
            fileinfo.SampleInterval=sampint*1e6; %microseconds
        end
        if ~isempty(bindef)
            obj.BinHdr.HdrDef = bindef;
        end
        obj.BinHdr.FileInfo=fileinfo;
        obj.FileInfo=fileinfo;        
        
        %Create text file header object
        if(obj.Debug),disp('Creating TxtHdr');end
        obj.TxtHdr = TextHeader(obj.FileID,obj.Permission,...
            obj.ByteOrder,obj.FileInfo.SegyRevision,obj.GUI);
        %Override guessed Text file header format
        if ~isempty(txtfmt)
            obj.TxtHdr.TxtFormat=txtfmt;
        end
        
        %Create extended text file header object
        if(obj.Debug),disp('Creating ExtTxtHdr');end
        obj.ExtTxtHdr = ExtendedTextHeader(obj.FileID,obj.Permission,...
            obj.ByteOrder,obj.FileInfo.SegyRevision,...
            obj.FileInfo.NumExtTxtHeaders,obj.GUI);
        %Override guessed Text file header format
        
        if ~isempty(txtfmt)
            obj.ExtTxtHdr.TxtFormat=txtfmt;
        end        
  
        %Create trace object
        if(obj.Debug),disp('Creating Trc');end
        
        if obj.FileInfo.TraceOneOffset==0
            obj.FileInfo.TraceOneOffset = ...
                obj.TxtHdr.SIZE...
                +obj.BinHdr.SIZE...
                +obj.TxtHdr.SIZE*double(obj.FileInfo.NumExtTxtHeaders);
        end
        
        obj.Trc = Trace(obj.FileID,obj.Permission,...
            nsamps,fmtcode,obj.ByteOrder,obj.FileInfo,obj.GUI);
        if~isempty(trcdef)
            obj.Trc.HdrDef=trcdef;
        end
        
        addlistener(obj,'FileNameChanged',@obj.listenFileName);
        addlistener(obj,'ByteOrderChanged',@obj.listenByteOrder);
        addlistener(obj,'PermissionChanged',@obj.listenPermission);
        
        if(obj.Debug),disp('SegyFile: end of constructor');end
    end
    
    function set.TxtHdr(obj, v)
        if isa(v,'TextHeader')
            obj.TxtHdr = v;
        else
            error('@SegyFile: TxtHdr must be a TextHeader object')
        end
    end
    
    function set.BinHdr(obj, v)
        if isa(v,'BinaryHeader')
            obj.BinHdr = v;
        else
            error('@SegyFile: BinHdr must be a BinaryHeader object')
        end
    end
    
    function set.ExtTxtHdr(obj, v)
        if isa(v,'ExtendedTextHeader')
            obj.ExtTxtHdr = v;
        else
            error('@SegyFile: ExtTxtHdr must be a ExtendedTextHeader object')
        end
    end
    
    function set.Trc(obj, v)
        if isa(v,'Trace')
            obj.Trc = v;
        else
            error('@SegyFile: Trc must be a Trace object')
        end
    end

    %Listeners
    function obj = listenFileName(obj, varargin)
        if(obj.Debug),disp('in File.listenFileName'); end
        if ~isempty(obj.FileName)
            obj.freopen();
            obj.TxtHdr.FileID = obj.FileID;
            obj.BinHdr.FileID = obj.FileID;
            obj.ExtTxtHdr.FileID = obj.FileID; 
            obj.Trc.FileID = obj.FileID;            
        end
    end
    
    function obj = listenByteOrder(obj, varargin)
        if(obj.Debug),disp('in File.listenByteOrder'); end
        if ~isempty(obj.FileName)
            obj.freopen();
            obj.TxtHdr.FileID = obj.FileID;
            obj.BinHdr.FileID = obj.FileID;
            obj.ExtTxtHdr.FileID = obj.FileID; 
            obj.Trc.FileID = obj.FileID;
        end
        obj.TxtHdr.ByteOrder = obj.ByteOrder;
        obj.BinHdr.ByteOrder = obj.ByteOrder;
        obj.ExtTxtHdr.ByteOrder = obj.ByteOrder;
        obj.Trc.ByteOrder = obj.ByteOrder;
    end
    
    function obj = listenPermission(obj, varargin)
        if(obj.Debug), disp('in File.listenPermission'); end
        if ~isempty(obj.FileName)
            obj.freopen();
            obj.TxtHdr.FileID = obj.FileID;
            obj.BinHdr.FileID = obj.FileID;
            obj.ExtTxtHdr.FileID = obj.FileID; 
            obj.Trc.FileID = obj.FileID;
        end
        obj.TxtHdr.Permission = obj.Permission;
        obj.BinHdr.Permission = obj.Permission;
        obj.ExtTxtHdr.Permission = obj.Permission;
        obj.Trc.Permission = obj.Permission;
    end
        
end % end methods

end % end classdef