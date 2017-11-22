classdef Trace < File
    %classdef Trace
    %
    % Class for SEGY Trace Headers
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
        SIZE    = 240.0;  %bytes
    end
    
    properties
        FileInfo          = [];
        HdrDef            = {};
        ApplyCoordScalars = true; % logical true or false. true applies scalars to header words
    end
    
    properties
        SegyRevision      = 1;
    end
    
    properties (Dependent)
        HdrFieldNames; %HdrDef col 1
        HdrDataTypes;  %HdrDef col 2
        HdrStartBytes; %HdrDef col 3
        HdrScalars;    %HdrDef col 4
        HdrLongName;   %HdrDef col 5
        FormatCodeType;
        BytesPerSample;
    end
    
    methods
        function obj = Trace(filename,permission,nsamps,segfmt,byteorder,fileinfo,gui)
            if nargin <1 || isempty(filename)
                filename = -1;
            end
            if nargin <2 || isempty(permission)
                permission = [];                                              
            end
            if nargin <2 || isempty(nsamps)
                nsamps = [];                                              
            end
            if nargin <2 || isempty(segfmt)
                segfmt = [];                                              
            end            
            if nargin <3 || isempty(byteorder)
                byteorder = [];
            end
            if nargin <4 || isempty(fileinfo)
                fileinfo = [];
            end
            if nargin <5
                gui = 1;
            end
            
            %Call superclass constructor
            obj = obj@File(filename,permission,byteorder,gui);
%             if obj.FileID < 1 %File.openFile failed
%                 error(['@Trace: Unable to open file ' filename ' for ' permission]);
%             end

            if isempty(fileinfo)
                fileinfo.SegyRevision       = 1;
                fileinfo.FormatCode         = 5; % 4-byte IEEE
                fileinfo.SamplesPerTrace    = 0;
                fileinfo.TracesPerRec       = 0;
                fileinfo.SampleInterval     = 1000; % 1 ms
                fileinfo.FixedTrcLength     = 0; %fixed trace length is true
                fileinfo.NumExtTxtHeaders   = 0;
                fileinfo.ExtTracesPerRec    = 0;
                fileinfo.ExtSamplesPerTrace = 0;
                fileinfo.ExtSampleInterval  = 0;
                fileinfo.IntegerConstant    = 0;
                fileinfo.NumExtTrcHeaders   = 0;
                fileinfo.TraceOneOffset     = 0;                               
            end
            
            obj.FileInfo = fileinfo;                       
            obj.HdrDef = obj.newDefinition(obj.FileInfo.SegyRevision);
            
%             if obj.fsize < obj.FileInfo.TraceOneOffset+obj.SIZE && ~strncmpi(obj.Permission,'r',1)
%                 error('@Trace: File does not appear to contain any data. Is it a Shortcut or Link?');
                       
            %elseif obj.FileID>0 && obj.fsize > obj.FileInfo.TraceOneOffset+obj.SIZE
            if obj.FileID>0 && obj.fsize > obj.FileInfo.TraceOneOffset+obj.SIZE
                if isempty(segfmt)
                    obj.guessFormatCode();
                end
                
                obj.FileInfo.TraceSize = double(obj.FileInfo.SamplesPerTrace)...
                    *double(obj.BytesPerSample)...
                    +obj.SIZE;
                obj.FileInfo.TracesInFile = ...
                    (obj.fsize()-obj.FileInfo.TraceOneOffset)/obj.FileInfo.TraceSize;
            else
                obj.FileInfo.TraceSize = 0;
                obj.FileInfo.TracesInFile = 0;
            end
            
            if obj.FileInfo.TracesInFile-fix(obj.FileInfo.TracesInFile)>0
                
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    a = mm_yesnodlg(['File may not contain fixed length traces! '...
                        'Assume fixed length and attempt to continue?'],...
                        'Warning!',...
                        'No',obj.GUI);
                elseif obj.GUI
                    disp('File may not contain fixed length traces!')
                    a = input('Assume fixed length and attempt to continue (y/n)? ','s');
                else
                    a='n';
                end
                
                if strncmpi(a,'n',1)
                    obj.FileID = -1;
                    return
                else
                    obj.FileInfo.TracesInFile = floor(obj.FileInfo.TracesInFile);
                end
            else
                obj.FileInfo.FixedTrcLength = 1;
            end            
        end        
        
        function set.HdrDef(obj,v)
            if iscell(v)
                obj.check(v);
                obj.HdrDef=v;
            else
                error('@Trace: HdrDef must be a valid cell array')
            end
        end
        
        function set.ApplyCoordScalars(obj,v)
            if islogical(v)
                obj.ApplyCoordScalars=v;
            else
                error('@Trace: ApplyCoordScalars must be logical (true/false)')
            end
        end
        
        function set.FileInfo(obj,v)
            if isstruct(v)
                obj.FileInfo=v;
            else
                error('@Trace: FileInfo must be a struct')
            end
        end
        
        %Get functions
        function nb = get.BytesPerSample(obj)
            switch obj.FileInfo.FormatCode
                case 1 %ibm floating point
                    nb = 4.0;
                case 2 %4-byte int
                    nb = 4.0;
                case 3 %2-byte int
                    nb = 2.0;
                    %             case 4 %4-byte fixed-point with gain (obsolete)
                    %                 nb = 4.0;
                case 5 %4-byte IEEE
                    nb = 4.0;
                case 6 %8-byte IEEE
                    nb = 8.0;
                case 7 %3-byte int
                    nb = 3.0;
                case 8 %1-byte int
                    nb = 1.0;
                case 9 %8-byte int
                    nb = 8.0;
                case 10 %uint32
                    nb = 4.0;
                case 11 %uint16
                    nb = 2.0;
                case 12 %uint64
                    nb = 8.0;
                case 15 %uint24
                    nb = 3.0;
                case 16 %uint8
                    nb = 1.0;
                otherwise
                    error(['@Trace: FormatCode ''' num2str(obj.FileInfo.FormatCode) ''' not supported ']);
            end
        end %end get.bytesPerSample
        
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
        
        function st = get.HdrScalars(obj)
            %converts cols 1 and 4 of obj.HdrDef into a struct
            idx = cellfun(@(X) ~isempty(X),obj.HdrDef(:,4));
            st = cell2struct(obj.HdrDef(idx,4),obj.HdrDef(idx,1));
        end %end get.HdrScalars
        
        function st = get.HdrLongName(obj)
            %converts cols 1 and 5 of obj.HdrDef into a struct
            st = cell2struct(obj.HdrDef(:,5),obj.HdrDef(:,1));
        end %end get.HdrLongName
        
        function dtype = get.FormatCodeType(obj)
            %for use with fread/fwrite
            switch obj.FileInfo.FormatCode
                case 1 %ibm floating point
                    dtype = 'ibm32';
                case 2 %4-byte int
                    dtype = 'int32';
                case 3 %2-byte int
                    dtype = 'int16';
%                 case 4 %4-byte fixed-point with gain (obsolete)
%                     dtype = 'fixed32';
                case 5 %4-byte IEEE
                    dtype = 'ieee32';
                case 6 %8-byte IEEE
                    dtype = 'ieee64';
                case 7 %3-byte int
                    dtype = 'int24';
                case 8 %1-byte int
                    dtype = 'int8';
                case 9 %8-byte int
                    dtype = 'int64';
                case 10 %uint32
                    dtype = 'uint32';
                case 11 %uint16
                    dtype = 'uint16';
                case 12 %uint64
                    dtype = 'uint64';
                case 15 %uint24
                    dtype = 'uint24';
                case 16 %uint8
                    dtype = 'uint8';
                otherwise
                    error(['@Trace: FormatCode ''' num2str(obj.FileInfo.FormatCode) ''' not supported ']);
            end
        end %end get.FormatCodeType
        
        function st = applyCoordinateScalars(obj, st)
            %apply coordinate scalars
            if obj.ApplyCoordScalars
                fnames = fieldnames(obj.HdrScalars);
                for ii = 1:length(fnames)
                    if isfield(st,fnames{ii}) && isfield(st,obj.HdrScalars.(fnames{ii}))
                        hw = single(st.(fnames{ii})); %get header word values
                        sc = single(st.(obj.HdrScalars.(fnames{ii}))); %get scalar values
                        
                        ps_idx = sc > 0;
                        ns_idx = sc < 0;
                        
                        hw(ps_idx) = hw(ps_idx)./sc(ps_idx); %divide if scalar is positive
                        hw(ns_idx) = hw(ns_idx).*abs(sc(ns_idx)); %multiply if scalar is negative
                        
                        st.(fnames{ii})=hw;
                    end
                end
            end
        end
        
        function st=removeCoordinateScalars(obj, st)
            %remove coordinate scalars
            if obj.ApplyCoordScalars
                fnames = fieldnames(obj.HdrScalars);
                for ii = 1:length(fnames)
                    if isfield(st,fnames{ii}) && isfield(st,obj.HdrScalars.(fnames{ii}))
                        hw = double(st.(fnames{ii})); %get header word values
                        sc = double(st.(obj.HdrScalars.(fnames{ii}))); %get scalar values
                        
                        ps_idx = sc > 0;
                        ns_idx = sc < 0;
                        
                        hw(ps_idx) = hw(ps_idx).*sc(ps_idx); %multiply if scalar is positive
                        hw(ns_idx) = hw(ns_idx)./abs(sc(ns_idx)); %divide if scalar is negative
                        
                        st.(fnames{ii})=hw;
                    end
                end
            end
        end
        
    end%end methods
    
    methods (Hidden)
        
        function [trcdata,trchead]=readBoth(obj,trcrange)
            
            ntrace = length(trcrange);
            %Yeah. Nice try-don't do this: ToDo; a subroutine that works
            %for custom HdrDef
%             [trcdata,trchead] = obj.new(ntrace,obj.FileInfo.SamplesPerTrace); %pre-allocate memory
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
%                 pos=get(gcf,'position');%get current figure position
%                 xc=pos(1)+.5*pos(3);%center of current figure
%                 yc=pos(2)+.5*pos(4);%center of current figure
                h = waitbar(0,['Reading ' num2str(ntrace) ' Trace(s): '],...
                    'Name','SEG-Y Input',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%                 posw=get(h,'position');%position of waitbar
%                 x0=xc-.5*posw(3);%adjusted position of waitbar
%                 y0=yc-.5*posw(4);%adjusted position of waitbar
%                 set(h,'position',[x0 y0 posw(3:4)]);%adjust position of waitbar
                setappdata(h,'canceling',0);
                mm_adjust(h);
            end
            for ii = 1:ntrace
                
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    %check for cancel button press
                    if getappdata(h,'canceling')
                        delete(h);
                        break
                    end
                    %update graphics progress bar
                    waitbar(ii/ntrace,h);
                    
                elseif obj.GUI
                    %update text progress bar
                    textbar(['Reading ' num2str(ntrace) ' Trace(s): '],...
                        30,ii,ntrace);
                end
                
                %fseek to start of trace
                obj.fseek(trcrange(ii),'bof');
                
                %read trace header
                for jj = 1:length(obj.HdrDef)
                    fieldname = obj.HdrDef{jj,1};
                    datatype = obj.HdrDef{jj,2};
                    trchead.(fieldname)(ii) = obj.fread(1,datatype);
                end
                
                %read trace data
                v = obj.fread(obj.FileInfo.SamplesPerTrace,obj.FormatCodeType);
                
                if strcmp(obj.FormatCodeType,'ibm32')                    
                    v = single(v); %convert to single-precision
                    idx = isinf(v); %numbers abs(v) > realmax('single') are now Inf
                    if sum(sum(isinf(v)))
                        v(idx)=sign(v(idx)).*realmax('single');
                    end                    
                end
                trcdata(:,ii) = v;
            end
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
                %close waitbar
                delete(h);
            end
            
            %remove coordinate scalars
            trchead = obj.removeCoordinateScalars(trchead);
        end
        
        function trchead = readHeader(obj,trcrange)
            ntrace   = length(trcrange);

            %Yeah. Nice try-don't do this: ToDo; a subroutine that works
            %for custom HdrDef
%             trchead  = obj.new(ntrace,obj.FileInfo.SamplesPerTrace,[],'headers'); %pre-allocate memory
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
%                 pos=get(gcf,'position');%get current figure position
%                 xc=pos(1)+.5*pos(3);%center of current figure
%                 yc=pos(2)+.5*pos(4);%center of current figure
                h = waitbar(0,['Reading Headers From ' num2str(ntrace) ' Trace(s): '],...
                    'Name','SEG-Y Input',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%                 posw=get(h,'position');%position of waitbar
%                 x0=xc-.5*posw(3);%adjusted position of waitbar
%                 y0=yc-.5*posw(4);%adjusted position of waitbar
%                 set(h,'position',[x0 y0 posw(3:4)]);%adjust position of waitbar
                setappdata(h,'canceling',0);
                mm_adjust(h);
            end
            
            for ii = 1:ntrace
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    %check for cancel button press
                    if getappdata(h,'canceling')
                        delete(h);
                        break
                    end
                    %update graphics progress bar
                    waitbar(ii/ntrace,h);
                elseif obj.GUI
                    %update text progress bar
                    textbar(['Reading Headers From ' num2str(ntrace) ' Trace(s): '],...
                        30,ii,length(trcrange))
                end
                
                %fseek to start of trace
                obj.fseek(trcrange(ii),'bof');
                
                %read trace header
                for jj = 1:length(obj.HdrDef)
                    fieldname = obj.HdrDef{jj,1};
                    datatype = obj.HdrDef{jj,2};
                    trchead.(fieldname)(ii) = obj.fread(1,datatype);
                end
            end
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
                %close waitbar
                delete(h);
            end
            
            %remove coordinate scalars
            trchead = obj.removeCoordinateScalars(trchead);
        end
        
        function trcdata = readData(obj,trcrange)
            ntrace = length(trcrange);
            trcdata = obj.new(ntrace,obj.FileInfo.SamplesPerTrace,[],'data'); %pre-allocate memory
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
%                 pos=get(gcf,'position');%get current figure position
%                 xc=pos(1)+.5*pos(3);%center of current figure
%                 yc=pos(2)+.5*pos(4);%center of current figure
                h = waitbar(0,['Reading Data From ' num2str(ntrace) ' Trace(s): '],...
                    'Name','SEG-Y Input',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%                 posw=get(h,'position');%position of waitbar
%                 x0=xc-.5*posw(3);%adjusted position of waitbar
%                 y0=yc-.5*posw(4);%adjusted position of waitbar
%                 set(h,'position',[x0 y0 posw(3:4)]);%adjust position of waitbar
                setappdata(h,'canceling',0);
                mm_adjust(h);
            end
            
            for ii = 1:ntrace
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    %check for cancel button press
                    if getappdata(h,'canceling')
                        delete(h);
                        break
                    end
                    %graphics progress bar
                    waitbar(ii/ntrace,h);
                elseif obj.GUI
                    %text progress bar
                    textbar(['Reading Data From ' num2str(ntrace) ' Trace(s): '],...
                        30,ii,ntrace)
                end
                
                %fseek to start of trace data
                obj.fseek(trcrange(ii)+obj.SIZE,'bof');
                
                %read trace data
                v = obj.fread(obj.FileInfo.SamplesPerTrace,obj.FormatCodeType);
                
                if strcmp(obj.FormatCodeType,'ibm32')                    
                    v = single(v); %convert to single-precision
                    idx = isinf(v); %numbers abs(v) > realmax('single') are now Inf
                    if sum(sum(isinf(v)))
                        v(idx)=sign(v(idx)).*realmax('single');
                    end                    
                end
                trcdata(:,ii) = v;
            end
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
                %close waitbar
                delete(h);
            end
        end
        
        function hw = readHeaderWord(obj,trcrange,whattoread)
            hw = readHeaderWordIgnoreScalars(obj,trcrange,whattoread);
            
            if obj.ApplyCoordScalars %remove coordinate scalars
                %Get header word scalar name that should be applied
                try
                    hwsname = obj.HdrScalars.(whattoread);
                catch ex
                    %most likely we're here because the header word does not need to have a scalar applied
                    %disp(ex.message)
                    return
                end
                
                sc = readHeaderWordIgnoreScalars(obj,trcrange,hwsname);
                ps_idx = sc > 0;
                ns_idx = sc < 0;
                
                hw=double(hw);
                sc=double(sc);
                
                hw(ps_idx) = hw(ps_idx).*sc(ps_idx); %multiply if scalar is positive
                hw(ns_idx) = hw(ns_idx)./abs(sc(ns_idx)); %divide if scalar is negative
            end
        end
        
        function hw=readHeaderWordIgnoreScalars(obj,trcrange,whattoread)
            
            %check to see if whattoread is a valid field
            if ~isfield(obj.HdrDataTypes,whattoread)
                error(['@Trace: Header word ''' whattoread ''' not found in obj.HdrDef']);
            end
            
            datatype = obj.HdrDataTypes.(whattoread);
            hwoffset = obj.HdrStartBytes.(whattoread);
            
            %pre-allocate memory
            ntrace = length(trcrange);
            hw=obj.newHeaderWord(length(trcrange),datatype);
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
%                 pos=get(gcf,'position');%get current figure position
%                 xc=pos(1)+.5*pos(3);%center of current figure
%                 yc=pos(2)+.5*pos(4);%center of current figure
                h = waitbar(0,['Reading Header Word ''' whattoread ''' from ' num2str(ntrace) ' Trace(s): '],...
                    'Name','SEG-Y Input',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%                 posw=get(h,'position');%position of waitbar
%                 x0=xc-.5*posw(3);%adjusted position of waitbar
%                 y0=yc-.5*posw(4);%adjusted position of waitbar
%                 set(h,'position',[x0 y0 posw(3:4)]);%adjust position of waitbar
                setappdata(h,'canceling',0);
                mm_adjust(h);
            end
            
            for ii = 1:ntrace
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    %check for cancel button press
                    if getappdata(h,'canceling')
                        delete(h);
                        break
                    end
                    %graphics progress bar
                    waitbar(ii/ntrace,h);
                elseif obj.GUI
                    %text progress bar
                    textbar(['Reading Header Word From ' num2str(ntrace) ' Trace(s): '],...
                        30,ii,ntrace)
                end
                
                %fseek to start of header word
                obj.fseek(trcrange(ii)+hwoffset,'bof');
                
                %read header word
                hw(ii) = obj.fread(1,datatype);
            end
            
            if isempty(obj.GUI) || isa(obj.GUI,'handle')
                %close waitbar
                delete(h);
            end
        end
        
        function th = newHeader(obj,ntrace,nsamp,sampint)
            %sampint in seconds
            if isempty(ntrace)
                ntrace=1;
            end
            if isempty(nsamp)
                nsamp=1;
            end           
            if isempty(sampint)
                sampint = obj.FileInfo.SampleInterval; %microseconds
            else
                sampint = sampint*1e6; %microseconds
            end

            segyrev = obj.FileInfo.SegyRevision;

            hdrdef = obj.newDefinition(segyrev);
            for ii = 1:length(hdrdef)
                th.(hdrdef{ii,1}) = obj.newHeaderWord(ntrace,hdrdef{ii,2});
            end

            hw = obj.byte2word(114); %Samples this trace; stored in byte 115; SEG-Y standard
            th.(hw)(:) = nsamp; %microseconds
            hw = obj.byte2word(116); %Sample rate this trace; stored in byte 115; SEG-Y standard
            th.(hw)(:) = sampint; %microseconds
        end
        
        function td = newData(obj,ntrace,nsamp)
            if isempty(ntrace)
                ntrace=1;
            end
            if isempty(nsamp)
                obj.FileInfo.SamplesPerTrace;
            end
            fct = obj.FormatCodeType;
            
            if isequal(fct,'ibm32')
                fct = 'single';
            elseif isequal(fct,'ieee32')
                fct = 'single';
            elseif isequal(fct,'ieee64')
                fct = 'double';
            elseif isequal(fct,'int24')
                fct = 'int32';
            elseif isequal(fct,'uint24')
                fct = 'uint32';
            end

            td = zeros(nsamp,ntrace,fct);
        end
        
    end %end methods (Hidden)
    
    methods (Static)
        td = newDefinition(segyrev,numexthdrs);
        hw = newHeaderWord(ntrace,datatype);
        [d,f,t] = struct2double(s);
        s = double2struct(d,f,t);        
    end % end static methods
    
    
end %end classdef