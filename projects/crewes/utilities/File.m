classdef File < handle
    %File class
    %  Open a file for file operations such as reading or writing
    %
    % Usage:
    %  File(filename,permission,byteorder,gui)
    %
    % Where:
    %   filename      = string containing full file name
    %   permission    = string cotaining file permissions to use with fopen
    %                     (optional), default is 'r' (see help fopen)
    %   byteorder     = string containing 'b' or 'l' for big- or little-endian
    %                     (optional), default is 'n' (see help fopen)
    %   gui           = gui=0: no prompts,
    %                       1: text prompts
    %         handle or empty: gui prompts. GUI calling this function should catch
    %                          errors and display using mm_errordlg;
    %
    % Example:
    %   sf = File('file.sgy','r','b')
    %   sf.fread(400, 'single') %system fread with error checking,
    %     equivalent to:
    %   fread(fid, 400, 'single=>single')
    %
    % Additional data types for fread and fwrite:
    %    'ibm32'  see ibm2num and num2ibm
    %    'ieee32' == 'single' with fopen() using 'ieee-le' or 'ieee-be'
    %    'ieee64' == 'double' with fopen() using 'ieee-le' or 'ieee-be'
    %    'int24'  see int24_2num and num2int24
    %    'uint24' see uint24_2num and num2uint24
    %
    %  Kevin Hall, 2009, 2017
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
        FileName=[];  % FileName corresponding to FileID
        FileID=[];    % File handle output from fopen
        Permission=[];% File Permission, eg. 'r', 'w', 'a'
        ByteOrder=[]; % Byte Order, 'l' for little endian, 'b' for big endian
        GUI=1;        % GUI flag: 0: no prompts, 1: text prompts (default), figure handle: gui prompts
    end
    
    properties (Hidden)
        Debug = false;
    end
    
    %events
    events
        FileNameChanged; % Notifies listeners that the FileName has changed
        ByteOrderChanged; % Notifies listeners that the ByteOrder has changed
        PermissionChanged; % Notifies listeners that the Permission has changed
        GUIchanged; % Notifies listeners that the GUI has changed
    end
    
    %methods
    methods
        
        %Constructor
        function obj = File(filename, permission, byteorder, gui)
            if(obj.Debug), disp('File: Start of constructor'); end
            % Constructor for File object
            % check number of input arguments
            narginchk(0,4)
            
            if nargin <1 || isempty(filename)
                filename=-1;
            end
            % Check inputs
            if nargin <2 || isempty(permission)
                permission='r';
            end
            if nargin <3 || isempty(byteorder)
                [~,~,e] = computer;
                byteorder = lower(e);
            end
            if nargin <4
                gui = true; %text warnings and prompts
            end
            
            obj.Permission = permission;
            obj.ByteOrder = byteorder;
            obj.GUI = gui;

            %set object values
            if ischar(filename)
                obj.FileName = filename; %string
            else
                obj.FileID = filename; %fid
            end
            
            % Open file
            if ~isempty(obj.FileName)
                obj = obj.fopen();
            end
            
            % Add Listeners
            addlistener(obj,'FileNameChanged',@obj.listenFileName);
            addlistener(obj,'ByteOrderChanged',@obj.listenByteOrder);
            addlistener(obj,'PermissionChanged',@obj.listenPermission);
        
            if(obj.Debug), disp('File: End of constructor'); end
        end %end constructor
        

        
        %Set functions
        function set.FileName(obj, v)
            if ischar(v)
                obj.FileName = v;
                
            else
                error('@File: FileName must be char');
            end
            notify(obj,'FileNameChanged');
        end
        
        function set.FileID(obj, v)
            if isnumeric(v)
                obj.FileID = v;
            else
                error('@File: FileID must be numeric');
            end
        end        
        
        function set.ByteOrder(obj, v)
            if ischar(v)
                switch lower(v)
                    case 'native' %- local machine format - the default
                        obj.ByteOrder='n';
                    case 'n'
                        obj.ByteOrder=v;
                    case'ieee-le' %- IEEE floating point with big-endian
                        obj.ByteOrder='l';
                    case 'l'
                        obj.ByteOrder=v;
                    case'ieee-be' %- IEEE floating point with big-endian
                        obj.ByteOrder='b';
                    case 'b'
                        obj.ByteOrder=v;
                    case'ieee-le.l64' %- IEEE floating point with little-endian
                        obj.ByteOrder='a';
                    case 'a'
                        obj.ByteOrder=v;
                    case'ieee-be.l64' %- IEEE floating point with big-endian byte
                        obj.ByteOrder='s';
                    case 's'
                        obj.ByteOrder=v;
                    otherwise
                        error(['@File: Unknown ByteOrder ''' v ''' see help fopen'])
                end
            elseif (isempty(v))
                obj.ByteOrder='n';
            else
                error('@File: ByteOrder must be empty or a character string')
            end            
            %disp('notifying obj that byteorder changed')
            notify(obj,'ByteOrderChanged');
        end
        
        function set.Permission(obj, v)
            if (ischar(v))
                if strcmp(v,'r') || strcmpi(v,'w') || strcmpi(v,'a') ||...
                        strcmp(v,'r+') || strcmp(v,'w+') || strcmp(v,'a+')
                    obj.Permission=v;
                    
                else
                    error(['@File: Permission ''' v ''' is invalid: see ''help fopen'''])
                end
            elseif (isempty(v))
                obj.Permission = 'r';
            else
                error('@File: Permission must be a character string')
            end
            notify(obj,'PermissionChanged');
        end
        
        function set.GUI(obj, v)
            if isnumeric(v)
                obj.GUI = logical(v);
                notify(obj,'GUIchanged');
            elseif islogical(v) || isempty(v) || isa(v,'handle')
                obj.GUI = v;
                notify(obj,'GUIchanged');                
            else
                error('@File: GUI must be numeric, empty, or a figure handle')
            end
        end
        
        %Open and close file
        function obj = fopen(obj)
            %'fopen()' with error checking
            if(obj.Debug),disp('in File.open');end
            if isempty(obj.FileName)
                %                 error('@File: Cannot open file without a filename');
                return
            end

            if strncmpi(obj.Permission,'w',1) && exist(obj.FileName,'file')
                if isempty(obj.GUI) || isa(obj.GUI,'handle')
                    a = mm_yesnodlg(...
                        ['File ''' obj.FileName ...
                        ''' exists and will be overwritten! '...
                        'Continue?'], ...
                        'Warning!','No');
                elseif islogical(obj.GUI)
                    if obj.GUI
                        disp(['File ''' obj.FileName ''' exists and will be overwritten!'])
                        a = input('Continue (y/n)? ','s');
                    else
                        a = 'y';
                    end
                end
                if strncmpi(a,'n',1)
                    obj.FileID = -1;
                    return
                end
            end
            
            obj.FileID = fopen(obj.FileName, obj.Permission, obj.ByteOrder);
            
            if (obj.FileID == -1)
                error(['@File: Unable to open ''' obj.FileName ''' for ''' obj.Permission '''']);
            end
        end
        
        function obj = fclose(obj)
            %'fclose()' with error checking
            if(obj.Debug), disp('in File.close'); end
            
            if obj.FileID>0 && ~isempty(obj.FileName)
                if(obj.Debug), disp('   closing'); end
                try
                    fclose(obj.FileID);
                    obj.FileID=-1;
                catch
                    obj.FileID=-1;
                end
            end
        end
        
        function obj = freopen(obj)
            %'freopen()' with error checking
            if(obj.Debug), disp('in File.reopen'); end
            obj.fclose();
            obj.fopen();
            
        end
        
        function v = fsize(obj)
            %File size in bytes
            fpos = obj.ftell();    %get current position
            obj.fseek(0,'eof');   %seek to end of file
            v = obj.ftell();       %get current position == filesize in bytes
            obj.fseek(fpos,'bof');%go back to where we began
        end
        function v = ftell(obj)
            %'ftell()' with error checking
            v = ftell(obj.FileID);
            if v <0
                error(['@File: ' ferror(obj.FileID)])
            end
        end
        function fseek(obj,nbytes,from)
            %'fseek()' with error checking
            if fseek(obj.FileID,nbytes,from) <0
                error(['@File: ' ferror(obj.FileID)])
            end
        end
        function v = fread(obj,c,d)
            %'fread()' with error checking and additional data types
            %
            %function v = fread(obj,c,d)
            %
            % Where:
            %   c = number of items to read
            %   d = datatype, eg. 'int8', 'uint8',....,'single','double'
            % Additional data types:
            %   'int24', 'uint24', 'ibm32', 'ieee32', 'ieee64'
            switch d
                case 'int24'
                    c = c*3;
                    p = 'uint8';
                case 'uint24'
                    c = c*3;
                    p = 'uint8';
                case 'ibm32'
                    p = 'uint32';
                case 'ieee32'
                    p = 'single';
                case 'ieee64'
                    p = 'double';
                otherwise
                    p=d;
            end
            
            [v,a] = fread(obj.FileID,c,[p '=>' p]);
            if ~isequal(a,c)
                error(['@File: ' ferror(obj.FileID)])
            end
            
            switch d
                case 'int24'
                    v = int32(int24_2num(v,obj.ByteOrder));
                case 'uint24'
                    v = uint32(uint24_2num(v,obj.ByteOrder));
                case 'ibm32'
                    v = ibm2num(v);                    
            end
        end
        
        function fwrite(obj,v,d)
            %'fwrite()' with error checking and addition datatypes
            %
            %function v = fwrite(obj,v,d)
            %
            % Where:
            %   v = scalar or vector to fwrite
            %   d = datatype, eg. 'int8', 'uint8',....,'single','double'
            % Additional data types:
            %   'int24', 'uint24', 'ibm32', 'ieee32', 'ieee64'            
            [m,n]=size(v);
            c = m*n;
            
            switch d
                case 'int24'
                    c = c*3;
                    p = 'uint8';
                    v = num2int24(v,obj.ByteOrder);
                case 'uint24'
                    c = c*3;
                    p = 'uint8';
                    v = num2uint24(v,obj.ByteOrder);
                case 'ibm32'
                    p = 'uint32';
                    v = num2ibm(v);
                case 'ieee32'
                    p = 'single';
                case 'ieee64'
                    p = 'double';
                otherwise
                    p=d;
            end

            a = fwrite(obj.FileID,v,p);
            if ~isequal(c,a)
                error(['@File: ' ferror(obj.FileID)])
            end
        end
        
        function fprintf(obj,varargin)
           %'fprintf()'
           
           %count = fprintf(obj.FileID,varargin{:});
           fprintf(obj.FileID,varargin{:});
           %error checking??            
        end        
        
        %Listeners
        function obj = listenFileName(obj, varargin)
            %FileName has changed, freopen file for file operations
            if ~isempty(obj.FileName)
               obj.freopen(); 
            end
        end
        function obj = listenByteOrder(obj, varargin)
            %ByteOrder has changed, freopen file for file operations
            
            %disp('In @File/listenByteOrder')
            if ~isempty(obj.FileName)
               obj.freopen(); 
            end        
        end
        function obj = listenPermission(obj, varargin)
            %Permission has changed, freopen file for file operations
            if ~isempty(obj.FileName)
               obj.freopen(); 
            end        
        end
    end %end static methods
    
    methods (Access=private)
        %Destructor
        function delete(obj)
            obj.fclose();
        end
    end %end methods (Access=private)
    
end % end classdef
