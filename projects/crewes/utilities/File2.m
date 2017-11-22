classdef File2 < handle
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
%   fread(sf.FileID, 3200, '*uchar') %system fread with no error checking
%    OR
%   File.read(sf.FileID, 3200, '*uchar') % fread with error checking
%
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
    FileName;  % filename corresponding to fid
    ByteOrder; % see 'help fopen'
    GUI;       % GUI flag 0 no prompts,1 text prompts,h gui prompts
    FileID;    
end

properties (Hidden)
   Permission ='r'; % see 'help fopen'    
   FirstOpen=true; %true or false
end

%events
events
    FileNameChanged;
    ByteOrderChanged;
    PermissionChanged;
end

%methods
methods
    
%Constructor
function obj = File2(filename, permission, byteorder, gui)
    % Constructor for File object
    % check number of input arguments
    narginchk(1,4)
        
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
    
    %set object values
    if ischar(filename)
       obj.FileName = filename;
    elseif isnumeric(filename)
       obj.FileID = filename;
       obj.FileName = [];    
    end
    
    obj.Permission = permission;
    obj.ByteOrder = byteorder;
    
    obj.GUI = gui;
    
    % Add listeners
    addlistener(obj,'FileNameChanged',@obj.listenFileName);
    addlistener(obj,'ByteOrderChanged',@obj.listenByteOrder);
    addlistener(obj,'PermissionChanged',@obj.listenPermission);    
    
    % Open file
    if ischar(filename)
        obj = obj.fopen();
    end
end %end constructor

%Set functions
function set.FileName(obj, v)
    if (ischar(v))
        obj.FileName = v;
        notify(obj,'FileNameChanged');
    elseif isempty(v)
        obj.FileName = [];        
    else
        error('@File: FileName must be a character string');
    end
end

function set.ByteOrder(obj, v) 
    if (ischar(v))
        v = lower(v);
        switch(v)
            case 'native' %- local machine format - the default
                obj.ByteOrder='n';
                notify(obj,'ByteOrderChanged');
            case 'n'
                obj.ByteOrder=v;
                notify(obj,'ByteOrderChanged');
            case'ieee-le' %- IEEE floating point with big-endian
                obj.ByteOrder='l';
                notify(obj,'ByteOrderChanged');
            case 'l'
                obj.ByteOrder=v;
                notify(obj,'ByteOrderChanged');               
            case'ieee-be' %- IEEE floating point with big-endian
                obj.ByteOrder='b';
                notify(obj,'ByteOrderChanged');
            case 'b'
                obj.ByteOrder=v;
                notify(obj,'ByteOrderChanged');
            case'ieee-le.l64' %- IEEE floating point with little-endian
                obj.ByteOrder='a';
                notify(obj,'ByteOrderChanged');
            case 'a'
                obj.ByteOrder=v;
                notify(obj,'ByteOrderChanged');                 
            case'ieee-be.l64' %- IEEE floating point with big-endian byte
                obj.ByteOrder='s';
                notify(obj,'ByteOrderChanged');
            case 's'
                obj.ByteOrder=v;
                notify(obj,'ByteOrderChanged');
            otherwise
                error(['@File: Unknown ByteOrder ''' v ''' see help fopen'])
        end
    elseif (isempty(v))
        obj.ByteOrder='n';
        notify(obj,'ByteOrderChanged');
    else
       error('@File: ByteOrder must be a character string') 
    end
end

function set.Permission(obj, v)
    if (ischar(v))
        if strcmp(v,'r') || strcmpi(v,'w') || strcmpi(v,'a') ||...
                        strcmp(v,'r+') || strcmp(v,'w+') || strcmp(v,'a+')
            obj.Permission=v;
            notify(obj,'PermissionChanged');
        else
            error(['@File: Permission ''' v ''' is invalid: see ''help fopen'''])
        end
    elseif (isempty(v))
        obj.Permission = 'r';
    else
       error('@File: Permission must be a character string') 
    end
end

function set.GUI(obj, v)
    if isnumeric(v)
        obj.GUI = logical(v);
    elseif islogical(v)
        obj.GUI = v;
    elseif isempty(v) || isa(v,'handle')
        obj.GUI = v;
    else
        error('@File: GUI must be numeric, empty, or a figure handle')
    end
end

%Open and close file
function obj = fopen(obj)
%     disp('in File.open')
    if obj.FirstOpen
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
        obj.FirstOpen = false;
    end
    
%    obj.FileName,obj.Permission, obj.ByteOrder
    
    obj.FileID = fopen(obj.FileName, obj.Permission, obj.ByteOrder);
    
    if (obj.FileID == -1)
        error(['@File: Unable to open ''' obj.FileName ''' for ''' obj.Permission '''']);
    end
end

function obj = fclose(obj)
%     disp('in File.close')
    try
        obj.FileID = fclose(obj.FileID);
    catch
        obj.FileID = -1;
    end
end

function obj = freopen(obj)
%     disp('in File.reopen')
    obj.fclose();
    obj.fopen();
end

%Listeners
function obj = listenFileName(obj, varargin)
%     disp('in File.listenFileName')
    obj.freopen();
end

function obj = listenByteOrder(obj, varargin)
%      disp('in File.listenByteOrder')
    obj.freopen();
end

function obj = listenPermission(obj, varargin)
%     disp('in File.listenPermission')
    obj.freopen();
end

end % end methods

methods
    function v = fsize(obj)
        fpos = obj.ftell();    %get current position      
        obj.fseek(0,'eof');   %seek to end of file
        v = obj.ftell();       %get current position == filesize in bytes
        obj.fseek(fpos,'bof');%go back to where we began
    end
    function v = ftell(obj)
        v = ftell(obj.FileID);
        if v <0
            error(['@File: ' ferror(obj.FileID)])
        end       
    end
    function fseek(obj,nbytes,from)
        if fseek(obj.FileID,nbytes,from) <0
            error(['@File: ' ferror(obj.FileID)])
        end
    end
    function v = fread(obj,c,d)
        %function v = read(obj,c,d)
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
        end
                
        [v,a] = fread(obj.FileID,c,[p '=>' p]);
        if ~isequal(a,c)
            error(['@File: ' ferror(obj.FileID)])            
        end
        
        switch d
            case 'int24'
                v = int24_2num(v,obj.ByteOrder);
            case 'uint24'
                v = uint24_2num(v,obj.ByteOrder);
            case 'ibm32'
                v = ibm2num(v);
        end                      
    end
    
    function fwrite(obj,v,p)
        [m,n]=size(v);
        c = m*n;
        switch d
            case 'int24'
                c = c*3;
                p = 'uint8';
                v = num2int24(v);
            case 'uint24'
                c = c*3;
                p = 'uint8';
                v = num2uint24(v);
            case 'ibm32'
                p = 'uint32';
                v = num2ibm(v);
            case 'ieee32'
                p = 'single';
            case 'ieee64'
                p = 'double';
        end       
        
        a = fwrite(obj.FileID,v,p);
        if ~isequal(a,c)
            error(['@File: ' ferror(obj.FileID)])
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
