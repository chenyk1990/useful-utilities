function bh = new ( obj, nsamp, sampint, segfmt, segrev, ndattrc, nauxtrc, ensfold, srtcode, meassys, nexthdrs )
%
% function bh = new ( nsamp, sampint, ...    %required inputs
%   segfmt, segrev, ndattrc, nauxtrc, ...   %optional inputs
%   ensfold, srtcode, nexthdrs )            %optional inputs
%
% nsamp    = number of samples per trace (mandatory)
% sampint  = sample interval (s) (mandatory)
% segfmt   = SEG-Y format code; in range 1:5,8 (4 is obsolete!) (mandatory)
%            1=4-byte IBM float, 2=4-byte integer, 3=2-byte integer, 4 is
%            not supported by this code, 5=4-byte IEEE float (default), 8=1-byte
%            integer
% segrev   = SEG-Y revision (mandatory) 0='rev0', 1='rev1' (default),
%            2='rev2'
%            NOTE: rev1 will be represented by the number 256 in the struct
% ndattrc  = number of data traces per ensemble (mandatory for pre-stack)
% nauxtrc  = number of aux traces per ensemble (mandatory for pre-stack)
% ensfold  = ensemble fold (eg. CMP fold) (recommended) 1=default
% srtcode  = trace sort code; in range -1:9 (recommended)
%            -1=other, 0=unknown (default), 1=no sort, 2=common depth point,
%            3=single fold continuous, 4=horizontal stack, 5=common source 
%            point, 6=common receiver point, 7=common offset point,
%            8=common mid point, 9=common converstion point
% meassys  = measurement sytem (recommended) 1=meters (default), 2=feet
% nexthdrs = number of extended textual headers; in range -1:N (mandatory 
%            for SEG-Y Rev 1) 0=default
%
% Creates a new SEG-Y binary file header struct
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

narginchk(3, 11)

%% Check input values; set defaults if necessary
if nargin < 4 || isempty(segfmt)
    segfmt = 5; %mandatory; 1:5,8
end
if nargin < 5 || isempty(segrev)
    segrev = 1; %used this script only
end
if nargin < 6 || isempty(ndattrc)
    ndattrc = 0; %mandatory for pre-stack data
end
if nargin < 7 || isempty(nauxtrc)
    nauxtrc = 0; %mandatory for pre-stack data
end
if nargin < 8 || isempty(ensfold)
    ensfold = 1; %recommended
end
if nargin < 9 || isempty(srtcode)
    srtcode = 0; %recommended; -1:9
end
if nargin < 10 || isempty(meassys)
    meassys = 1; %recommended; 1:2
end
if nargin < 11 || isempty(nexthdrs)
    nexthdrs = 0; %mandatory for SEG-Y rev1
end


%% Sanity checks
%sample interval
if sampint <0
    error('@BinaryHeader/new: Sample interval must be positive');
end
%number samples per trace
if nsamp <0
    error('@BinaryHeader/new: Number of samples per trace must be positive');
end
%segy revision number
if segrev==0
    %do nothing?
elseif segrev==1
    segrev=1; %trust us; this is represents a 1.0 on disk
elseif segrev==2
    segrev=2; %trust us; this represents a 2.0 on disk
else
    error(['@BinaryHeader/new: Unknown SEG-Y revision number: ' num2str(segrev)])
end 
%format code
if segfmt<1 || segfmt==4 || segfmt==7 || (segfmt>12 && segfmt<16) || segfmt> 16
    error(['@BinaryHeader/new: Format code ' num2str(segfmt) ' is not supported by this code']);
elseif segrev==0 && segfmt >4
    error(['@BinaryHeader/new: Format code ' num2str(segfmt) ' is not defined in SEG-Y rev 0']);
elseif segrev==1 && (segfmt>5 && segfmt<8 || segfmt >8)
    error(['@BinaryHeader/new: Format code ' num2str(segfmt) ' is not defined in SEG-Y rev 1']);    
end

% num data traces
if ndattrc <0
    error('@BinaryHeader/new: Number of data traces must be positive');
end
%num aux traces
if nauxtrc <0
    error('@BinaryHeader/new: Number of auxilliary traces must be positive');
end
%ensemble fold
if ensfold <0
    error('@BinaryHeader/new: Ensemble fold must be positive');
end
%trace sort code
if segrev==0 && (srtcode < 0 || srtcode > 4)
    error(['@BinaryHeader/new: Trace sort code ' num2str(strcode) ' is not supported in SEG-Y rev 0']);
elseif (srtcode < -1 || srtcode > 9)
    error(['@BinaryHeader/new: Trace sort code ' num2str(strcode) ' is not supported in SEG-Y rev 1']);
end
%measurement system
if meassys <1 || meassys >2
    error(['@BinaryHeader/new: Unkown measurement system ' num2str(meassys)]);
end
%fixed trace length flag
if segrev==0
    fixlenflag=0;
else
    fixlenflag=1; %refusing to deal with variable trace lengths!
end
%number of extended textual file headers
if segrev==0
    nexthdrs=0; %extended textual file headers are not a thing in rev. 0
elseif nexthdrs < -1
    error('@BinaryHeader/new: Number of Extended Textual File Headers cannot be less than -1');
end
                        
%% Create binary file header struct from header definition
for ii = 1:length(obj.HdrDef)
    %initialize struct
    switch obj.HdrDef{ii,2}
        case 'uint8'
            bh.(obj.HdrDef{ii,1}) = uint8(0);
        case 'int8'
            bh.(obj.HdrDef{ii,1}) = int8(0);            
        case 'uint16'
            bh.(obj.HdrDef{ii,1}) = uint16(0);            
        case 'int16'
            bh.(obj.HdrDef{ii,1}) = int16(0);
         case 'uint24'
            bh.(obj.HdrDef{ii,1}) = uint24(0);            
        case 'int24'
            bh.(obj.HdrDef{ii,1}) = int24(0);           
        case 'uint32'
            bh.(obj.HdrDef{ii,1}) = uint32(0);             
        case 'int32'           
            bh.(obj.HdrDef{ii,1}) = int32(0); 
         case 'uint64'
            bh.(obj.HdrDef{ii,1}) = uint64(0);             
        case 'int64'
            bh.(obj.HdrDef{ii,1}) = int64(0);
         case 'single'
            bh.(obj.HdrDef{ii,1}) = single(0.0);           
        case 'ieee32'
            bh.(obj.HdrDef{ii,1}) = single(0.0);
        case 'ibm32'
            bh.(obj.HdrDef{ii,1}) = single(0.0);
         case 'double'
            bh.(obj.HdrDef{ii,1}) = double(0.0);            
        case 'ieee64'
            bh.(obj.HdrDef{ii,1}) = double(0.0);            
        otherwise
            error('@BinaryHeader/new: Unknown data type')
    end
    %set struct values
    switch obj.HdrDef{ii,3} %start byte in header (= SEG-Y standard byte 
                            %location-3201)
        case 12 %data traces per ensemble (mandatory for pre-stack)
            bh.(obj.HdrDef{ii,1}) = ndattrc;
        case 14 %aux traces per ensemble (mandatory for pre-stack)
            bh.(obj.HdrDef{ii,1}) = nauxtrc;
        case 16 %sample rate in microseconds (mandatory)
            bh.(obj.HdrDef{ii,1}) = sampint*1e6;
        case 20 %number of samples per trace (mandatory)
            bh.(obj.HdrDef{ii,1}) = nsamp;
        case 24 %format code (mandatory)
            bh.(obj.HdrDef{ii,1}) = segfmt;
        case 26 %ensemble fold (recommended)
            bh.(obj.HdrDef{ii,1}) = ensfold;
        case 28 %trace sorting code (recommended)
            bh.(obj.HdrDef{ii,1}) = srtcode;
        case 54 %measurement system (recommended)
            bh.(obj.HdrDef{ii,1}) = meassys;
        case 300 %SEG-Y revision number (mandatory)
            bh.(obj.HdrDef{ii,1}) = segrev;
        case 302 %Fixed length trace flag (mandatory)
            bh.(obj.HdrDef{ii,1}) = fixlenflag;
        case 304 %Number of extended textual file headers (mandatory)
            bh.(obj.HdrDef{ii,1}) = nexthdrs;
    end
end

end