function write(obj,th,td)
% function write(obj,th,td)
%
% Writes td and th to disk where:
%  th = trace header struct such as returned by Trace.new
%       or by Trace.read
%  td = numeric matrix with samples in the rows and traces in the
%       columns
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
%

%check input arguments
%narginchk(3,3)
nargin

if ~isnumeric(td)
    error('@Trace/write: Trace data must be numeric')
end

if ~isstruct(th)
    error('@Trace/write: Trace header must be a struct')
end


%% Make certain file headers have already been written to disk
fs = obj.fsize();

if isequal(fs,0)
    error('@Trace/write: File is too small to contain any SEG-Y file headers');
elseif fs <obj.OFFSET
    error(['@Trace/write: File headers are less than ' num2str(obj.FileInfo.TraceOneOffset) ' bytes']);
end

%% Make certain th and td have the same number of traces
%get some basic info about the td array
[ndatsamp, ndattrc] = size(td);
nthdrtrc = size(th.(obj.byte2word(0)),2);

if ~isequal(ndattrc,nthdrtrc)
    error('@Trace/write: Trace data and trace headers contain a different number of traces');
end

%% Apply coordinate scalars
th = obj.applyCoordinateScalars(th); %apply scalars

%% Make certain number of data samples per trace in trace headers matches binary header
nbhdrsamp = obj.FileInfo.SamplesPerTrace;

%check number of data samples and truncate or pad if necessary to match
%binary header
if ndatsamp>nbhdrsamp
    if isempty(obj.GUI) || isa(obj.GUI,'handle')
        a = mm_yesnodlg(...
            ['@Trace/write: Trace data has more samples per trace than expected. '...
            'Truncate trace data?'], ...
            'Warning!','Yes',obj.GUI);
    elseif obj.GUI
        disp('@Trace/write: Trace data has more samples per trace than expected')
        a = input('Truncate trace data (y/n)? ','s');
    else
        a='y';       
    end
    if strncmpi(a,'y',1)
        td = td(1:nbhdrsamp,:);
    else
        if isempty(obj.GUI) || isa(obj.GUI,'handle')
            mm_warndlg('@Trace/write: No traces have been written to disk', ...
                'Warning!', obj.GUI);
        elseif obj.GUI
            disp('@Trace/write: No traces have been written to disk')            
        end
        return
    end
    
elseif ndatsamp<nbhdrsamp
    if isempty(obj.GUI) || isa(obj.GUI,'handle')
        a = mm_yesnodlg(...
            ['@Trace/write: Trace data has fewer samples per trace than expected. '...
            'Pad trace data with zeros?'], ...
            'Warning!','Yes',obj.GUI);
    elseif obj.GUI
        disp('@Trace/write: Trace data has fewer samples per trace than expected')
        a = input('Pad trace data with zeros (y/n)? ','s');
    else
        a='y';        
    end
        
    if strncmpi(a,'y',1)
        w = whos('td');
        td = [td; zeros(nbhdrsamp-ndatsamp,ndattrc,w.class)];
    else
        if isempty(obj.GUI) || isa(obj.GUI,'handle')
            mm_warndlg('@Trace/write: No traces have been written to disk', ...
                'Warning!', obj.GUI);
        elseif obj.GUI
            disp('@Trace/write: No traces have been written to disk')            
        end
        return
    end
end

%overwrite number of samples in this trace if trace headers matches binary header
th.(obj.byte2word(114))(1:end) = nbhdrsamp; %stored in byte 115; SEG-Y standard

%% Make certain sample interval matches this trace in trace headers 
thdrdt = th.(obj.byte2word(116))(:); %stored in byte 117; SEG-Y standard
bhdrdt = obj.FileInfo.SampleInterval; %read from bin hdr already written to disk

if ~isequal(sum(thdrdt/bhdrdt),ndattrc)
    if isempty(obj.GUI) || isa(obj.GUI,'handle')
        a = mm_yesnodlg(...
            ['@Trace/write: Trace header sample interval differs from binary file header. '...
            'Update trace headers?'], ...
            'Warning!','Yes',obj.GUI);
    elseif obj.GUI
        disp('@Trace/write: Trace header sample interval differs from binary file header')
        a = input('Update trace headers [y/n]? ','s');
    else
        a='y';
    end
    
    if strncmpi(a,'y',1)
        %overwrite sample interval in this trace in trace headers
        fieldname = obj.byte2word(116);
        th.(fieldname)(1:end) = bhdrdt;       
    else
        if isempty(obj.GUI) || isa(obj.GUI,'handle')
            mm_warndlg('@Trace/write: No traces have been written to disk', ...
                'Warning!', obj.GUI);
        elseif obj.GUI
            disp('@Trace/write: No traces have been written to disk')            
        end
        return
    end
end

%% renumber trace sequence number within SEG-Y file trace header
fieldname = obj.byte2word(5); %stored in byte 5; SEG-Y standard
th.(fieldname) = 1:ndattrc +obj.FileInfo.TracesInFile;

%% Final check of trace headers and trace data to write
obj.check(th);
obj.check(td);

% disp(['Trace.write datatype: ' datatype]);
% disp(['Trace.write Trace 1, sample1: ' num2str(td(1,1))]);
    
%% Write traces to disk
%We are appending traces to end of file
obj.fseek(0,'eof');

%Write trace headers and trace data
% 
if isempty(obj.GUI) || isa(obj.GUI, 'handle')
    h = waitbar(0,['Writing ' num2str(ndattrc) ' Trace(s): '],...
        'Name','SEG-Y output',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    mm_adjust(h);
    setappdata(h,'canceling',0);
end

for jj = 1:ndattrc
       
    if isempty(obj.GUI) || isa(obj.GUI, 'handle')
        %check for cancel button press
        if getappdata(h,'canceling')
            delete(h);
            break
        end
        %update progress bar
        waitbar(jj/ndattrc,h);
    elseif obj.GUI
        %text progress bar
        %text progress bar
        textbar(['Writing ' num2str(ndattrc) ' Trace(s): '],30,jj,ndattrc);
    end
            
    %write trace hdr; hdrdef is in obj
    for ii = 1:length(obj.HdrDef)
        v = th.(obj.HdrDef{ii})(jj);
        
        datafmt = obj.HdrDef{ii,2};       
        obj.fwrite(v,datafmt);
        
    end %end for loop for header word writing

	obj.fwrite(td(:,jj),obj.FormatCodeType);

end % end for jj

if isempty(obj.GUI) || isa(obj.GUI, 'handle')
    %close waitbar
    delete(h);
end
    
end %end function test_write