function obj = guessFormatCode(obj)
%
% function fc = guessFormatCode(obj)
%
% Algorithm: If SEG-Y FormatCode = 1 read the first trace from the file as
% uint32, then use ibm2num and num2ibm to convert to single and back to 
% uint32. If the result does not match the input then trace data may 
% actually be IEEE floating point.
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

%Format code 1 is ambiguous, could be IBM or IEEE floats

if isequal(obj.FileInfo.FormatCode,1)
    
    %fseek to start of first trace
    obj.fseek(obj.FileInfo.TraceOneOffset+obj.SIZE, 'bof');
    %fread first trace data as uint32
    trcin     = obj.fread(obj.FileInfo.SamplesPerTrace, 'uint32');
    %convert to single assuming data is IBM floating point
    trcdbl    = ibm2num(trcin);
    %convert back to uint32
    trcout    = num2ibm(trcdbl);
    
%     plot(trcin-trcout)
%     sum(trcin-trcout)

    %Test
    if sum(trcin-trcout) %trcin and trcout are identical
        if isempty(obj.GUI) || isa(obj.GUI,'handle')
            mm_warndlg(['@Trace/guessFormatCode: Binary Header Format Code is 1 (4-byte IBM), '...
                'but trace data appear to be Format Code 5 (4-byte IEEE)'],...
                'Warning!',obj.GUI);
        elseif obj.GUI
            warning(['@Trace/guessFormatCode: Binary Header Format Code is 1 (4-byte IBM), '...
                      'but trace data appear to be Format Code 5 (4-byte IEEE)']);
        end
        %Update FormatCode to 5
        obj.FileInfo.FormatCode = 5;
    end    
end

end %end function