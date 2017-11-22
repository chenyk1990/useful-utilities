function checkArray(obj, td)
%
% function checkArray(obj, td)
%
% CHECKARRAY checks the trace data values in td to see if they can be
% stored in a SEG-Y file using the current format code; ie. MIN < td < MAX
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

switch(obj.FileInfo.FormatCode)
    case 1
%         datfmt  = 'ibm32';
%         if ~check_fp_range(td,'single')
%         %It turns out that Vista and ProMAX do not correctly read IBM
%         %floats bigger or smaller than IEEE max and min, so DO NOT check_fp_range(td,'ibm32')
%             error('@Trace/checkArray: Trace data contains values that cannot be stored as ''ibm32''')
%         end
    case 2
%         datfmt  = 'int32';
        if ~check_int_range(td,'int32')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''int32''')
        end
    case 3
%         datfmt  = 'int16';
        if ~check_int_range(td,'int16')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''int16''')
        end
    case 4        
        error('@Trace/checkArray: SEG-Y Format code 4 (4-byte fixed point with gain) is not supported');
    case 5
%         datfmt  = 'single';
        if ~check_fp_range(td,'ieee32')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''ieee32''')
        end
    case 6
%         datfmt  = 'double';
        if ~check_fp_range(td,'ieee64')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''ieee64''')
        end
    case 7
%         datfmt  = 'int24';
        if ~check_int_range(td,'int24')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''int24''')
        end                
    case 8
%         datfmt  = 'int8';
        if ~check_int_range(td,'int8')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''int8''')
        end
    case 9
%         datfmt  = 'int64';
        if ~check_int_range(td,'int64')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''int8''')
        end
    case 10
%         datfmt  = 'uint32';
        if ~check_int_range(td,'uint32')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''uint32''')
        end
    case 11
%         datfmt  = 'uint16';
        if ~check_int_range(td,'uint16')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''uint16''')
        end        
    case 12
%         datfmt  = 'uint64';
        if ~check_int_range(td,'uint64')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''uint64''')
        end        
    case 15
%         datfmt  = 'uint24';
        if ~check_int_range(td,'uint24')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''uint24''')
        end
    case 16
%         datfmt  = 'uint8';
        if ~check_int_range(td,'uint8')
            error('@Trace/checkArray: Trace data contains values that cannot be stored as ''uint8''')
        end        
    otherwise
        error(['@Trace/checkArray: Unknown SEG-Y Format code: ' num2str(obj.FormatCode)])
end

end %end function