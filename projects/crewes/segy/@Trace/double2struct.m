function s = double2struct(d,f,t)
%
% function s = double2struct(d,f,t)
% Where:
%     s = struct such as produced by Trace/read for the trace headers
%     d = all values stored in s as a 2D matrix of doubles
%     f = all fieldnames in s as a cell array
%     t = cell array of datatypes (eg. 'int8', 'uint16')
%
% Authors: Kevin Hall, 2017
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.
%

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

[m,n] = size(d);
d = mat2cell(d,ones(1,m),n);

d = cellfun(@(X,Y) convertdatatype(X,Y),d,t,'UniformOutput',false);


s = cell2struct(d,f);

end


function data = convertdatatype(data,datatype)
    allswell = true;
    switch datatype
        case 'uint8'
            if check_int_range(data,datatype);
                data = uint8(data);
            else
                allswell = false;
            end
        case 'uint16'
            if check_int_range(data,datatype);
                data = uint16(data);
            else
                allswell = false;
            end
        case 'uint32'
            if check_int_range(data,datatype);
                data = uint32(data);
            else
                allswell = false;
            end
        case 'uint64'
            if check_int_range(data,datatype);
                data = uint64(data);
            else
                allswell = false;
            end                
        case 'int8'
            if check_int_range(data,datatype);
                data = int8(data);
            else
                allswell = false;
            end
        case 'int16'
            if check_int_range(data,datatype);
                data = int16(data);
            else
                allswell = false;
            end
        case 'int32'
            if check_int_range(data,datatype);
                data = int32(data);
            else
                allswell = false;
            end
        case 'int64'
            if check_int_range(data,datatype);
                data = int64(data);
            else
                allswell = false;
            end                       
        case 'single'
            if check_fp_range(data,datatype);
                data = single(data);
            else
                allswell = false;
            end             
        case 'double'
            if check_fp_range(data,datatype);
                data = double(data);
            else
                allswell = false;
            end 
        otherwise
            error('@Trace/double2struct: Unknown data type');
    end
    if ~allswell
        error(['@Trace/double2struct: Data contains NaN or Inf or is out of range for storage '
               'as ' datatype]);
    end
end
