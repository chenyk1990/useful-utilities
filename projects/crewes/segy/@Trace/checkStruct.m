function checkStruct(obj, th)
%
%function checkStruct(obj, th)
%
% CHECKSTRUCT compares a binary header struct to the current binary header
% definition in the object for:
% 1. fieldname in the current definition exists in the struct
% 2. value in struct.fieldname are scalar (or vector) and match the number
%    of traces
% 3. value in struct.fieldname can be stored as the datatype in the
%    definition. eg. MIN < struct.fieldname < MAX
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

%make sure th is a struct
if ~isstruct(th)
    error('@Trace/checkStruct: Trace Header to write to disk must be a struct')
end

ndefs   = length(obj.HdrFieldNames);
nfields = 0;
oldntr  = [];

for ii = 1:ndefs
    %check if header definition field exists in struct
    fieldname = obj.HdrFieldNames{ii};
    datatype  = obj.HdrDataTypes.(fieldname);
    if ~isfield(th,fieldname)
        error(['@Trace/checkStruct: Trace header struct must contain fieldname ''' ...
            fieldname ''''])
    else
        nfields = nfields+1;
    end

    %check size of header word
    [nrow,ntr] = size(th.(fieldname));
    if ~isequal(nrow,1)
        error(['@Trace/checkStruct: Trace header value for fieldname ''' ...
            fieldname ...
            ''' must be a row vector'])
    end

    %check number of traces
    if ~isempty(oldntr) && ~isequal(oldntr,ntr)
        error(['@Trace/checkStruct: Trace header value for fieldname ''' ...
            fieldname ...
            ''' contains data for a different number of traces than expected'])       
    else
        oldntr = ntr;
    end
    
    %check to see if struct value is in range for datatype
    if isempty(strfind(datatype,'int'))
        if ~check_fp_range(th.(fieldname),datatype)
            error(['@Trace/checkStruct: Trace header value for fieldname ''' ...
                fieldname ...
                ''' cannot be stored as datatype '''...
                datatype ''''])
        end
    else
        if ~check_int_range(th.(fieldname),datatype)
            error(['@Trace/checkStruct: Trace header value for fieldname ''' ...
                fieldname ...
                ''' cannot be stored as datatype '''...
                datatype ''''])
        end
    end
    
end %end for

end %end function checkStruct
    

