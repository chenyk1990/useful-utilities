function finfo = readFileInfo(obj)
%
% function finfo = readFileInfo(obj)
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
% END TERMS OF USE LICENSE()

finfo.SegyRevision       = 1;
finfo.FormatCode         = 5; % 4-byte IEEE
finfo.SamplesPerTrace    = 0;
finfo.TracesPerRec       = 0;
finfo.SampleInterval     = 1000; % 1 ms
finfo.FixedTrcLength     = 0; %fixed trace length is true
finfo.NumExtTxtHeaders   = 0;
finfo.ExtTracesPerRec    = 0;
finfo.ExtSamplesPerTrace = 0;
finfo.ExtSampleInterval  = 0;
finfo.IntegerConstant    = 0;
finfo.NumExtTrcHeaders   = 0;
finfo.TraceOneOffset     = 0;

if obj.FileID > 0 && obj.fsize >= obj.SIZE
    finfo.SegyRevision    = obj.readSegyRevision();
    finfo.FormatCode      = obj.readFormatCode();
    finfo.SamplesPerTrace = obj.readSamplesPerTrace();
    finfo.TracesPerRec    = obj.readTracesPerRec();
    finfo.SampleInterval  = obj.readSampleInterval();

    if finfo.SegyRevision >0
        finfo.FixedTrcLength    = obj.readFixedTrcLength();
        finfo.NumExtTxtHeaders  = obj.readNumExtHeaders();
    end
    
    if finfo.SegyRevision >1
        finfo.ExtTracesPerRec    = obj.readExtTracesPerRec();
        finfo.ExtSamplesPerTrace = obj.readExtSamplesPerTrace();
        finfo.ExtSampleInterval  = obj.readExtSampleInterval();
        finfo.IntegerConstant    = obj.readIntegerConstant();
        finfo.NumExtTrcHeaders   = obj.readNumExtTrcHeaders();
        finfo.TraceOneOffset     = obj.readTraceOneOffset();
        finfo.NumDataTrailers    = obj.readNumDataTrailers();
    end
end
