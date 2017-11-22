function bd = newDefinition(segyrev)
%
% function bd = newDefinition(segyrev)
%
% newDefinition creates a cell array 'bd' where:
%
%  The binary file header definition is N rows by 4 columns where:
%    column1: Short header name (no spaces; can be used as fieldname in struct)
%    column2: Long header name (can have spaces)
%    column3: Byte location within binary header (for use with fseek)
%    column4: Format (eg. uint8, uint16, uint32, int8, int16, int32, 
%             ieee32, ibm32, single, double)
%
% Authors: Kevin Hall, 2016, 2017
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

if nargin < 1 || isempty(segyrev)
    segyrev=1;
end

bd(1,:)={'JobID','int32',0,'Job identification number'};
bd(2,:)={'ReelNum','int32',4,'Line number'};
bd(3,:)={'LineNum','int32',8,'Reel number'};
bd(4,:)={'DataTrcPerEns','int16',12,'*Number of data traces per ensemble'};
bd(5,:)={'AuxTrcPerEns','int16',14,'*Number of auxiliary traces per ensemble'};
bd(6,:)={'SampleRate','uint16',16,'*Sample interval in microseconds '};
bd(7,:)={'OrigSampleRate','uint16',18,'Sample interval of field data in microseconds '};
bd(8,:)={'SampPerTrc','uint16',20,'*Number of samples per data trace'};
bd(9,:)={'OrigSampPerTrc','uint16',22,'Number of samples recorded in field data per data trace'};
bd(10,:)={'FormatCode','uint16',24,'*Data sample format code'};
bd(11,:)={'EnsembleFold','int16',26,'*Ensemble fold '};
bd(12,:)={'TrcSortCode','int16',28,'*Trace sorting code'};
bd(13,:)={'VertSumCode','int16',30,'Vertical sum code'};
bd(14,:)={'SweepStartFreq','int16',32,'Sweep frequency at start (Hz)'};
bd(15,:)={'SweepEndFreq','int16',34,'Sweep frequency at end (Hz)'};
bd(16,:)={'SweepLength','int16',36,'Sweep length (ms)'};
bd(17,:)={'SweepType','int16',38,'Sweep type code'};
bd(18,:)={'SweepTrcNum','int16',40,'Trace number of sweep channel'};
bd(19,:)={'SweepStartTaper','int16',42,'Sweep trace taper length (ms)'};
bd(20,:)={'SweepEndTaper','int16',44,'Sweep trace taper length (ms)'};
bd(21,:)={'SweepTaperType','int16',46,'Taper type'};
bd(22,:)={'CorrelatedData','int16',48,'Correlated data traces'};
bd(23,:)={'BinaryGainRecvd','int16',50,'Binary gain recovered'};
bd(24,:)={'AmpRecMethod','int16',52,'Amplitude recovery method'};
bd(25,:)={'MeasurementSys','int16',54,'*Measurement system'};
bd(26,:)={'ImpSigPolarity','int16',56,'Impulse signal polarity'};
bd(27,:)={'VibeSigPolarity','int16',58,'Vibratory polarity code'};

if segyrev == 0
    bd = [bd; mkunassigned(1,60,60)];
    bd(88,:)={'SegyRevNumMaj','uint8',300,'*SEGY Format Revision Major Number'};
    bd(89,:)={'SegyRevNumMin','uint8',301,'*SEGY Format Revision Minor Number'};
    bd(90,:)={'Unassigned61','uint16',302,'Unassigned 61'};
    bd = [bd; mkunassigned(62,85,304)];    
elseif segyrev == 1
    bd = [bd; mkunassigned(1,60,60)];
    bd(88,:)={'SegyRevNumMaj','uint8',300,'*SEGY Format Revision Major Number'};
    bd(89,:)={'SegyRevNumMin','uint8',301,'*SEGY Format Revision Minor Number'};    
    bd(90,:)={'FixedLenTraces','int16',302,'*Fixed length trace flag'};
    bd(91,:)={'NumExtTextHdrs','int16',304,'*Number Extended Textual File Headers'};
    bd(92,:)={'Unassigned61','uint16',306,'Unassigned 61'};
    bd = [bd; mkunassigned(62,84,308)];    
elseif segyrev == 2
    bd(28,:)={'ExtDataTrcPerEns','uint32',60,'Extended Number of data traces per ensemble'};
    bd(29,:)={'ExtAuxTrcPerEns','uint32',64,'Extended Number of auxiliary traces per ensemble'};
    bd(30,:)={'ExtSampPerTrc','uint32',68,'Extended Number of samples per data trace'};
    bd(31,:)={'ExtSampleRate','ieee64',72,'Extended Sample interval in microseconds'};
    bd(32,:)={'ExtOrigSampleRate','ieee64',80,'Extended Sample interval of field data in microseconds '};
    bd(33,:)={'ExtOrigSampPerTrc','uint32',88,'Extended Number of samples recorded in field data per data trace'};
    bd(34,:)={'ExtEnsembleFold','uint32',92,'Extended Ensemble fold'};
    bd(35,:)={'IntegerConstant','uint32',96,'Integer Constant == 16909060'};   
    bd = [bd; mkunassigned(1,50,100)];
    bd(86,:)={'SegyRevNumMaj','uint8',300,'*SEGY Format Revision Major Number'};
    bd(87,:)={'SegyRevNumMin','uint8',301,'*SEGY Format Revision Minor Number'};
    bd(88,:)={'FixedLenTraces','int16',302,'*Fixed length trace flag'};
    bd(89,:)={'NumExtTextHdrs','int16',304,'*Number Extended Textual File Headers'};
    bd(90,:)={'NumExtTrcHdrs','uint32',306,'*Number Extended Trace Headers'};
    bd(91,:)={'TimeBasisCode','uint16',310,'Time Basis Code'};    
    bd(92,:)={'NumTrcInFile','uint64',312,'Number of traces in the file or stream'};
    bd(93,:)={'TraceOneOffset','uint64',320,'Byte offset of first trace from start of file or stream'};
    bd(94,:)={'NumDataTrailers','int32',328,'Number of binary or textual data trailers'};
    bd = [bd; mkunassigned(51,67,332)];
else
    error('@BinaryHeader/newDefinition: Unknown SEG-Y revision number');
end

end %end function

function ca = mkunassigned(widx1,widx2,sbyte)
    nrow = widx2-widx1+1;
    ca = cell(nrow,4);
    cidx=1;
    for ii=widx1:widx2
        ca(cidx,:)={sprintf('Unassigned%02d',ii),'uint32',sbyte,...
            sprintf('Unassigned %d',ii)};
        cidx=cidx+1;
        sbyte=sbyte+4;
    end
end


