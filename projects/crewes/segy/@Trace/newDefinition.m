function td = newDefinition(segyrev,numexthdr)
%
% function td = newDefinition(segyrev,numexthdr)
%
% newDefinition creates a cell array 'td' where:
%
%  The binary file header definition is N rows by 4 columns where:
%    column1: Short header name (no spaces; can be used as fieldname in struct)
%    column2: Long header name (can have spaces)
%    column3: Byte location within binary header (for use with fseek)
%    column4: Format (eg. uint8, uint16, uint32, int8, int16, int32, 
%             ibm32, ieee32, single, ieee64, double)
%
% Authors: Chad Hogan, 2008
%          Kevin Hall, 2016, 2017
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

if nargin < 2 || isempty(numexthdr)
    numexthdr = 0;
end

td(1,:)={'TrcNumLine','int32',0,'','*Trace sequence number within line'};
td(2,:)={'TrcNumReel','int32',4,'','Trace sequence number within SEGY file'};
td(3,:)={'FieldRecNum','int32',8,'','*Original field record number'};
td(4,:)={'TrcNumFile','int32',12,'','*Trace number within the original field record'};
td(5,:)={'SourcePoint','int32',16,'','Energy source point number'};
td(6,:)={'EnsembleNum','int32',20,'','Ensemble number'};
td(7,:)={'TrcNumCDP','int32',24,'','Trace number within the ensemble'};
td(8,:)={'TrcCode','int16',28,'','*Trace identification code'};
td(9,:)={'VertFold','int16',30,'','Number of vertically summed traces yielding this trace'};
td(10,:)={'HorzFold','int16',32,'','Number of horizontally stacked traces yielding this trace'};
td(11,:)={'DataUse','int16',34,'','Data use'};
td(12,:)={'SrcRecOffset','int32',36,'','Distance from center of the source point to the center of the receiver group'};
td(13,:)={'GroupElev','int32',40,'ElevScalar','Receiver group elevation'};
td(14,:)={'SrcElev','int32',44,'ElevScalar','Surface elevation at source'};
td(15,:)={'SrcDepth','int32',48,'ElevScalar','Source depth below surface'};
td(16,:)={'GroupDatumElev','int32',52,'ElevScalar','Datum elevation at receiver group'};
td(17,:)={'SrcDatumElev','int32',56,'ElevScalar','Datum elevation at source'};
if segyrev < 2
    td(18,:)={'SrcWaterDepth','int32',60,'ElevScalar','Water depth at source'};
    td(19,:)={'GroupWaterDepth','int32',64,'ElevScalar','Water depth at group'};
else
    td(18,:)={'SrcWaterColHeight','int32',60,'ElevScalar','Water column height at source'};
    td(19,:)={'GroupWaterColHeight','int32',64,'ElevScalar','Water column height at group'};
end
td(20,:)={'ElevScalar','int16',68,'','Scalar to be applied to all elevations and depths'};
td(21,:)={'CoordScalar','int16',70,'','Scalar to be applied to all coordinates'};
td(22,:)={'SrcX','int32',72,'CoordScalar','Source coordinate X'};
td(23,:)={'SrcY','int32',76,'CoordScalar','Source coordinate Y'};
td(24,:)={'GroupX','int32',80,'CoordScalar','Group coordinate X'};
td(25,:)={'GroupY','int32',84,'CoordScalar','Group coordinate Y'};
td(26,:)={'CoordUnits','int16',88,'','Coordinate units'};
td(27,:)={'WeatherVel','int16',90,'','Weathering velocity'};
td(28,:)={'SubWeatherVel','int16',92,'','Subweathering velocity'};
if segyrev==0
    td(29,:)={'SrcUphole','int16',94,'','Uphole time at source in milliseconds'};
    td(30,:)={'GroupUphole','int16',96,'','Uphole time at group in milliseconds'};
    td(31,:)={'SrcStatic','int16',98,'','Source static correction in milliseconds'};
    td(32,:)={'GroupStatic','int16',100,'','Group static correction in milliseconds'};
    td(33,:)={'TotalStatic','int16',102,'','Total static applied in milliseconds'};
    td(34,:)={'LagTimeA','int16',104,'','Lag time A'};
    td(35,:)={'LagTimeB','int16',106,'','Lag Time B'};
    td(36,:)={'RecTimeDelay','int16',108,'','Delay recording time'};
    td(37,:)={'MuteStart','int16',110,'','Mute start time in milliseconds'};
    td(38,:)={'MuteEnd','int16',112,'','Mute end time in milliseconds'};
else
    td(29,:)={'SrcUphole','int16',94,'TimeScalar','Uphole time at source in milliseconds'};
    td(30,:)={'GroupUphole','int16',96,'TimeScalar','Uphole time at group in milliseconds'};
    td(31,:)={'SrcStatic','int16',98,'TimeScalar','Source static correction in milliseconds'};
    td(32,:)={'GroupStatic','int16',100,'TimeScalar','Group static correction in milliseconds'};
    td(33,:)={'TotalStatic','int16',102,'TimeScalar','Total static applied in milliseconds'};
    td(34,:)={'LagTimeA','int16',104,'TimeScalar','Lag time A'};
    td(35,:)={'LagTimeB','int16',106,'TimeScalar','Lag Time B'};
    td(36,:)={'RecTimeDelay','int16',108,'TimeScalar','Delay recording time'};
    td(37,:)={'MuteStart','int16',110,'TimeScalar','Mute start time in milliseconds'};
    td(38,:)={'MuteEnd','int16',112,'TimeScalar','Mute end time in milliseconds'};
end
td(39,:)={'SampThisTrc','uint16',114,'','*Number of samples in this trace'};
td(40,:)={'SampRateThisTrc','uint16',116,'','*Sample interval in microseconds for this trace'};
td(41,:)={'GainType','int16',118,'','Gain type of field instruments'};
td(42,:)={'InstGainConst','int16',120,'','Instrument gain constant (dB)'};
td(43,:)={'InstEarlyGain','int16',122,'','Instrument early or initial gain (dB)'};
td(44,:)={'Correlated','int16',124,'','Correlated'};
td(45,:)={'SweepStartFreq','int16',126,'','Sweep frequency at start (Hz)'};
td(46,:)={'SweepEndFreq','int16',128,'','Sweep frequency at end (Hz)'};
td(47,:)={'SweepLength','int16',130,'','Sweep length in milliseconds'};
td(48,:)={'SweepType','int16',132,'','Sweep type'};
td(49,:)={'SweepStartTaper','int16',134,'','Sweep trace taper length at start in milliseconds'};
td(50,:)={'SweepEndTaper','int16',136,'','Sweep trace taper length at end in milliseconds'};
td(51,:)={'SweepTaperType','int16',138,'','Taper type'};
td(52,:)={'AliasFiltFreq','int16',140,'','Alias filter frequency (Hz)'};
td(53,:)={'AliasFileSlope','int16',142,'','Alias filter slope (dB/octave)'};
td(54,:)={'NotchFiltFreq','int16',144,'','Notch filter frequency (Hz)'};
td(55,:)={'NotchFiltSlope','int16',146,'','Notch filter slope (dB/octave)'};
td(56,:)={'LowCutFreq','int16',148,'','Low-cut frequency (Hz)'};
td(57,:)={'HighCutFreq','int16',150,'','High-cut frequency (Hz)'};
td(58,:)={'LowCutSlope','int16',152,'','Low-cut slope (dB/octave)'};
td(59,:)={'HighCutSlope','int16',154,'','High-cut slope (dB/octave)'};
td(60,:)={'Year','int16',156,'','Year data recorded'};
td(61,:)={'Day','int16',158,'','Day of year'};
td(62,:)={'Hour','int16',160,'','Hour of day'};
td(63,:)={'Minute','int16',162,'','Minute of hour'};
td(64,:)={'Second','int16',164,'','Second of minute'};
td(65,:)={'TimeBasisCode','int16',166,'','Time basis code'};
td(66,:)={'TrcWeightFactor','int16',168,'','Trace weighting factor'};
td(67,:)={'GroupNumRollSw1','int16',170,'','Geophone group number of roll switch position one'};
td(68,:)={'GroupNumTrc1','int16',172,'','Geophone group number of trace number one within original field record'};
td(69,:)={'GroupNumTrcN','int16',174,'','Geophone group number of last trace within original field record'};
td(70,:)={'GapSize','int16',176,'','Gap size'};
td(71,:)={'Overtravel','int16',178,'','Over travel associated with taper at beginning or end of line'};

if segyrev == 0
    td = [td; mkunassigned(1,15,180)];
elseif segyrev == 1 || segyrev == 2    
    td(72,:)={'CdpX','int32',180,'CoordScalar','X coordinate of ensemble (CDP) position of this trace (scalar in Trace Header bytes 71-72 applies).'};
    td(73,:)={'CdpY','int32',184,'CoordScalar','Y coordinate of ensemble (CDP) position of this trace (scalar in bytes Trace Header 71-72 applies).'};
    td(74,:)={'InlineNum','int32',188,'','For 3-D poststack data, this field should be used for the in-line number.'};
    td(75,:)={'XlineNum','int32',192,'','For 3-D poststack data, this field should be used for the cross-line number.'};
    td(76,:)={'ShotPointNum','int32',196,'ShotPointScalar','Shotpoint number — This is probably only applicable to 2-D poststack data.'};
    td(77,:)={'ShotPointScalar','int16',200,'','Scalar to be applied to the shotpoint number in Trace Header bytes 197-200 to give the real value.'};
    td(78,:)={'TrcValMeasUnit','int16',202,'','Trace value measurement unit'};
    td(79,:)={'TransConstMant','int32',204,'','Transduction Constant mantissa'};
    td(80,:)={'TransConstExp','int16',208,'','Transduction Constant exponent'};
    td(81,:)={'TransUnits','int16',210,'','Transduction Units'};
    td(82,:)={'DeviceID','int16',212,'','Device/Trace Identifier'};
    td(83,:)={'TimeScalar','int16',214,'','Scalar to be applied to times'};
    td(84,:)={'SrcTypeOrient','int16',216,'','Source Type/Orientation'};
    td(85,:)={'VerticalSrcDir','int16',218,'','Vertical source energy direction with respect to the source orientation'};
    td(86,:)={'XlineSrcDir','int16',220,'','Cross-line source energy direction with respect to the source orientation'};
    td(87,:)={'InlineSrcDir','int16',222,'','In-line source energy direction with respect to the source orientation'};    
    td(88,:)={'SrcMeasureMant','int32',224,'','Source Measurement mantissa'};
    td(89,:)={'SrcMeasureExp','int16',228,'','Source Measurement exponent'};    
    td(90,:)={'SrcMeasureUnit','int16',230,'','Source Measurement Unit'};
    if segyrev == 1    
        td = [td; mkunassigned(1,2,232)];
    elseif segyrev == 2
        td(91,:)={'TrcHdrName','uint64',232,'','Trace Header Name'};
        if numexthdr
            td(92,:)={'ExtTrcNumLine','uint64',240,'','Extended Trace sequence number within line'};
            td(93,:)={'ExtTrcNumReel','uint64',248,'','Extended Trace sequence number within SEGY file'};
            td(94,:)={'ExtFieldRecNum','int64',256,'','Extended Original field record number'};
            td(95,:)={'ExtEnsembleNum','int64',264,'','Extended Ensemble number'};
            td(96,:)={'ExtGroupElev','ieee64',272,'','Extended Receiver group elevation'};
            td(97,:)={'ExtGroupDepth','ieee64',280,'','Extended Receiver group depth below surface'};            
            td(98,:)={'ExtSrcElev','ieee64',288,'','Extended Surface elevation at source'};
            td(99,:)={'ExtSrcDepth','ieee64',296,'','Extended Source depth below surface'};
            td(100,:)={'ExtGroupDatumElev','ieee64',304,'','Extended Datum elevation at receiver group'};
            td(101,:)={'ExtSrcDatumElev','ieee64',312,'','Extended Datum elevation at source'};
            td(102,:)={'ExtSrcWaterColHeight','ieee64',320,'','Extended Water column height at source'};
            td(103,:)={'ExtGroupWaterColHeight','ieee64',328,'','Extended Water column height at group'};
            td(104,:)={'ExtSrcX','ieee64',336,'','Extended Source coordinate X'};
            td(105,:)={'ExtSrcY','ieee64',344,'','Extended Source coordinate Y'};
            td(106,:)={'ExtGroupX','ieee64',352,'','Extended Group coordinate X'};
            td(107,:)={'ExtGroupY','ieee64',360,'','Extended Group coordinate Y'};
            td(108,:)={'ExtSrcRecOffset','ieee64',368,'','Extended Distance from center of the source point to the center of the receiver group'};
            td(109,:)={'ExtSampThisTrc','uint32',376,'','Extended Number of samples in this trace'};
            td(110,:)={'NanoSeconds','int32',380,'','Nanoseconds to add to second of minute'};
            td(111,:)={'ExtSampRateThisTrc','ieee64',384,'','Extended Sample interval in microseconds for this trace'};
            td(112,:)={'CableNum','int32',392,'','Cable number or Recording device/sensor ID number'};
            td(113,:)={'LastTrcFlag','uint16',396,'','Last trace flag'};
            td(114,:)={'ExtCdpX','ieee64',398,'','X coordinate of ensemble (CDP) position of this trace'};
            td(115,:)={'ExtCdpY','ieee64',406,'','Y coordinate of ensemble (CDP) position of this trace'};
            td(116,:)={'Unassigned01','uint16',414,'','Unassigned 1'};
            td = [td; mkunassigned(2,15,416)];
            td(131,:)={'ExtTrcHdrName','uint64',472,'','Extended Trace Header Name'};
        end
    end  
else
    error('@Trace/newDefinition: Unknown SEG-Y revision number');
end

end %end function

function ca = mkunassigned(widx1,widx2,sbyte)
    nrow = widx2-widx1+1;
    ca = cell(nrow,5);
    cidx=1;
    for ii=widx1:widx2
        ca(cidx,:)={sprintf('Unassigned%02d',ii),'uint32',sbyte,'',...
            sprintf('Unassigned %d',ii)};
        cidx=cidx+1;
        sbyte=sbyte+4;
    end
end