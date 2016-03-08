function [FIsec, FIyyyy, EQsec, Omarker] = getFileAndEQseconds(F,eqin,offset)
%calculate start times of the files in seconds after midnight, January 1st
%this works for SAC files created with rdseed4.5.1
%eg: F = '1993.159.23.15.09.7760.IU.KEV..BHN.D.SAC'
% if your filnames contains no julian day, please use command
% dayofyear (in Splitlab/Tools)
% 


% Windows user can try a renamer , for example 1-4aren (one-for all renamer)
% http://www.1-4a.com/rename/ perhaps this adress is still valid

global config


if config.UseHeaderTimes | strcmp(config.FileNameConvention, '*.e; *.n; *.z')
    for k=1:size(F,1)
        workbar(k/size(F,1),'Reading header')
        try
            sac = rsac([config.datadir filesep F(k,:)]);
        catch
            sac = rsacsun([config.datadir filesep F(k,:)]);
        end
        [FIyyyy(k), FIddd(k), FIHH(k), FIMM(k), FISS(k)] =...
            lh(sac, 'NZYEAR','NZJDAY','NZHOUR','NZMIN', 'NZSEC'); 
        Omarker(k) = lh(sac, 'O');
    end
    
     Omarker(Omarker == -12345) = 0;  %verify, if O-marker is set   
     FIsec  =  FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400 + Omarker;

    fclose all;



else % USE FILENAME
    switch config.FileNameConvention
        case 'RDSEED'
            % RDSEED format '1993.159.23.15.09.7760.IU.KEV..BHN.D.SAC' 
            FIyyyy = str2num(F(:,1:4));
            FIddd  = str2num(F(:,6:8));
            FIHH   = str2num(F(:,10:11));
            FIMM   = str2num(F(:,13:14));
            FISS   = str2num(F(:,16:17));
            FIMMMM = str2num(F(:,18:22));
            FIsec  = FIMMMM + FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400;



        case 'SEISAN'
            % SEISAN format '2003-05-26-0947-20S.HOR___003_HORN__BHZ__SAC'
            FIyyyy = str2num(F(:,1:4));
            FImonth= str2num(F(:,6:7));
            FIdd   = str2num(F(:,9:10));
            FIHH   = str2num(F(:,12:13));
            FIMM   = str2num(F(:,14:15));
            FISS   = str2num(F(:,17:18));

            FIddd = dayofyear(FIyyyy',FImonth',FIdd')';%julian Day
            FIsec =  FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400;


        case 'YYYY.JJJ.hh.mm.ss.stn.sac.e'
            %  Format: 1999.136.15.25.00.ATD.sac.z
            FIyyyy = str2num(F(:,1:4));
            FIddd  = str2num(F(:,6:8));
            FIHH   = str2num(F(:,10:11));
            FIMM   = str2num(F(:,13:14));
            FISS   = str2num(F(:,16:17));
            FIsec  = FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400;

        case 'YYYY.MM.DD-hh.mm.ss.stn.sac.e';
            % Format: 2003.10.07-05.07.15.DALA.sac.z
            FIyyyy = str2num(F(:,1:4));
            FImonth= str2num(F(:,6:7));
            FIdd   = str2num(F(:,9:10));
            FIHH   = str2num(F(:,12:13));
            FIMM   = str2num(F(:,15:16));
            FISS   = str2num(F(:,18:19));

            FIddd = dayofyear(FIyyyy',FImonth',FIdd')';%julian Day
            FIsec  =  FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400;
            
        case 'YYYY_MM_DD_hhmm_stnn.sac.e';
            % Format: 2005_03_02_1155_pptl.sac (LDG/CEA data)
            FIyyyy = str2num(F(:,1:4));
            FImonth= str2num(F(:,6:7));
            FIdd   = str2num(F(:,9:10));
            FIHH   = str2num(F(:,12:13));
            FIMM   = str2num(F(:,14:15));
            
            FIddd = dayofyear(FIyyyy',FImonth',FIdd')';%julian Day
            FIsec = FIMM*60 + FIHH*3600 + (FIddd)*86400;

        case 'stn.YYMMDD.hhmmss.e'
            % Format: fp2.030723.213056.X (BroadBand OBS data)
            FIyyyy = 2000 + str2num(F(:,5:6));%only two-digit year identifier => add 2000, assuming no OBS data before 2000
            FImonth= str2num(F(:,7:8));
            FIdd   = str2num(F(:,9:10));
            FIHH   = str2num(F(:,12:13));
            FIMM   = str2num(F(:,14:15));
            FISS   = str2num(F(:,16:17));

            FIddd = dayofyear(FIyyyy',FImonth',FIdd')';%julian Day
            FIsec = FISS + FIMM*60 + FIHH*3600 + (FIddd)*86400;
    end
    
    Omarker = zeros(size(FIsec));
end

%% get earthquake origin times
for a=1:length(eqin);
    EQsec(a) = eqin(a).date(6) + eqin(a).date(5)*60 + eqin(a).date(4)*3600 + eqin(a).date(7)*86400;
end

EQsec = EQsec + offset;

