%% function to process interpolated data, instrument type by instrument type
% 'ProcessedData' contains interpolated observations grouped by instrument
% 'ProfilesBSNE' contains information about BSNE flux profiles
% 'InstrumentMetadata' contains metadata about instrument time intervals
% Output 'ProcessedData' is built on 'InterpolatedData',
% but it adds/removes fields to account for interpolation
% Also, 'ProcessedData' contains added column 'TotalFlux,' which is comes
% from interpolation of Wenglor data
% Dependencies: NONE
% Used by: Processing_Master
 
function [ProcessedData] = ProcessData(InterpolatedData, WeightBSNE, InstrumentMetadata, GrainSize_BSNE, GrainSize_Surface)

%% PARAMETERS
%set time increment for TotalFlux calculations
TotalFlux_dt = duration(0,0,0.04); %set to be 0.04 seconds (25 Hz), which was used for all data analysis
InstrumentHeightCorrection_m = 0; %amount by which to adjust instrument heights to account for scour beneath Wenglors

%% GO THROUGH INSTRUMENT TYPES, DETERMINE ABSOLUTE HEIGHTS BASED ON DISTANCE SENSOR READINGS
InstrumentTypes = fieldnames(InterpolatedData); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type 
    Instruments = fieldnames(InterpolatedData.(InstrumentType)); %get list of instruments
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        Instrument = Instruments{j} %get current instrument
        N_Intervals = length(InterpolatedData.(InstrumentType).(Instrument)); %how many time intervals for this instrument
        RefHeight_List = zeros(N_Intervals,1); %create list of reference heights (m)
        RefHeightErr_List = zeros(N_Intervals,1); %create list of reference heights errors (m)
        
        for k = 1:N_Intervals
            StartTime = InterpolatedData.(InstrumentType).(Instrument)(k).StartTime; %get start time
            EndTime = InterpolatedData.(InstrumentType).(Instrument)(k).EndTime; %get end time
            StartHeight_m = InterpolatedData.(InstrumentType).(Instrument)(k).StartHeight_m; %get start height
            EndHeight_m = InterpolatedData.(InstrumentType).(Instrument)(k).EndHeight_m; %get end height
            HeightRef = InterpolatedData.(InstrumentType).(Instrument)(k).HeightRef; %get reference variable for height
            HeightErr_m = InterpolatedData.(InstrumentType).(Instrument)(k).HeightErr_m; %get height error for instrument
            
            %compute reference height
            %if it is 0, RefHeight_m and RefHeightErr_m are 0
            if strcmp(HeightRef,'0')
                RefHeight_m = 0;
                RefHeightErr_m = 0;
            
            %otherwise extract data from distance sensor specified by this string
            else
               
                %extract distance sensor heights (mm) for this interval
                z_interval_mm = ExtractVariableTimeInterval(InterpolatedData.Distance.(HeightRef),StartTime,EndTime,'z','int','int');
                
                %if values are contained in this interval, compute mean of them 
                if ~isempty(z_interval_mm)
                    RefHeight_m = mean(z_interval_mm)/1000+InstrumentHeightCorrection_m; %compute mean, convert to meters.  Add in height correction value
                    RefHeightErr_m = std(z_interval_mm)/1000; %error is standard deviation, convert to meters
                    
                %if no values contained in interval, use values from previous interval
                else
                    RefHeight_m = RefHeight_List(k-1);
                    RefHeightErr_m = RefHeightErr_List(k-1);
                end
            end

            %assign values to lists
            RefHeight_List(k) = RefHeight_m;
            RefHeightErr_List(k) = RefHeightErr_m;
            
            %compute instrument height
            z_Instrument = RefHeight_m+mean([StartHeight_m; EndHeight_m]);
            
            %compute error as combination of measurement error, reference
            %height error, and difference of start / end heights
            z_Err = sqrt(HeightErr_m^2 + RefHeightErr_m^2 + (StartHeight_m-EndHeight_m).^2);
            
            %assign height to structured array
            InterpolatedData.(InstrumentType).(Instrument)(k).InstrumentHeight.z = z_Instrument;
            InterpolatedData.(InstrumentType).(Instrument)(k).InstrumentHeight.err = z_Err;
            InterpolatedData.(InstrumentType).(Instrument)(k).InstrumentHeight.Units = 'm';
        end
        
    %reorder fields to move height column to beginning of list
    N_fields = length(fieldnames(InterpolatedData.(InstrumentType).(Instrument)));
    InterpolatedData.(InstrumentType).(Instrument) = orderfields(InterpolatedData.(InstrumentType).(Instrument), [N_fields,1:(N_fields-1)]);        
    
    end
end


%% EXTRACT BSNE FLUX PROFILES BASED ON BSNE WEIGHTS
ProfilesBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE,GrainSize_Surface); %process BSNE data to get profile values


%% GO THROUGH BSNE PROFILES TO CALIBRATE WENGLOR COUNTS TO HEIGHT-SPECIFIC FLUXES
%initialize ProcessedWenglors from raw (interpolated) Wenglor data
ProcessedWenglors = InterpolatedData.Wenglor;

%get names of Wenglors
WenglorNames = fieldnames(ProcessedWenglors);
N_Wenglors = length(WenglorNames);

%inialize flux and calibration entries for Wenglor structured array
for i=1:N_Wenglors
    
    %compute number of intervals for Wenglor
    N_Intervals = length(ProcessedWenglors.(WenglorNames{i}));
    
    %go through each interval, initialize flux column
    for j=1:N_Intervals
        
        %get number of time steps for flux calc
        N_t = length(ProcessedWenglors.(WenglorNames{i})(j).t.int);
        
        %assign 'flux' structured array to new column of master stuctured array
        ProcessedWenglors.(WenglorNames{i})(j).flux = struct('qz',zeros(N_t,1),'units','g/m^2/s','qzPerCount',zeros(N_t,1));
    end
    
    %rearrange fields
    N_fields = length(fieldnames(ProcessedWenglors.(WenglorNames{i})));
    ProcessedWenglors.(WenglorNames{i}) = orderfields(ProcessedWenglors.(WenglorNames{i}), [1:3,N_fields,4:(N_fields-1)]);        
end

% % ** Generate matrix of Wenglor calibration values ** %
% WenglorCountsMatrix = zeros(45,9)*NaN;
% WenglorCalibrationMatrix = zeros(45,9)*NaN;
% FluxPredictionMatrix = zeros(45,9)*NaN;
% % ** *
% ** NEED TO DETERMINE WHY CALIBRATION VALUE TENDS TO INCREASE TOWARD TOP
% OF PROFILE (SOMEHOW TOP WENGLORS ARE UNDERCOUNTING RELATIVE TO
% PREDICTIONS -- IS IT A GRAIN SIZE EFFECT?)

%go through BSNE profiles for calibration
for i=1:length(ProfilesBSNE);

    %get start and end time of BSNE profile
    StartTimeBSNE_Primary = ProfilesBSNE(i).StartTime %print on command window to track
    EndTimeBSNE_Primary = ProfilesBSNE(i).EndTime;

    %get BSNE flux information for time interval
    q0_BSNE = ProfilesBSNE(i).q0; %qz at 0, (g/m^2/s)
    zbar_BSNE = ProfilesBSNE(i).zbar.mid; %get mean BSNE height, (m)
    
    %create enlarged time intervals for calibrating Wenglor times outside of BSNE profile times
    if i == 1 %earliest start time is beginning of first BSNE day
        StartTimeBSNE_Enlarged = ProfilesBSNE(i).Date;
        EndTimeBSNE_Enlarged = mean([ProfilesBSNE(i).EndTime, ProfilesBSNE(i+1).StartTime]);
    elseif i == length(ProfilesBSNE) %latest start time is end of last BSNE day
        StartTimeBSNE_Enlarged = mean([ProfilesBSNE(i-1).EndTime, ProfilesBSNE(i).StartTime]);
        EndTimeBSNE_Enlarged = ProfilesBSNE(i).Date + 1;
    else %middle start and end times are spaced halfway between BSNE end time before and start time after
        StartTimeBSNE_Enlarged = mean([ProfilesBSNE(i-1).EndTime, ProfilesBSNE(i).StartTime]);
        EndTimeBSNE_Enlarged = mean([ProfilesBSNE(i).EndTime, ProfilesBSNE(i+1).StartTime]);
    end

    %go through each Wenglor for calibration
    for j = 1:N_Wenglors
        
        %get Wenglor counts and interval number (there should be only one) for primary interval
        [PrimaryInterval_n, PrimaryInterval_t, PrimaryInterval_IntervalNumber, ~] = ...
            ExtractVariableTimeInterval(ProcessedWenglors.(WenglorNames{j}),StartTimeBSNE_Primary,EndTimeBSNE_Primary,'n','int','int');
        
        %perform calculations assuming single set of Wenglor counts within interval
        if length(PrimaryInterval_IntervalNumber)==1

            %extract information about Wenglor for primary interval
            z_Wenglor = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).InstrumentHeight.z; %get Wenglor height for primary interval
            ErrorCode_Wenglor = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).ErrorCode; %get ErrorCode for Wenglor

            %compute expected flux at Wenglor height
            qz_pred = q0_BSNE.*exp(-z_Wenglor/zbar_BSNE); %(g/m^2/s)

            %compute counts per second for Wenglor during BSNE time interval
            W_CountsPerSecond = sum(PrimaryInterval_n)/seconds(EndTimeBSNE_Primary - StartTimeBSNE_Primary);

            %get time increment for Wenglor observations
            W_dt = seconds(mode(diff(PrimaryInterval_t)));

            %compute conversion factor from Wenglor counts to flux, "qzPerCount"
            qzPerCount = qz_pred/(W_CountsPerSecond*W_dt);

%             % ** Generate matrix of Wenglor calibration values ** %
%             WenglorCountsMatrix(i,j) = W_CountsPerSecond;
%             WenglorCalibrationMatrix(i,j) = qzPerCount;
%             FluxPredictionMatrix(i,j) = qz_pred;
%             % ** *
            
            %get interval numbers for enlarged interval subject to this calibration (there may be more than one)
            [~, ~, EnlargedInterval_IntervalNumbers, EnlargedInterval_IntervalIndices] = ...
                ExtractVariableTimeInterval(ProcessedWenglors.(WenglorNames{j}),StartTimeBSNE_Enlarged,EndTimeBSNE_Enlarged,'n','int','int');

            %go through each of these intervals and perform calibration
            N_IntervalNumbers = length(EnlargedInterval_IntervalNumbers);

            for k = 1:N_IntervalNumbers

                %get list of counts for calibration
                n_list = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).n.int(EnlargedInterval_IntervalIndices{k});

                %convert to flux values
                qz_list = qzPerCount*n_list;

                %assign fluxes and calibration factor to structured array
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).flux.qz(EnlargedInterval_IntervalIndices{k}) = qz_list;
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).flux.qzPerCount(EnlargedInterval_IntervalIndices{k}) = qzPerCount;
            end
        end
    end
end

%% COMPUTE TOTAL FLUXES FROM CALIBRATED WENGLORS

%get indices of valid Wenglor times in metadata table
ind_MetadataWenglors = find(strcmp(InstrumentMetadata.InstrumentType,'Wenglor')&(InstrumentMetadata.ErrorCode==0));

%get raw start and end times for total flux calcs (i.e., start and end times of all possible intervals without any sorting)
TotalFlux_RawStartTimes = InstrumentMetadata.StartTime(ind_MetadataWenglors);
TotalFlux_RawEndTimes = InstrumentMetadata.EndTime(ind_MetadataWenglors);
TotalFlux_RawTimes = unique(union(TotalFlux_RawStartTimes,TotalFlux_RawEndTimes));

%compute start times of overlapping intervals
[TotalFlux_OverlappingStartTimes, TotalFlux_OverlappingEndTimes] = CombineIntervals(TotalFlux_RawStartTimes, TotalFlux_RawEndTimes);
N_OverlappingIntervals = length(TotalFlux_OverlappingStartTimes);

%now, subdivide these into all distinctive intervals for computing total flux
TotalFlux_StartTimes = []; %initialize list of all distinctive start times
TotalFlux_EndTimes = []; %initialize list of all distinctive end times
for i = 1:N_OverlappingIntervals
    OverlappingInterval_AllStartTimes = unique(TotalFlux_RawTimes(...
        TotalFlux_RawTimes>=TotalFlux_OverlappingStartTimes(i)&...
        TotalFlux_RawTimes<TotalFlux_OverlappingEndTimes(i)));
    TotalFlux_StartTimes = [TotalFlux_StartTimes; OverlappingInterval_AllStartTimes];
    OverlappingInterval_AllEndTimes = unique(TotalFlux_RawTimes(...
        TotalFlux_RawTimes>TotalFlux_OverlappingStartTimes(i)&...
        TotalFlux_RawTimes<=TotalFlux_OverlappingEndTimes(i)));
    TotalFlux_EndTimes = [TotalFlux_EndTimes; OverlappingInterval_AllEndTimes];
end

%extract associated dates for entries
[y,m,d] = ymd(TotalFlux_StartTimes);
TotalFlux_Dates = datetime(y,m,d);

%extract associated sites for entries
N_TotalFluxIntervals = length(TotalFlux_StartTimes);
TotalFlux_Sites = cell(N_TotalFluxIntervals,1);
for i=1:N_TotalFluxIntervals
    TotalFlux_Sites{i} = char(unique(InstrumentMetadata.Site(InstrumentMetadata.Date==TotalFlux_Dates(i))));
end

%initialize structured array for total flux from profile
TotalFlux = struct(...
    't',struct('calc',[]),...
    'Q',struct('calc',[],'Units','g/m/s'),...
    'qz',struct('calc',[],'qzPerCount',[],'WenglorNames',[],'Units','g/m^2/s'),...
    'z',struct('calc',[],'Units','m'),...
    'Site', cellstr(TotalFlux_Sites),...
    'Date', num2cell(TotalFlux_Dates),...
    'StartTime', num2cell(TotalFlux_StartTimes),...
    'EndTime', num2cell(TotalFlux_EndTimes));

%figure out which Wenglors are contained in each interval
for i = 1:N_Wenglors
    %***NEED TO FIX THIS -- EXTRA WENGLORS ARE SLIPPING THROUGH, IT SEEMS
    ind_MetadataWenglor = find(strcmp(InstrumentMetadata.Instrument,WenglorNames{i})&(InstrumentMetadata.ErrorCode==0));
    Wenglor_StartTimes = InstrumentMetadata.StartTime(ind_MetadataWenglor);
    Wenglor_EndTimes = InstrumentMetadata.EndTime(ind_MetadataWenglor);
    for j = 1:N_TotalFluxIntervals
        if ~isempty(find((Wenglor_StartTimes<=TotalFlux(j).StartTime)&(Wenglor_EndTimes>=TotalFlux(j).EndTime)))
            TotalFlux(j).qz.WenglorNames = [TotalFlux(j).qz.WenglorNames, {WenglorNames{i}}];
        end
    end
end

%fill in times and initialize components of structured array
for i = 1:N_TotalFluxIntervals
    %print data for script tracking purposes
    TotalFlux(i).Date
    
    %calculate times
    t_calc = (TotalFlux(i).StartTime:TotalFlux_dt:TotalFlux(i).EndTime)';
    N_t = length(t_calc);
    TotalFlux(i).t.calc = t_calc;
    
    %initialize total fluxes
    TotalFlux(i).Q.calc = zeros(N_t,1);
    
    %get information about Wenglors in profile
    N_qz = length(TotalFlux(i).qz.WenglorNames);
    TotalFlux(i).z.calc = zeros(1,N_qz);

    %initialize partial fluxes
    TotalFlux(i).qz.calc = zeros(N_t,N_qz)*NaN;
    TotalFlux(i).qz.qzPerCount = zeros(N_t,N_qz)*NaN;
    
    %fill in values
    for j = 1:N_qz
        WenglorName = TotalFlux(i).qz.WenglorNames{j}; %get specific Wenglor name
        
        %extract flux information
        [qz, qz_t, qz_IntervalNumber, qz_ind] = ...
            ExtractVariableTimeInterval(ProcessedWenglors.(WenglorName),TotalFlux(i).StartTime,TotalFlux(i).EndTime,'flux','qz','int');
        
        %line things up with times in flux table
        qz_StartInd(j) = find(t_calc==qz_t(1)); %index in table of first qz from Wenglor, keep track of it
        qz_EndInd(j) = find(t_calc==qz_t(end)); %index in table of last qz from Wenglor, keep track of it
        
        %submit qz's into flux table
        TotalFlux(i).qz.calc(qz_StartInd(j):qz_EndInd(j),j) = qz;
        
        %submit calibration values into flux table
        TotalFlux(i).qz.qzPerCount(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).flux.qzPerCount(cell2mat(qz_ind));
        
        %get information about Wenglor height
        TotalFlux(i).z.calc(j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).InstrumentHeight.z;
    end
    
    %sort values by z, compute order of data columns from this
    [z_sorted, z_sort_ind] = sort(TotalFlux(i).z.calc);

    %get indices of rows that are contained in data to keep
    row_keep_ind = max(qz_StartInd):min(qz_EndInd);
    
    %modify matrices based on sorting and determination of rows to keep
    TotalFlux(i).z.calc = z_sorted;
    TotalFlux(i).qz.WenglorNames = TotalFlux(i).qz.WenglorNames(z_sort_ind);
    TotalFlux(i).t.calc = TotalFlux(i).t.calc(row_keep_ind);
    TotalFlux(i).Q.calc = TotalFlux(i).Q.calc(row_keep_ind);
    TotalFlux(i).qz.calc = TotalFlux(i).qz.calc(row_keep_ind, z_sort_ind);
    TotalFlux(i).qz.qzPerCount = TotalFlux(i).qz.qzPerCount(row_keep_ind, z_sort_ind);
    
    %For Q integration, compute dz values based on unique heights
    dz_profile = zeros(N_qz,1); %dz for qz used in integration for Q
    z_profile = TotalFlux(i).z.calc; %get Wenlgor heights for profile
    z_unique = unique(TotalFlux(i).z.calc); %get unique heights (in case of repeats)
    
    %get bounding bottom and top z's for Q integration dz values
    for j = 1:length(z_unique)
        z_ind = find(z_profile==z_unique(j)); %get indices of heights
        
        %if lowest, use 0 as bottom height
        if j == 1;
            zbot = 0;
            ztop = sqrt(z_unique(j)*z_unique(j+1));
        
        %if highest, set top height so that Wenglor height is at geometric mean of its box
        elseif j == length(z_unique);
            zbot = sqrt(z_unique(j-1)*z_unique(j)); 
            ztop = z_unique(j).^2/zbot;
        
        %otherwise, dz range is set at geometric mean of heights above and below 
        else
            zbot = sqrt(z_unique(j-1)*z_unique(j));
            ztop = sqrt(z_unique(j)*z_unique(j+1));
        end
        
        %dz is the range of z values, divided by the number of Wenglors at this height
        dz_profile(z_ind) = (ztop-zbot)/length(z_ind);
    end
    
    %Perform integration to get Q
    TotalFlux(i).Q.calc = TotalFlux(i).qz.calc*dz_profile;
    %** ISSUE: THIS INTEGRATION IS GIVING TOO-LOW VALUES FOR Q, BECAUSE IT
    %MISSES FLUXES NEAR BOTTOM OF PROFILE
end


%% RENAME "INTERPOLATEDDATA" TO "PROCESSEDDATA", ADD TOTAL FLUX AND REVISED WENGLOR INFORMATION TO STRUCTURED ARRAY

%now that everything is done, rename "InterpolatedData" as "ProcessedData"
ProcessedData = InterpolatedData;

%add BSNE data to processedData structured array
ProcessedData.BSNE = ProfilesBSNE;

%add total flux to processed data structured array
ProcessedData.TotalFlux = TotalFlux;

%add revised Wenglor structured array to processed data structured array
ProcessedData.Wenglor = ProcessedWenglors;

%add grain-size data to processed data structured array
ProcessedData.GrainSize = struct('Surface',GrainSize_Surface,'BSNE',GrainSize_BSNE);