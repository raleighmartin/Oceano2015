%% function to read in raw data from logger files, assign metadata, and aggregate by instrument
% Dependencies: GetDataTypes, GetTimes, ReadTable
% Used by: Processing_Master

function [RawData] = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData)

%% Initialize structured array
RawData = struct; %initialize array
InstrumentTypes = unique(InstrumentMetadata.InstrumentType); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type
    RawData.(InstrumentType) = struct; %initialize nested structured array
    Instruments = unique(InstrumentMetadata.Instrument(strcmp(InstrumentMetadata.InstrumentType,InstrumentType))); %get list of instruments for this type
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        %get current instrument
        Instrument = Instruments{j};
        
        %find indices of InstrumentMetadata associated with this instrument
        InstrumentMetadata_ind = find(strcmp(InstrumentMetadata.Instrument,Instrument)); 
        N_Intervals = length(InstrumentMetadata_ind); %get number of intervals
        
        %populate metadata values for instrument into structured array
        RawData.(InstrumentType).(Instrument) = struct(...
            'Site', cellstr(InstrumentMetadata.Site(InstrumentMetadata_ind)),...
            'Date', num2cell(InstrumentMetadata.Date(InstrumentMetadata_ind)),...
            'StartTime', num2cell(InstrumentMetadata.StartTime(InstrumentMetadata_ind)),...
            'EndTime', num2cell(InstrumentMetadata.EndTime(InstrumentMetadata_ind)),...
            'InstrumentType', cellstr(InstrumentMetadata.InstrumentType(InstrumentMetadata_ind)),...    
            'Instrument', cellstr(InstrumentMetadata.Instrument(InstrumentMetadata_ind)),...    
            'StartHeight_m', num2cell(InstrumentMetadata.StartHeight_m(InstrumentMetadata_ind)),...
            'EndHeight_m', num2cell(InstrumentMetadata.EndHeight_m(InstrumentMetadata_ind)),...
            'HeightErr_m', num2cell(InstrumentMetadata.HeightErr_m(InstrumentMetadata_ind)),...
            'HeightRef', cellstr(InstrumentMetadata.HeightRef(InstrumentMetadata_ind)),...
            'Longitudinal_m', num2cell(InstrumentMetadata.Longitudinal_m(InstrumentMetadata_ind)),...
            'Spanwise_m', num2cell(InstrumentMetadata.Spanwise_m(InstrumentMetadata_ind)),...
            'AngleErr_deg', num2cell(InstrumentMetadata.AngleErr_deg(InstrumentMetadata_ind)),...
            'ErrorCode', num2cell(InstrumentMetadata.ErrorCode(InstrumentMetadata_ind)));
    end
end

%% Go through all data tables
N_tables = length(LoggerTables.Table); % How many data tables to import?

for i=1:N_tables
    %get info about filepath for table, display to show progress of processing
    tableFilepath = strcat(...
        folder_LoggerRawData,...
        LoggerTables.Site{i},'_',...
        datestr(LoggerTables.Date(i),'yyyy-mm-dd'),'_',...
        LoggerTables.Logger{i},'_',...
        LoggerTables.Table{i},'.DAT')
        
    %get information about variables for table
    tableVariableList = strsplit(LoggerTables.Variables{i},', '); %determine variables in logger file
    
    %get additional information about table
    tableTypeList = GetDataTypes(tableVariableList,InstrumentVariables); %get type for each variable
    tableTimes = GetTimes(LoggerTables.Site(i),LoggerTables.Date(i),LoggerTables.Logger(i),LoggerTimes);
  
    %run script to read data from logger file, output as element in cell array of logger data
    [tableTimestamp, tableRecords, tableVariables] = ReadTable(tableFilepath,tableVariableList,tableTypeList,tableTimes);
    
    %find instrument metadata intervals associated with this logger table
    Variables = fieldnames(tableVariables); %get list of variables in data table
    N_Variables = length(Variables); %how many variables in table
    for j = 1:N_Variables
        Variable = char(Variables(j)); %get specific variable name
        Variable_ind = find(strcmp(InstrumentVariables.VarNameSpecific,Variable)); %get index for this variable within InstrumentVariables
        VarNameGeneral = char(InstrumentVariables.VarNameGeneral(Variable_ind)); %get the general variable name for this variable
        Units = char(InstrumentVariables.Units(Variable_ind)); %get the units for this variable
        CalibrationFactor = InstrumentVariables.CalibrationFactor(Variable_ind); %get calibration factor for this variable
        CalibrationIntercept = InstrumentVariables.CalibrationIntercept(Variable_ind); %get calibration factor for this variable
        InstrumentType = char(InstrumentVariables.InstrumentType(Variable_ind)); %get the instrument type for this variable
        Instrument = char(InstrumentVariables.Instrument(Variable_ind)); %get the instrument name for this variable
        Instrument_ind = find(strcmp(InstrumentMetadata.Instrument,Instrument)); %get indices in InstrumentMetadata associated with this instrument
        Intervals_ind = intersect(find(InstrumentMetadata.StartTime>=tableTimes.StartTime(1)),find(InstrumentMetadata.EndTime<=tableTimes.EndTime(end))); %get indices of InstrumentMetadata that occur within table interval
        InstrumentIntervals_ind = intersect(Instrument_ind,Intervals_ind); %get combined indices for instrument within time interval
        N_InstrumentIntervals = length(InstrumentIntervals_ind); %how many intervals in table
        for k = 1:N_InstrumentIntervals
            IntervalStartTime = InstrumentMetadata.StartTime(InstrumentIntervals_ind(k)); %get start time of instrument interval
            IntervalEndTime = InstrumentMetadata.EndTime(InstrumentIntervals_ind(k)); %get end time of instrument interval
            tableInd = find(isbetween(tableTimestamp,IntervalStartTime,IntervalEndTime)); %get indices of table within interval

            %get index of structured array for inserting values
            ind_struct = find(Instrument_ind == InstrumentIntervals_ind(k));

            %times and record numbers are general and not specific to the variable.  if these already exist, overwrite (should be same)
            RawData.(InstrumentType).(Instrument)(ind_struct).t.raw = tableTimestamp(tableInd); %assign times
            RawData.(InstrumentType).(Instrument)(ind_struct).t.Record = tableRecords(tableInd); %assign record numbers

            %add values to structered array at this index
            RawData.(InstrumentType).(Instrument)(ind_struct).(VarNameGeneral).Units = Units; %assign raw values for variable
            RawData.(InstrumentType).(Instrument)(ind_struct).(VarNameGeneral).CalibrationFactor = CalibrationFactor; %assign calibration factor for variable
            RawData.(InstrumentType).(Instrument)(ind_struct).(VarNameGeneral).CalibrationIntercept = CalibrationIntercept; %assign calibration intercept for variable
            RawData.(InstrumentType).(Instrument)(ind_struct).(VarNameGeneral).raw = tableVariables.(Variable)(tableInd); %assign raw values for variable
        end
    end
end

%% reorder fields so that variable information is first 
InstrumentTypes = fieldnames(RawData); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type
    Instruments = fieldnames(RawData.(InstrumentType));
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        Instrument = Instruments{j}; %get current instrument
        N_Fields = length(fieldnames(RawData.(InstrumentType).(Instrument)));
        RawData.(InstrumentType).(Instrument) = orderfields(RawData.(InstrumentType).(Instrument), [15:N_Fields,1:14]);
    end
end