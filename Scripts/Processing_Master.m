%% Master script for reading in raw logger data and producing useful processed data files
% Dependencies: 'CellStrFind', 'CombineIntervals', 'ExtractVariableTimeInterval',
% 'GetDataTypes', 'GetRawData', 'GetTimes', 'InterpolateData', 'IntervalsStartsWithin',
% 'IntervalsWithin', 'MaxOverlappingInterval', 'ParseBSNE', 'ParseInstrumentMetadata',
% 'ParseInstrumentVariables', 'ParseLoggerTables', 'ParseLoggerTimes', 
% 'ParseTimestamp', 'ProcessBSNEs', 'ProcessData', 'qzCalc', 'ReadTable', 
% 'RoundTimeMin'

%% 0. Initialize
clear all;

%% 1. Assign folder and files with data and metadata
folder_Metadata = '../Metadata/'; %folder with metadata
folder_LoggerRawData = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Raw/'; %folder with logger raw data tables
folder_DataOutput = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
file_LoggerTables = 'LoggerTables.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata.xlsx'; %file with logger tables
file_WeightBSNE = 'WeightBSNE.xlsx'; %file with BSNE weights and times
Metadata_Path = strcat(folder_DataOutput,'Metadata_Oceano'); %get path to saving meta data
RawData_Path = strcat(folder_DataOutput,'RawData_Oceano'); %get path to saving raw data
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_Oceano'); %get path to saving interpolated data
ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_Oceano'); %get path to saving processed data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_Metadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_Metadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_Metadata,file_InstrumentVariables); %extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_Metadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file
WeightBSNE = ParseBSNE(folder_Metadata,file_WeightBSNE); %extract info about BSNE weights from .xlsx spreadsheet
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables','InstrumentMetadata','WeightBSNE','-v7.3');

%% 3. Import and aggregate data from logger files
RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
%save(RawData_Path,'RawData','-v7.3'); %save raw data

%% 4. Perform interpolation on raw data
%load(RawData_Path); %load raw data
InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors
%save(InterpolatedData_Path,'InterpolatedData','-v7.3');

%% 5. Process BSNE profiles and interpolated data
%load(InterpolatedData_Path); %load interpolated data
ProcessedData = ProcessData(InterpolatedData, WeightBSNE, InstrumentMetadata); %process data
save(ProcessedData_Path,'ProcessedData','-v7.3'); %save processed data