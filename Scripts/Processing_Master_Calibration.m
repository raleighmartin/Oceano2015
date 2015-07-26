%% Master script for reading in raw logger data and producing useful processed data files
% specifically for analysis of calibration
% Dependencies: 'GetDataTypes', 'GetRawData', 'GetTimes',
% 'InterpolateData', 'ParseInstrumentMetadata', 'ParseInstrumentVariables',
% 'ParseLoggerTables', 'ParseLoggerTimes', 'ParseTimestamp', 'ReadTable',
% 'RoundTimeMin'

%% 0. Initialize
clear all;

%% 1. Assign folder and files with data and metadata
folder_Metadata = '../Metadata/'; %folder with metadata
folder_LoggerRawData = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Raw/'; %folder with logger raw data tables
folder_DataOutput = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
file_LoggerTables = 'LoggerTables_OceanoCalibration.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes_OceanoCalibration.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata_OceanoCalibration.xlsx'; %file with logger tables
Metadata_Path = strcat(folder_DataOutput,'Metadata_OceanoCalibration'); %get path to saving meta data
RawData_Path = strcat(folder_DataOutput,'RawData_OceanoCalibration'); %get path to saving raw data
ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_OceanoCalibration'); %get path to saving processed data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_Metadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_Metadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_Metadata,file_InstrumentVariables); %extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_Metadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables','InstrumentMetadata','-v7.3');

%% 3. Import and aggregate data from logger files
RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
save(RawData_Path,'RawData','-v7.3'); %save raw data

%% 4. Perform interpolation on raw data, save as processed data
%load(RawData_Path); %load raw data
InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors
ProcessedData = InterpolatedData; %Processed Data is simply interpolated data
save(ProcessedData_Path,'ProcessedData','-v7.3');