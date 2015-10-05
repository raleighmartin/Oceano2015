%% Master script for reading in raw logger data and producing useful processed data files
% Dependencies: 'CellStrFind', 'CombineIntervals', 'ExtractVariableTimeInterval',
% 'GetDataTypes', 'GetRawData', 'GetTimes', 'InterpolateData', 'IntervalsStartsWithin',
% 'IntervalsWithin', 'MaxOverlappingInterval', 'ParseBSNE', 'ParseInstrumentMetadata',
% 'ParseInstrumentVariables', 'ParseLoggerTables', 'ParseLoggerTimes', 
% 'ParseTimestamp', 'ParseGrainSizeMetadata','ProcessBSNEs', 'ProcessData',
% 'GetGrainSize', 'qzCalc', 'ReadTable','RoundTimeMin','gsimport',
% 'qz_profilefit', 't_95'

%% 0. Initialize
clear all;

%% 1. Assign folder and files with data and metadata
folder_Metadata = '../Metadata/'; %folder with metadata
folder_LoggerRawData = '../../../Google Drive/Data/AeolianFieldwork/Raw/Datalogger/'; %folder with logger raw data tables
folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_GrainSize = '../../../Google Drive/Data/AeolianFieldwork/Raw/GrainSize/'; %folder with grain size data
file_LoggerTables = 'LoggerTables.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata.xlsx'; %file with logger tables
file_WeightBSNE = 'WeightBSNE.xlsx'; %file with BSNE weights and times
file_GrainSizeMetadata = 'GrainSizeMetadata.xlsx'; %file with grain size metadata information
Metadata_Path = strcat(folder_DataOutput,'Metadata_Oceano'); %get path to saving meta data
RawData_Path = strcat(folder_DataOutput,'RawData_Oceano'); %get path to saving raw data
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_Oceano'); %get path to saving interpolated data
ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_Oceano'); %get path to saving processed data
GrainSize_Path = strcat(folder_DataOutput,'GrainSize_Oceano'); % get path to saving grain size data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_Metadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_Metadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_Metadata,file_InstrumentVariables); %extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_Metadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file
WeightBSNE = ParseBSNE(folder_Metadata,file_WeightBSNE); %extract info about BSNE weights from .xlsx spreadsheet
[GrainSizeMetadata_Surface, GrainSizeMetadata_BSNE] = ParseGrainSizeMetadata(folder_Metadata,file_GrainSizeMetadata); %extract info about instruments from .xlsx file
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables',...
    'InstrumentMetadata','WeightBSNE','GrainSizeMetadata_Surface',...
    'GrainSizeMetadata_BSNE','-v7.3');
% 
% %% 3. Import and aggregate data from grain size files
% GrainSize_Surface = GetGrainSize(GrainSizeMetadata_Surface,folder_GrainSize);
% GrainSize_BSNE = GetGrainSize(GrainSizeMetadata_BSNE,folder_GrainSize);
% save(GrainSize_Path,'GrainSize_Surface','GrainSize_BSNE');
% 
% %% 4. Import and aggregate data from logger files
% RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
% save(RawData_Path,'RawData','-v7.3'); %save raw data
% 
% %% 5. Perform interpolation on raw data
% %load(RawData_Path); %load raw data
% InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors
% save(InterpolatedData_Path,'InterpolatedData','-v7.3');
% 
% % 6. Process BSNE profiles and interpolated data
% %load(InterpolatedData_Path); %load interpolated data
% ProcessedData = ProcessData(InterpolatedData, WeightBSNE, InstrumentMetadata, GrainSize_BSNE, GrainSize_Surface); %process data
% save(ProcessedData_Path,'ProcessedData','-v7.3'); %save processed data