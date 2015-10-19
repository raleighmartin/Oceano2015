%% Master script for reading in raw logger data and producing useful interpolated data files, but without interpretation

%% 0. Initialize
clearvars;

%% 1. Assign folder and files with data and metadata
folder_SiteMetadata = '../Metadata/'; %folder with site metadata
folder_LoggerRawData = '../../../Google Drive/Data/AeolianFieldwork/Raw/Datalogger/'; %folder with logger raw data tables
folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Scripts = '../../AeolianFieldworkAnalysis/Scripts/'; %folder with general scripts for data analysis
path(path,folder_Scripts); %add this path to locations for running functions

%site-specific information
file_LoggerTables = 'LoggerTables_OceanoCalibration.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes_OceanoCalibration.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables_OceanoCalibration.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata_OceanoCalibration.xlsx'; %file with logger tables
Metadata_Path = strcat(folder_DataOutput,'Metadata_OceanoCalibration'); %get path to saving meta data
RawData_Path = strcat(folder_DataOutput,'RawData_OceanoCalibration'); %get path to saving raw data
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_OceanoCalibration'); %get path to saving interpolated data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_SiteMetadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_SiteMetadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_SiteMetadata,file_InstrumentVariables); %!!extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_SiteMetadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file

%% 3. Save all metadata together
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables',...
   'InstrumentMetadata','-v7.3');

% %% 4. Import and aggregate data from logger files
% RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
% save(RawData_Path,'RawData','-v7.3'); %save raw data

%% 5. Perform interpolation on raw data
load(RawData_Path); %load raw data if necessary
InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors
save(InterpolatedData_Path,'InterpolatedData','-v7.3'); %save data

%% restore function path to default value
restoredefaultpath;