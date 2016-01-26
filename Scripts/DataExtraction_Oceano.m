%% Master script for reading in raw logger data and producing useful interpolated data files, but without interpretation

%% 0. Initialize
clearvars;

%% 1. Assign folder and files with data and metadata
folder_SiteMetadata = '../Metadata/'; %folder with site metadata
folder_GeneralMetadata = '../../AeolianFieldworkAnalysis/Metadata/'; %folder with general metadata
folder_LoggerRawData = '../../../Google Drive/Data/AeolianFieldwork/Raw/Datalogger/'; %folder with logger raw data tables
folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_GrainSize = '../../../Google Drive/Data/AeolianFieldwork/Raw/GrainSize/'; %folder with grain size data
file_InstrumentCalibration = 'InstrumentCalibration.xlsx'; %file with instrument calibration values
folder_Functions = '../../AeolianFieldworkAnalysis/Scripts/Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%site-specific information
file_LoggerTables = 'LoggerTables_Oceano.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes_Oceano.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables_Oceano.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata_Oceano.xlsx'; %file with logger tables
file_WeightBSNE = 'WeightBSNE_Oceano.xlsx'; %file with BSNE weights and times
file_GrainSizeMetadata = 'GrainSizeMetadata_Oceano.xlsx'; %file with grain size metadata information
Metadata_Path = strcat(folder_DataOutput,'Metadata_Oceano'); %get path to saving meta data
RawData_Path = strcat(folder_DataOutput,'RawData_Oceano'); %get path to saving raw data
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_Oceano'); %get path to saving interpolated data
GrainSize_Path = strcat(folder_DataOutput,'GrainSize_Oceano'); % get path to saving grain size data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_SiteMetadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_SiteMetadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_SiteMetadata,file_InstrumentVariables); %!!extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_SiteMetadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file
InstrumentCalibration = ParseInstrumentCalibration(folder_GeneralMetadata,file_InstrumentCalibration); %!!extract info about instrument calibrations from .xlsx file
WeightBSNE = ParseBSNE(folder_SiteMetadata,file_WeightBSNE); %extract info about BSNE weights from .xlsx spreadsheet
[GrainSizeMetadata_Surface, GrainSizeMetadata_BSNE] = ParseGrainSizeMetadata(folder_SiteMetadata,file_GrainSizeMetadata); %extract info about instruments from .xlsx file

%% 3. Import and aggregate data from grain size files
GrainSize_Surface = GetGrainSize(GrainSizeMetadata_Surface,folder_GrainSize);
GrainSize_BSNE = GetGrainSize(GrainSizeMetadata_BSNE,folder_GrainSize);

%% 4. Save all metadata together
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables',...
   'InstrumentMetadata','InstrumentCalibration','WeightBSNE',...
   'GrainSize_Surface','GrainSize_BSNE',...
   'GrainSizeMetadata_Surface','GrainSizeMetadata_BSNE','-v7.3');

%% 5. Import and aggregate data from logger files
RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
save(RawData_Path,'RawData','-v7.3'); %save raw data

%% 6. Perform interpolation on raw data
%load(RawData_Path); %load raw data if necessary
InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors

%% 7. Apply calibration factor to data based on unique identifier
InterpolatedData = ParseCalibrationFactor(InterpolatedData,InstrumentVariables,InstrumentCalibration);
save(InterpolatedData_Path,'InterpolatedData','-v7.3'); %save data

%% restore function path to default value
restoredefaultpath;