% Script for taking interpolated data and making useful basic calculations

%% 0. Initialize
clearvars;

%% 1. Assign folder and files with data
folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Scripts = '../../AeolianFieldworkAnalysis/Scripts/'; %folder with general scripts for data analysis
path(path,folder_Scripts); %add this path to locations for running functions

%load metadata and interpolated data, create path for processed data
Metadata_Path = strcat(folder_DataOutput,'Metadata_OceanoCalibration'); %get path to metadata
load(Metadata_Path);
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_OceanoCalibration'); %get path to interpolated data
load(InterpolatedData_Path);
ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_OceanoCalibration'); %create path to processed data

%% 1. Process instrument heights
ProcessedData = ProcessInstrumentHeights(InterpolatedData);

%% 2. Move extraneous metadata fields out of file
ProcessedData = MoveExtraneousMetadataFields(ProcessedData);

%% 3. Save processed data
save(ProcessedData_Path,'ProcessedData','-v7.3'); %save data

%% Restore function path to default value
restoredefaultpath;