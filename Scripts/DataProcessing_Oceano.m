% Script for taking interpolated data and making useful basic calculations
% fit profile to BSNEs

%% 0. Initialize
clearvars;

%% 1. Assign folder and files with data
folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Functions = '../../AeolianFieldworkAnalysis/Scripts/Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%load metadata and interpolated data, create path for processed data
Metadata_Path = strcat(folder_DataOutput,'Metadata_Oceano'); %get path to metadata
load(Metadata_Path);
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_Oceano'); %get path to interpolated data
load(InterpolatedData_Path);
ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_Oceano'); %create path to processed data

%% 1. Process instrument heights
ProcessedData = ProcessInstrumentHeights(InterpolatedData);

%% 2. Process BSNE profiles
FluxBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE);
ProcessedData.FluxBSNE = FluxBSNE;

%% 3. Process Wenglors
[ProcessedWenglors, FluxWenglor] = ProcessWenglors(ProcessedData, FluxBSNE, InstrumentMetadata);
ProcessedData.Wenglor = ProcessedWenglors;
ProcessedData.FluxWenglor = FluxWenglor;

%% 4. Move extraneous metadata fields out of file
ProcessedData = MoveExtraneousMetadataFields(ProcessedData);

%% 5. Save processed data
save(ProcessedData_Path,'ProcessedData','-v7.3'); %save data

%% Restore function path to default value
restoredefaultpath;