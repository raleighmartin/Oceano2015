%% Function to read in metadata from 'file_InstrumentMetadata' (a .xlsx spreadsheet) that
% describes information about continuous time intervals of data collection for each
% instrument.  Information includes start times, end times, and spatial information.
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: RoundTimeMin
% Used by: Processing_Master

function InstrumentMetadata = ParseInstrumentMetadata(folder_InstrumentMetadata,file_InstrumentMetadata)

%% Get file location
filePath = strcat(folder_InstrumentMetadata,file_InstrumentMetadata);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Site = {raw(:,1)};
Date = round(cell2mat(raw(:,2))); %round this to get precise time
StartTime = cell2mat(raw(:,3));
EndTime = cell2mat(raw(:,4));
InstrumentType = {raw(:,5)};
Instrument = {raw(:,6)};
StartHeight_m = cell2mat(raw(:,7));
EndHeight_m = cell2mat(raw(:,8));
HeightErr_m = cell2mat(raw(:,9));
HeightRef = {cellfun(@num2str,raw(:,10),'UniformOutput',0)}; %convert to strings
Longitudinal_m = cell2mat(raw(:,11));
Spanwise_m = cell2mat(raw(:,12));
AngleErr_deg = cell2mat(raw(:,13));
ErrorCode = cell2mat(raw(:,14));

%% Convert dates to MATLAB datetime format
StartTime = datetime(Date+StartTime,'ConvertFrom','excel');
EndTime = datetime(Date+EndTime,'ConvertFrom','excel');
Date = datetime(Date,'ConvertFrom','excel');

%% Correct for roundoff error (assuming time resolution of minutes)
StartTime = RoundTimeMin(StartTime);
EndTime = RoundTimeMin(EndTime);

%% Create structured array with spreadsheet information
InstrumentMetadata = struct(...
    'Site',Site,...
    'Date',Date,...
    'StartTime',StartTime,...
    'EndTime',EndTime,...
    'InstrumentType',InstrumentType,...
    'Instrument',Instrument,...
    'StartHeight_m',StartHeight_m,...
    'EndHeight_m',EndHeight_m,...
    'HeightErr_m',HeightErr_m,...
    'HeightRef',HeightRef,...
    'Longitudinal_m',Longitudinal_m,...
    'Spanwise_m',Spanwise_m,...
    'AngleErr_deg',AngleErr_deg,...
    'ErrorCode',ErrorCode);
end