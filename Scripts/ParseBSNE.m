%% Function to read in data from 'file_WeightBSNE' (a .xlsx spreadsheet) that
% gives measured weights of sediment collected in BSNE traps during prescribed intervals
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: NONE, RoundTimeMin
% Used by: Processing_Master

function WeightBSNE = ParseBSNE(folder_InstrumentMetadata,file_WeightBSNE)

%% Get file location
filePath = strcat(folder_InstrumentMetadata,file_WeightBSNE);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Site = {raw(:,1)};
Date = round(cell2mat(raw(:,2))); %round this to get precise time
StartTime = cell2mat(raw(:,3));
EndTime = cell2mat(raw(:,4));
NameBSNE = {raw(:,5)};
Weight_g = cell2mat(raw(:,6));
WeightErr_g = cell2mat(raw(:,7));
BottomHeight_cm = cell2mat(raw(:,8));
BottomHeightErr_cm = cell2mat(raw(:,9));
HeightBSNE_cm = cell2mat(raw(:,10));
WidthBSNE_cm = cell2mat(raw(:,11));
Longitudinal_m = cell2mat(raw(:,12));
Spanwise_m = cell2mat(raw(:,13));
ErrorCode = cell2mat(raw(:,14));
StartTime_GrainSize = cell2mat(raw(:,15));
EndTime_GrainSize = cell2mat(raw(:,16));

%% Convert dates to MATLAB datetime format
StartTime = datetime(Date+StartTime,'ConvertFrom','excel');
EndTime = datetime(Date+EndTime,'ConvertFrom','excel');
StartTime_GrainSize = datetime(Date+StartTime_GrainSize,'ConvertFrom','excel');
EndTime_GrainSize = datetime(Date+EndTime_GrainSize,'ConvertFrom','excel');
Date = datetime(Date,'ConvertFrom','excel');

%% Correct for roundoff error (assuming time resolution of minutes)
StartTime = RoundTimeMin(StartTime);
EndTime = RoundTimeMin(EndTime);
StartTime_GrainSize = RoundTimeMin(StartTime_GrainSize);
EndTime_GrainSize = RoundTimeMin(EndTime_GrainSize);

%% Create structured array with spreadsheet information
WeightBSNE = struct(...
    'Site',Site,...
    'Date',Date,...
    'StartTime',StartTime,...
    'EndTime',EndTime,...
    'NameBSNE',NameBSNE,...
    'Weight_g',Weight_g,...
    'WeightErr_g',WeightErr_g,...
    'BottomHeight_cm',BottomHeight_cm,...
    'BottomHeightErr_cm',BottomHeightErr_cm,...
    'HeightBSNE_cm',HeightBSNE_cm,...
    'WidthBSNE_cm',WidthBSNE_cm,...
    'Longitudinal_m',Longitudinal_m,...
    'Spanwise_m',Spanwise_m,...
    'ErrorCode',ErrorCode,...
    'StartTime_GrainSize',StartTime_GrainSize,...
    'EndTime_GrainSize',EndTime_GrainSize);