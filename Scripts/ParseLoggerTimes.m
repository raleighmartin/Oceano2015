%% Function to read in metadata from 'file_LoggerTimes' (a .xlsx spreadsheet) that
% describes information about times to read in from logger raw data files and time offsets
% among loggers determined from timing pulses. 
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: RoundTimeMin
% Used by: Processing_Master

function LoggerTimes = ParseLoggerTimes(folder_LoggerMetadata,file_LoggerTimes)

%% Get file location
filePath = strcat(folder_LoggerMetadata,file_LoggerTimes);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Site = {raw(:,1)};
Date = round(cell2mat(raw(:,2))); %round this to nearest day
Logger = {raw(:,3)};
StartTime = cell2mat(raw(:,4));
EndTime = cell2mat(raw(:,5));
Delay_Seconds = cell2mat(raw(:,6));

%% Convert times to MATLAB datetime format
StartTime = datetime(Date+StartTime,'ConvertFrom','excel');
EndTime = datetime(Date+EndTime,'ConvertFrom','excel');
Date = datetime(Date,'ConvertFrom','excel');

%% Correct for roundoff error (assuming time resolution of minutes)
StartTime = RoundTimeMin(StartTime);
EndTime = RoundTimeMin(EndTime);

%% Create structured array with spreadsheet information
LoggerTimes = struct(...
    'Site',Site,...
    'Date',Date,...
    'Logger',Logger,...
    'StartTime',StartTime,...
    'EndTime',EndTime,...
    'Delay_Seconds',Delay_Seconds);

end