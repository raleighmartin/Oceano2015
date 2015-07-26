%% Function to read in metadata from 'file_LoggerTables' (a .xlsx spreadsheet) that
% describes information about each raw logger data file
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: NONE
% Used by: Processing_Master

function LoggerTables = ParseLoggerTables(folder_LoggerMetadata,file_LoggerTables)

%% Get file location
filePath = strcat(folder_LoggerMetadata,file_LoggerTables);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Site = {raw(:,1)};
Date = cell2mat(raw(:,2));
Logger = {raw(:,3)};
Table = {raw(:,4)};
Variables = {raw(:,5)};
ErrorCode = cell2mat(raw(:,6));

%% Convert dates to MATLAB datetime format
Date = datetime(Date,'ConvertFrom','excel');

%% Create structured array with spreadsheet information
LoggerTables = struct(...
    'Site',Site,...
    'Date',Date,...
    'Logger',Logger,...
    'Table',Table,...
    'Variables',Variables,...
    'ErrorCode',ErrorCode);

end