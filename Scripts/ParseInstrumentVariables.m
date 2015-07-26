%% Function to read in metadata from 'file_InstrumentVariables' (a .xlsx spreadsheet) that
% describes information about the variables associated with each instrument
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: NONE
% Used by: Processing_Master

function InstrumentVariables = ParseInstrumentVariables(folder_InstrumentMetadata,file_InstrumentVariables)
    
%% Get file location
filePath = strcat(folder_InstrumentMetadata,file_InstrumentVariables);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
InstrumentType = {raw(:,1)};
Instrument = {raw(:,2)};
VarNameSpecific = {raw(:,3)};
VarNameGeneral = {raw(:,4)};
DataType = {raw(:,5)};
Units = {raw(:,6)};
CalibrationFactor = cell2mat(raw(:,7));
CalibrationIntercept = cell2mat(raw(:,8));

%% Create structured array with spreadsheet information
InstrumentVariables = struct(...
    'InstrumentType',InstrumentType,...
    'Instrument',Instrument,...
    'VarNameSpecific',VarNameSpecific,...
    'VarNameGeneral',VarNameGeneral,...
    'DataType',DataType,...
    'Units',Units,...
    'CalibrationFactor',CalibrationFactor,...
    'CalibrationIntercept',CalibrationIntercept);

end