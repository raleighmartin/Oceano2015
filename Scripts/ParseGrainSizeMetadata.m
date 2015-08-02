%% Function to read in data from 'file_GrainSizeMetadata' (a .xlsx spreadsheet) that
% gives basic info about grain size samples for analysis
% Function outputs two structured arrays with file information about grain
% size samples for surface sand and BSNes
% Dependencies: NONE, RoundTimeMin
% Used by: Processing_Master

function [GrainSizeMetadata_Surface, GrainSizeMetadata_BSNE] = ...
    ParseGrainSizeMetadata(folder_GrainSizeMetadata,file_GrainSizeMetadata)

%% Get file location
filePath = strcat(folder_GrainSizeMetadata,file_GrainSizeMetadata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Surface grain size metadata %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, raw] = xlsread(filePath,'SurfaceSamples');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Filename = raw(:,1);
Site = raw(:,2);
Date = round(cell2mat(raw(:,3))); %round this to get precise time
Location = raw(:,4);
CollectionTime = cell2mat(raw(:,5));

%% Convert dates to MATLAB datetime format
CollectionTime = datetime(Date+CollectionTime,'ConvertFrom','excel');
Date = datetime(Date,'ConvertFrom','excel');

%% Correct for roundoff error (assuming time resolution of minutes)
CollectionTime = RoundTimeMin(CollectionTime);

%% Assign values to structured array
N_Samples = length(Filename);

% for i = 1:N_Samples
%     GrainSizeMetadata_Surface(i).Filename = Filename{i};
%     GrainSizeMetadata_Surface(i).Site = Site{i};
%     GrainSizeMetadata_Surface(i).Date = Date(i);
%     GrainSizeMetadata_Surface(i).Location = Location{i};
%     GrainSizeMetadata_Surface(i).CollectionTime = CollectionTime(i);
% end

GrainSizeMetadata_Surface = struct(...
    'Filename',cellstr(Filename),...
    'Site',cellstr(Site),...
    'Date',num2cell(Date),...
    'Location',cellstr(Location),...
    'CollectionTime',num2cell(CollectionTime));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import BSNE grain size metadata %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, raw] = xlsread(filePath,'BSNE');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
Filename = raw(:,1);
Site = raw(:,2);
Date = round(cell2mat(raw(:,3))); %round this to get precise time
NameBSNE = raw(:,4);
StartTime = cell2mat(raw(:,5));
EndTime = cell2mat(raw(:,6));
Notes = raw(:,7);

%% Convert dates to MATLAB datetime format
StartTime = datetime(Date+StartTime,'ConvertFrom','excel');
EndTime = datetime(Date+EndTime,'ConvertFrom','excel');
Date = datetime(Date,'ConvertFrom','excel');

%% Correct for roundoff error (assuming time resolution of minutes)
StartTime = RoundTimeMin(StartTime);
EndTime = RoundTimeMin(EndTime);

% %% Assign values to structured array
% N_Samples = length(Filename);
% 
% for i = 1:N_Samples
%     GrainSizeMetadata_BSNE(i).Filename = Filename;
%     GrainSizeMetadata_BSNE(i).Site = Site;
%     GrainSizeMetadata_BSNE(i).Date = Date;
%     GrainSizeMetadata_BSNE(i).NameBSNE = NameBSNE;
%     GrainSizeMetadata_BSNE(i).StartTime = StartTime;
%     GrainSizeMetadata_BSNE(i).EndTime = EndTime;
%     GrainSizeMetadata_BSNE(i).Notes = Notes;
% end

GrainSizeMetadata_BSNE = struct(...
    'Filename',cellstr(Filename),...
    'Site',cellstr(Site),...
    'Date',num2cell(Date),...
    'NameBSNE',cellstr(NameBSNE),...
    'StartTime',num2cell(StartTime),...
    'EndTime',num2cell(EndTime),...
    'Notes',cellstr(Notes));