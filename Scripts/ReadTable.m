%% function to parse logger data file at 'tableFilepath' for the variables
% 'tableVariableList', the types 'tableTypeList' and units 'tableUnitsList'
% at the times 'tableTimes'
% at the specific 'site' on the specific 'date'
% Output 'timestamp_allinteravls', 'record_allintervals' as lists of all timestamps and associated record numbers
% Output 'variables' as structured array of variables, with each element as a different variable from the data table
% Dependencies: ParseTimestamp
% Used by: GetRawData

function [timestamp_allintervals, record_allintervals, variables] = ReadTable(tableFilepath,tableVariableList,tableTypeList,tableTimes)

%% Initialize input
delimiter = ',';
startRow = 5;
endRow = inf;

%% Set format string
formatSpec = '%s%f'; %initialize format string
N_variables = length(tableVariableList); %how many variables?
for i = 1:N_variables
    if strcmp(tableTypeList(i),'float') %floats
        formatSpec = strcat(formatSpec,'%f');
    elseif strcmp(tableTypeList(i),'int') %integers
        formatSpec = strcat(formatSpec,'%d');
    else %set other as strings
        formatSpec = strcat(formatSpec,'%s');
    end
end
formatSpec = strcat(formatSpec,'%[^\n\r]');

%% Open the text file
fileID = fopen(tableFilepath,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'TreatAsEmpty','"NAN"','HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'TreatAsEmpty','"NAN"','HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Get time and record information
timestamp = ParseTimestamp(dataArray{1});
record = dataArray{2};

%% Adjust times based on recorded delay
dt = mode(diff(timestamp)); %get dt for calculating adjustments
timestamp_allintervals = []; %initialize list of all timestamps
ind_timestamp_allintervals = []; %initialize list of indices for all timestamps
for i = 1:tableTimes.N_Intervals
    
    %determine delay to apply to apply to timestamps
    Delay_duration = duration(0,0,tableTimes.Delay_Seconds(i)); %compute duration of delay
    Delay_steps = round(Delay_duration/dt); %how many timesteps for this delay?
    Delay_duration = Delay_steps*dt; %convert back to duration to ensure that value is exact
    
    %compute adjusted timestamps
    timestamp_adjusted = timestamp+Delay_duration;
    
    %get indices of adjusted timestamps within range
    ind_timestamp_interval = find((timestamp_adjusted>=tableTimes.StartTime(i))&...
        (timestamp_adjusted<tableTimes.EndTime(i)));
    
    %get adjusted timestamps within range
    timestamp_interval = timestamp_adjusted(ind_timestamp_interval);
    
    %add to lists of timestamps and timestamp indices
    timestamp_allintervals = [timestamp_allintervals; timestamp_interval];
    ind_timestamp_allintervals = [ind_timestamp_allintervals; ind_timestamp_interval];
end

%% get record numbers for all intervals
record_allintervals = record(ind_timestamp_allintervals);

%% create structured array with variable values

variables = struct; %initialize structured array

% add variables to structured array
for i = 1:N_variables
    variableName = tableVariableList{i}; %get variable name
    values = dataArray{i+2}; %get values for variable
    values_allintervals = values(ind_timestamp_allintervals); %limit values to indices for timestamps
    
    %assign variable to structured array
    variables.(variableName) = values_allintervals;
end