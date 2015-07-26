%% function to extract data for a variable given a time interval
% INPUTS:
% 'InstrumentData' is data from an individual instrument
% 'StartTime' and 'EndTime', given in datetime format, describe the interval of interest
% 'Variable' is the name of the variable for the given instrument to be examined
% 'Subvariable' is an optional argument describing which values to use. Default is "int".
% 'TimeSubvariable' is an optional argument specifying which times to use.  Default is "int".
% variable to find the values.  If none given, assume that we are not usin
% OUTPUTS:
% 'IntervalVariable' are the resulting values of the variable over the time interval
% 'IntervalTimes' gives the times associated with the observations
% 'IntervalNumbers' gives the numbers of the intervals within the list
% 'IntervalIndices' gives the specific index of the value within its interval
% DEPENDENCIES: NONE
% USED BY: ProcessWenglors, GetVelocityProfiles, StressFluxBSNE,
% WindFluxAnalysis

function [IntervalValues, IntervalTimes, IntervalNumbers, IntervalIndices] = ...
    ExtractVariableTimeInterval(InstrumentData,StartTime,EndTime,Variable,Subvariable,TimeSubvariable)

%set subvariable values as 'int' if none given
if nargin<5
    Subvariable = 'int';
end

if nargin<6
    TimeSubvariable = 'int';
end

%get number of time intervals associated with instrument data
N_Intervals = length(InstrumentData);

%initialize lists of variable values, times, and indices
IntervalValues = [];
IntervalTimes = [];
IntervalNumbers = [];
IntervalIndices = {}; %cell array, with each element corresponding to one interval number 

%go through each interval
for i = 1:N_Intervals
    
    %get indices of times for interval between given start and end times
    indices_StartTime = find(InstrumentData(i).t.(TimeSubvariable)>=StartTime);
    indices_EndTime = find(InstrumentData(i).t.(TimeSubvariable)<=EndTime);
    indices_Interval = intersect(indices_StartTime,indices_EndTime);
    
    %if any times exist, extract values
    if ~isempty(indices_Interval)
        IntervalNumbers = [IntervalNumbers; i]; %add to list of interval numbers
        IntervalIndices{length(IntervalIndices)+1} = indices_Interval; %add to list of indices
        IntervalValues = [IntervalValues; InstrumentData(i).(Variable).(Subvariable)(indices_Interval)]; %add to list of interval values
        IntervalTimes = [IntervalTimes; InstrumentData(i).t.(TimeSubvariable)(indices_Interval)]; %add to list of interval times
    end
end
