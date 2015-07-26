%% Given 'StartTimes' and 'EndTimes', combine these into minimal set of
% intervals with 'CombinedStartTimes' and 'CombinedEndTimes' for which
% there are overlapping times from instruments
% Dependencies: MaxOverlappingInterval
% Used by: TimeIntervalUnion

function [CombinedStartTimes, CombinedEndTimes] = CombineIntervals(StartTimes, EndTimes)

%get start and end time for first interval
index_FirstStartTime = find(StartTimes==min(StartTimes), 1); %get index for interval corresponding to first start time
FirstStartTime = StartTimes(index_FirstStartTime); %get first start time
FirstEndTime = EndTimes(index_FirstStartTime); %get end time of the interval with first start time

%get max overlapping interval of start and end times
[CombinedStartTimes, CombinedEndTimes] = MaxOverlappingInterval(FirstStartTime, FirstEndTime, StartTimes, EndTimes);

%get indices of all remaining start times after last combined end time
index_LaterStartTimes = find(StartTimes>CombinedEndTimes(end));

%continue if there are additional later start times
while isempty(index_LaterStartTimes)==0
    NextStartTime = min(StartTimes(index_LaterStartTimes)); %get earliest of remaining start times
    index_NextStartTime = find(StartTimes==NextStartTime, 1); %get index for interval corresponding to this start time
    NextEndTime = EndTimes(index_NextStartTime); %get end time for this interval
    [NextCombinedStartTime, NextCombinedEndTime] = MaxOverlappingInterval(NextStartTime, NextEndTime, StartTimes, EndTimes); %get max overlapping interval for this start and end time
    CombinedStartTimes = [CombinedStartTimes; NextCombinedStartTime]; %append combined start time to list
    CombinedEndTimes = [CombinedEndTimes; NextCombinedEndTime]; %append combined end time to list
    index_LaterStartTimes = find(StartTimes>CombinedEndTimes(end)); %get indices of all remaining start times after last combined end time
end