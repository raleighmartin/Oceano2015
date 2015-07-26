%% function to determine the maximum overlapping interval of the given
% intervals, beginning with an initial start and end time
% Dependencies: IntervalsWithin
% Used by: CombineIntervals

function [StartTimeMax, EndTimeMax] = MaxOverlappingInterval(StartTime, EndTime, StartTimesAll, EndTimesAll)

[StartTimesWithin, EndTimesWithin] = IntervalsWithin(StartTime,EndTime,StartTimesAll,EndTimesAll); %find all intervals contained within these start and end times
StartTimeMax = min(StartTimesWithin); %new start time is the earliest of start times of these contained intervals
EndTimeMax = max(EndTimesWithin); %new end time is the latest of end times of these contained intervals

%continue this procedure recursive process captures no further reduction in number of start times 
if length(StartTimesWithin)~=length(StartTimesAll)
    [StartTimeMax, EndTimeMax] = MaxOverlappingInterval(StartTimeMax, EndTimeMax, StartTimesWithin, EndTimesWithin);
end

end