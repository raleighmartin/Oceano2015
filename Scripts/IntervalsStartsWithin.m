%% Determine 'StartTimesWithin' and 'EndTimesWithin' of the intervals
% given by 'StartTimesAll' and 'EndTimesAll' that start within the interval
% from 'StartTime' to 'EndTime'. Also return 'indices_StartsWithin' of
% original 'StartTimesAll' list that end up in 'StartTimesWithin'
% Dependencies: NONE
% Used by: ProcessBSNEs

function [StartTimesWithin, EndTimesWithin, indices_StartsWithin] = IntervalsStartsWithin(StartTime,EndTime,StartTimesAll,EndTimesAll)

indices_StartsWithin = find(isbetween(StartTimesAll,StartTime,EndTime));
StartTimesWithin = StartTimesAll(indices_StartsWithin);
EndTimesWithin = EndTimesAll(indices_StartsWithin);

end