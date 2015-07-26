%% Determine 'StartTimesWithin' and 'EndTimesWithin' of the intervals
% given by 'StartTimesAll' and 'EndTimesAll' that start or end within the
% interval from 'StartTime' to 'EndTime'. Also return 'indices_Within' of
% original 'StartTimesAll' list that end up in 'StartTimesWithin'
% Dependencies: NONE
% Used by: MaxOverlappingInterval, GetVelocityProfiles

function [StartTimesWithin, EndTimesWithin, indices_Within] = IntervalsWithin(StartTime,EndTime,StartTimesAll,EndTimesAll)

indices_StartsWithin = find(isbetween(StartTimesAll,StartTime,EndTime));
indices_EndsWithin = find(isbetween(EndTimesAll,StartTime,EndTime));
indices_Within = union(indices_StartsWithin,indices_EndsWithin);
StartTimesWithin = StartTimesAll(indices_Within);
EndTimesWithin = EndTimesAll(indices_Within);

end