% Dependencies: NONE
% Used by: GetVelocityProfiles

function [BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(StartTimes, EndTimes, TimeInterval)

%information about intervals
N_Intervals = length(StartTimes);

BlockStartTimes = [];
BlockEndTimes = [];
for i = 1:N_Intervals
    StartTimesBlock = StartTimes(i):TimeInterval:(EndTimes(i)-TimeInterval);
    BlockStartTimes = [BlockStartTimes; StartTimesBlock']; %starting times of blocks
    EndTimesBlock = (StartTimes(i)+TimeInterval):TimeInterval:EndTimes(i);
    BlockEndTimes = [BlockEndTimes; EndTimesBlock']; %ending times of blocks
end