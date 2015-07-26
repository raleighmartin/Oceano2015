%% function to convert datetimes contained in 'OldTimes' into datetimes rounded to the nearest minute
% Dependencies: NONE
% Used by: ParseLoggerTimes, ParseInstrumentMetadata, ParseBSNE

function [NewTimes] = RoundTimeMin(OldTimes)

N_Times = length(OldTimes);
[Y,M,D,H,MN,S] = datevec(OldTimes);
for i = 1:N_Times
    if S(i)>30;
        MN(i)=MN(i)+1;
    end
end
S = zeros(N_Times,1); %set seconds to zero

NewTimes = datetime([Y,M,D,H,MN,S]);

end