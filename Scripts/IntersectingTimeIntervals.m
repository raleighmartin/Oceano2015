%function to get intersecting time intervals from two sets of time
%intervals

function [StartTimesIntersecting, EndTimesIntersecting] = ...
    IntersectingTimeIntervals(StartTimes1, EndTimes1, StartTimes2, EndTimes2)

%initialize list of intersecting start and end times
StartTimesIntersecting = [];
EndTimesIntersecting = [];

%go through each set 1 interval
N1 = length(StartTimes1); %get total number of set 1 intervals
for i=1:N1
    %get start times of set 2 within set 1 interval
    StartTimes2_Within = StartTimes2(isbetween(StartTimes2,StartTimes1(i),EndTimes1(i)));
    StartTimesIntersecting = [StartTimesIntersecting; StartTimes2_Within]; %add to list of intersecting start times
    
    %get end times of set 2 corresponding to start times of set 2 within set 1 interval
    EndTimes2_StartTimes2_Within = EndTimes2(isbetween(StartTimes2,StartTimes1(i),EndTimes1(i)));

    %see which of these end times are within set 1 interval and which are outside
    %add to list of intersecting end times accordingly
    for j=1:length(StartTimes2_Within)
        if(EndTimes2_StartTimes2_Within(j)<=EndTimes1(i))
            EndTimesIntersecting = [EndTimesIntersecting; EndTimes2_StartTimes2_Within(j)];
        else
            EndTimesIntersecting = [EndTimesIntersecting; EndTimes1(i)];
        end
    end
end

%go through each set 2 interval
N2 = length(StartTimes2); %get total number of set 2 intervals
for i=1:N2 
    %get start times of set 1 within set 2 interval
    StartTimes1_Within = StartTimes1(isbetween(StartTimes1,StartTimes2(i),EndTimes2(i)));
    StartTimesIntersecting = [StartTimesIntersecting; StartTimes1_Within]; %add to list of intersecting start times
    
    %get end times of set 1 corresponding to start times of set 1 within set 2 interval
    EndTimes1_StartTimes1_Within = EndTimes1(isbetween(StartTimes1,StartTimes2(i),EndTimes2(i)));

    %see which of these end times are within set 2 interval and which are outside
    %add to list of intersecting end times accordingly
    for j=1:length(StartTimes1_Within)
        if(EndTimes1_StartTimes1_Within(j)<=EndTimes2(i))
            EndTimesIntersecting = [EndTimesIntersecting; EndTimes1_StartTimes1_Within(j)];
        else
            EndTimesIntersecting = [EndTimesIntersecting; EndTimes2(i)];
        end
    end
end

%get rid of repeats
StartTimesIntersecting = unique(StartTimesIntersecting);
EndTimesIntersecting = unique(EndTimesIntersecting);