%% function to determine the start, end, and delay times for processing each logger table
% Dependencies: NONE
% Used by: LoggerImport

function tableTimes = GetTimes(Site, Date, Logger, LoggerTimes)

%find sites, dates, and loggers that match inputs
SiteMatch = strcmp(LoggerTimes.Site,Site);
DateMatch = LoggerTimes.Date==Date;
LoggerMatch = strcmp(LoggerTimes.Logger,Logger);

%get indices for time intervals meeting all criteria above
ind_intervals = find(SiteMatch.*DateMatch.*LoggerMatch); 

%how many intervals?
N_Intervals = length(ind_intervals);

%get start times, end times, delays and number of intervals
StartTime = LoggerTimes.StartTime(ind_intervals);
EndTime = LoggerTimes.EndTime(ind_intervals);
Delay_Seconds = LoggerTimes.Delay_Seconds(ind_intervals);

tableTimes = struct('StartTime',StartTime,'EndTime',EndTime,'Delay_Seconds',Delay_Seconds,'N_Intervals',N_Intervals);