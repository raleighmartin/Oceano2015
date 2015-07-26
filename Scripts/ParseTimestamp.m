%% function to convert string elements from 'timestamp_list' and convert them into datatimes
% 'timestamp_datetime'
% Dependencies: NONE
% Used by: ReadLogger

function timestamp_datetime = ParseTimestamp(timestamp_list)

timestamp_short = datetime(timestamp_list,'InputFormat','"yyyy-MM-dd HH:mm:ss"'); %get short format datetime

if max(cellfun(@length,timestamp_list))==21; %if they are all short format, it's easy
    timestamp_datetime = timestamp_short; 
else %if some are long format, need to parse things out
    timestamp_long = datetime(timestamp_list,'InputFormat','"yyyy-MM-dd HH:mm:ss.SS"');

    timestamp_short_ind = find(~isnat(timestamp_short)); %get indices of short format timestamps
    timestamp_long_ind = find(~isnat(timestamp_long)); %get indices of long format timestamps

    timestamp_datetime_unsorted = [timestamp_short(timestamp_short_ind); timestamp_long(timestamp_long_ind)]; %combine the datetimes with both formats
    [~, timestamp_sort_ind] = sort([timestamp_short_ind; timestamp_long_ind]); %get order of indices for sorting
    timestamp_datetime = timestamp_datetime_unsorted(timestamp_sort_ind); %sort out the list of datetimes
end