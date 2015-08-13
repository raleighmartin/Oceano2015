%% function to interpolate and flag error values in raw data
% 'RawData' contains raw observations grouped by instrument
% Output 'InterpolatedData' is built on 'RawData', but adds fields to account for interpolations 
% Dependencies: NONE
% Used by: Processing_Master
 
function [InterpolatedData] = InterpolateData(RawData, InstrumentVariables)

%go through instrument types, first perform basic processing common to all instruments
InstrumentTypes = fieldnames(RawData); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type 
    Instruments = fieldnames(RawData.(InstrumentType)); %get list of instruments
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        Instrument = Instruments{j} %get current instrument
        N_Intervals = length(RawData.(InstrumentType).(Instrument)); %how many time intervals for this instrument
        
        Variables = unique(InstrumentVariables.VarNameGeneral(strcmp(InstrumentVariables.Instrument,Instrument))); %get variables for instrument
        N_Variables = length(Variables);
        
        for k = 1:N_Intervals
            t_raw = RawData.(InstrumentType).(Instrument)(k).t.raw; %get raw times
            
            %deal with repeated times
            [t_unique, ind_unique] = unique(t_raw); %find unique times and indices
            t_err = t_raw(setdiff(1:length(t_raw),ind_unique)); %initialize list of error times with repeat times
            
            %deal with time gaps
            dt = mode(diff(t_unique)); %get expected timestep
            StartTime = min(t_unique); %initial time
            EndTime = max(t_unique); %final time
            t_int = (StartTime:dt:EndTime)'; %create list of interpolated times
            [~, ind_nongap] = intersect(t_int, t_unique); %record indices of nongap times
            
            %now, interpolate individual variables
            for l = 1:N_Variables
                values_raw = RawData.(InstrumentType).(Instrument)(k).(Variables{l}).raw; %get raw values
                values_unique = values_raw(ind_unique); %remove repeats
                values_int = zeros(size(t_int))*NaN; %initialize vector of interpolated values as list of NaNs
                values_int(ind_nongap) = values_unique; %use values from "unique" vector for nongap locations
                good_ind = intersect(find(values_int~=-999),find(isnan(values_int)==0)); %points without -999 or NaN are good
                error_ind = union(find(values_int==-999),find(isnan(values_int)==1)); %find error points based on -999 or NaN
                N_error = length(error_ind); %how many error points?
                t_err = [t_err; t_int(error_ind)]; %add to list of error times
                t_err = unique(t_err); %keep only unique times in list
                
                %go through each error point and interpolate
                values_int(1) = values_int(min(good_ind)); %set first point 
                values_int(end) = values_int(max(good_ind)); %set last point
                for m = 1:N_error
                    ind = error_ind(m);
                    if ind~=1&&ind~=length(values_int) %ignore first and last point (set above)
                        ind_prev = max([good_ind(good_ind<ind);1]); %get index of previous good point (or 1)
                        ind_next = min([good_ind(good_ind>ind);length(values_int)]); %get index of next good point (or last point)
                        value_prev = values_int(ind_prev); %get value of last good point
                        value_next = values_int(ind_next); %get value of next good point
                        values_int(ind) = ... %perform interpolation
                            value_prev*(ind_next-ind)/(ind_next-ind+1)+...
                            value_next/(ind_next-ind+1);
                    end
                end
                RawData.(InstrumentType).(Instrument)(k).(Variables{l}).int = values_int; %add interpolation field to variable values
            end
            RawData.(InstrumentType).(Instrument)(k).t = struct('raw',t_raw,'int',t_int,'err',t_err); %record interpolated times (with raw and error times) in structured array
        end 
    end
end
    
InterpolatedData = RawData; %now that everything is done, rename "RawData" as "ProcessedData"