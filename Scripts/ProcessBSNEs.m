%% function to take info about BSNE trap weights ('WeightBSNE') and aggregate into mass profiles
% contained in list of structured arrays, 'ProfilesBSNE'
% Dependencies: CellStrFind, IntervalsStartsWithin, qzCalc
% Used by: Processing_Master

function ProfilesBSNE = ProcessBSNEs(WeightBSNE)

%% get information about flux calculation time intervals, based on A1
indices_A1 = CellStrFind(WeightBSNE.NameBSNE,'A1'); %indices of A1 BSNE records
N_Intervals = length(indices_A1); %count intervals based on A1, which is in all profiles
StartTimes_A1 = WeightBSNE.StartTime(indices_A1); %get start times associated with A1
EndTimes_A1 = WeightBSNE.EndTime(indices_A1); %get end times associated with A1
Site_A1 = WeightBSNE.Site(indices_A1); %get sites from A1
Date_A1 = WeightBSNE.Date(indices_A1); %get dates from A1

%% initialize cell array of info for each interval
ProfilesBSNE = cell(N_Intervals,1);

%% compute heights of BSNEs (m) for all intervals, based on geometric mean height of trap opening
z_BSNE_m = sqrt(WeightBSNE.BottomHeight_cm.*(WeightBSNE.BottomHeight_cm+WeightBSNE.HeightBSNE_cm))/100;

%% get nonerror indices
indices_nonerror = find(WeightBSNE.ErrorCode==0);

%% go through each flux calculation time interval
for i = 1:N_Intervals;
    
    %look for all time intervals that start within A1 interval
    [StartTimesContained, EndTimesContained, indices_Intervals] = IntervalsStartsWithin(StartTimes_A1(i),EndTimes_A1(i),WeightBSNE.StartTime,WeightBSNE.EndTime);
    
    %generate minimal time intervals (i.e., include only portions that are mutually overlapping for all traps)
    StartTime = max(StartTimesContained);
    EndTime = min(EndTimesContained);
    
    %duration of each individual BSNE
    DurationBSNEs = EndTimesContained-StartTimesContained;
    
    %reduce intervals based on error code
    indices_Intervals = intersect(indices_Intervals,indices_nonerror);
    N_Intervals = length(indices_Intervals); %how many are there?
    
    %create lists for flux profiles
    name = WeightBSNE.NameBSNE(indices_Intervals); %get names of BSNEs for interval
    z = z_BSNE_m(indices_Intervals); %get heights of BSNEs for interval
    
    %Compute fluxes (g/m^2/s)
    qz = zeros(N_Intervals,1);
    for k = 1:N_Intervals
        qz(k) = qzCalc(WeightBSNE.Weight_g(indices_Intervals(k)),...
            WeightBSNE.HeightBSNE_cm(indices_Intervals(k)),...
            WeightBSNE.WidthBSNE_cm(indices_Intervals(k)),...
            DurationBSNEs(k));
    end
    
    %sort out values by height
    [z, sort_ind] = sort(z);
    qz = qz(sort_ind);
    name = name(sort_ind);
    
    %height-integrated flux from exponential fit
    fit_params = polyfit(z(qz~=0),log(qz(qz~=0)),1); %fit only to nonzero values of qz
    zbar = -1/fit_params(1); %m
    q0 = exp(fit_params(2)); %g/m^2/s
    Q = q0*zbar; %g/m/s

    %create structured array with units
    units = struct('z',{'m'}, 'qz',{'g/m^2/s'}, 'Q', {'g/m/s'}, 'zbar', {'m'}, 'q0', {'g/m^2/s'});
    
    %combine all info into structured array
    ProfilesBSNE{i} = struct('Site',Site_A1(i),'Date',Date_A1(i),'StartTime',StartTime,'EndTime',EndTime,'z',z,'qz',qz,'Q',Q,'zbar',zbar,'q0',q0,'units',units);

%     %plot it
%     figure(i);
%     semilogx(qz,z,'x');
%     text(qz,z,name);
%     title(datestr(StartTime))
%     xlabel('flux (g/m/s^2)');
%     ylabel('height (m)');
end

%% convert cell array to structured array
ProfilesBSNE = [ProfilesBSNE{1:end}]';