%% function to take info about BSNE trap weights ('WeightBSNE') and grain sizes ('GrainSize_BSNE')
% and aggregate this into mass and grain-size profiles
% Dependencies: IntervalsStartsWithin, qzCalc, qz_profilefit
% Used by: ProcessingData

function ProfilesBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE,GrainSize_Surface)

%% get information about flux calculation time intervals, based on A1
indices_A1 = find(strcmp(WeightBSNE.NameBSNE,'A1')); %indices of A1 BSNE records
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
    
    %get associated date for this profile
    Date = Date_A1(i);
    
    %duration of each individual BSNE
    DurationBSNEs = EndTimesContained-StartTimesContained;
    
    %reduce intervals based on error code
    indices_Intervals = intersect(indices_Intervals,indices_nonerror); %indices to BSNE non-error intervals to be included in profile
    N_ProfileBSNEs = length(indices_Intervals); %how many BSNEs are there for profile?
    
    %create lists for flux and grain-size profiles
    name = WeightBSNE.NameBSNE(indices_Intervals); %get names of BSNEs for interval
    z = z_BSNE_m(indices_Intervals); %get heights of BSNEs for interval (m)
    sigmaz = WeightBSNE.BottomHeightErr_cm(indices_Intervals)*100; %get z errors for interval (m)
    d_StartTime_profile = WeightBSNE.StartTime_GrainSize(indices_Intervals);
    d_EndTime_profile = WeightBSNE.EndTime_GrainSize(indices_Intervals);
    
    %get associated grain-size (d10, d50, d90) values for profile and info about combined samples
    d_10_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_50_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_90_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_IsCombined_profile = zeros(N_ProfileBSNEs,1)*NaN;

    for j = 1:N_ProfileBSNEs
        %find index in array associated with grain size sample for BSNE
        ind_GrainSize = intersect(...
            intersect(find([GrainSize_BSNE(:).StartTime]==d_StartTime_profile(j)),...
            find([GrainSize_BSNE(:).EndTime]==d_EndTime_profile(j))),...
            find(strcmp({GrainSize_BSNE(:).NameBSNE},name(j))));
        
        %get d10, d50, d90, and sample information for this (could add other values later if desired)
        if ~isempty(ind_GrainSize) %add only if a value is found, otherwise default is NaN
            d_10_profile(j) = GrainSize_BSNE(ind_GrainSize).d_10_mm;
            d_50_profile(j) = GrainSize_BSNE(ind_GrainSize).d_50_mm;
            d_90_profile(j) = GrainSize_BSNE(ind_GrainSize).d_90_mm;
            d_IsCombined_profile(j) = strcmp(GrainSize_BSNE(ind_GrainSize).Notes,'combined');
        end
    end
        
    %Compute fluxes (g/m^2/s) and associated errors, sigmaqz (g/m^2/s)
    qz = zeros(N_ProfileBSNEs,1);
    sigmaqz = zeros(N_ProfileBSNEs,1);
    for j = 1:N_ProfileBSNEs
        qz(j) = qzCalc(WeightBSNE.Weight_g(indices_Intervals(j)),...
            WeightBSNE.HeightBSNE_cm(indices_Intervals(j)),...
            WeightBSNE.WidthBSNE_cm(indices_Intervals(j)),...
            DurationBSNEs(j));
        sigmaqz(j) = qz(j).*...
            WeightBSNE.Weight_g(indices_Intervals(j))./...
            WeightBSNE.WeightErr_g(indices_Intervals(j)); %error is total flux times relative error in mass of sample
    end
    
    %sort out values by height
    [z, sort_ind] = sort(z);
    sigmaz = sigmaz(sort_ind)
    qz = qz(sort_ind);
    sigmaqz = sigmaqz(sort_ind)
    name = name(sort_ind);
    d_10_profile = d_10_profile(sort_ind);
    d_50_profile = d_50_profile(sort_ind);
    d_90_profile = d_90_profile(sort_ind);
    d_IsCombined_profile = d_IsCombined_profile(sort_ind);
    d_StartTime_profile = d_StartTime_profile(sort_ind);
    d_EndTime_profile = d_EndTime_profile(sort_ind);
    
    %height-integrated flux from exponential fit
    [zbar,q0,Q] = qz_profilefit(z, qz, sigmaz, sigmaqz);
 
    %put all surface and profile grain size values together into structured array
    %d_Surface = struct('d_10',
    d_Profile = struct('d_10',d_10_profile,'d_50',d_50_profile,'d_90',d_90_profile,'IsCombined',d_IsCombined_profile,'StartTime',d_StartTime_profile','EndTime',d_EndTime_profile','units',{'mm'});
    
    %create structured array with units
    units = struct('z',{'m'}, 'qz',{'g/m^2/s'}, 'Q', {'g/m/s'}, 'q0', {'g/m^2/s'});
    
    %combine all info into structured array
    ProfilesBSNE{i} = struct('Site',Site_A1(i),'Date',Date_A1(i),...
        'StartTime',StartTime,'EndTime',EndTime,...
        'name',{name'},'z',z,...
        'qz',qz,'zbar',zbar,'Q',Q,'q0',q0,...
        'd_Profile',d_Profile,...
        'units',units);

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