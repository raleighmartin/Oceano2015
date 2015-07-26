%% SCRIPT TO GENERATE VELOCITY PROFILES, INCLUDING DATA FROM CUP ANEMOMETERS
% Dependencies: 'CombineIntervals', 'CreateTimeBlocks', 'IntervalsWithin',
% 'MaxOverlappingInterval', 'reorient_anemometers_vanboxel2004'

%clear existing data and load processed data and metadata
clear all;
folder_ProcessedData = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/AnemometerProfiles/VelocityProfilesWithCups/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
load(Metadata_Path); %load meta data
load(ProcessedData_Path); %load processed data

%set physical parameters
rho_a = 1.23; %kg/m^3

%set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,30,0); %30 minutes

%set types of instruments for profiles
InstrumentTypes = {'Sonic','Cup'};
N_InstrumentTypes = length(InstrumentTypes);

%get indices of these instrument types in metadata table
ind_AllInstrumentTypes = []; %initialize list
for i = 1:N_InstrumentTypes
    ind_InstrumentType = find(strcmp(InstrumentMetadata.InstrumentType,InstrumentTypes{i})); %find this instrument type
    ind_AllInstrumentTypes = [ind_AllInstrumentTypes; ind_InstrumentType]; %add to list
end

%get start times and end times for cup anemometers and sonics
InstrumentStartTimes = InstrumentMetadata.StartTime(ind_AllInstrumentTypes);
InstrumentEndTimes = InstrumentMetadata.EndTime(ind_AllInstrumentTypes);

%combine time intervals based on start and end times for individual instruments
[CombinedStartTimes, CombinedEndTimes] = CombineIntervals(InstrumentStartTimes, InstrumentEndTimes);

%create time blocks based on combined time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(CombinedStartTimes, CombinedEndTimes, ProfileTimeInterval);
N_Blocks = length(BlockStartTimes);

%go through each block
for i = 1:N_Blocks
    BlockStartTime = BlockStartTimes(i);
    BlockEndTime = BlockEndTimes(i);
    
    %initialize profiles
    u_profile = []; %initialize profile of velocities
    u_profile_cal = []; %initialize profile of calibrated velocities
    z_profile = []; %initialize profile of heights
    theta_profile = []; %initialize profile of velocity angles
    tauRe_profile = []; %initialize profile of Reynolds stresses
    InstrumentName_profile = []; %initialize profile of instrument names
    
    %generate profiles
    for j = 1:N_InstrumentTypes
        InstrumentType = InstrumentTypes{j}; %set specific instrument type
        Instruments = fieldnames(ProcessedData.(InstrumentType));
        N_Instruments = length(Instruments);
        for k = 1:N_Instruments
            Instrument = Instruments{k}; %set specific instrument
            InstrumentCalibrationFactor = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'u')); %get calibration factor for instrument
            InstrumentCalibrationIntercept = InstrumentVariables.CalibrationIntercept(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'u')); %get calibration factor for instrument
            ind_InstrumentIntervalStartsBefore = find([ProcessedData.(InstrumentType).(Instrument).StartTime]<=BlockStartTime); %indices of instrument time intervals starting before block times
            ind_InstrumentIntervalEndsAfter = find([ProcessedData.(InstrumentType).(Instrument).EndTime]>=BlockEndTime); %indices of instrument time intervals ending after block times
            ind_InstrumentIntervalsContainingBlock = intersect(ind_InstrumentIntervalStartsBefore,ind_InstrumentIntervalEndsAfter); %indices of instrument time intervals containing block times
            if length(ind_InstrumentIntervalsContainingBlock)==1
                InstrumentInterval = ProcessedData.(InstrumentType).(Instrument)(ind_InstrumentIntervalsContainingBlock); %extract data for instrument time interval
                ind_TimesInsideBlock = find((InstrumentInterval.t.int>=BlockStartTime)&(InstrumentInterval.t.int<=BlockEndTime)); %get indices of times within block
                if strcmp(InstrumentType,'Cup')
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    u_bar = mean(InstrumentInterval.u.int(ind_TimesInsideBlock));
                    u_cal = (u-InstrumentCalibrationIntercept)./InstrumentCalibrationFactor; %for cups, include both intercept and calibration factor
                    u_bar_cal = mean(u_cal);
                    tauRe = NaN;
                    theta = NaN;
                elseif strcmp(InstrumentType,'Sonic')
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    v = InstrumentInterval.v.int(ind_TimesInsideBlock);
                    w = InstrumentInterval.w.int(ind_TimesInsideBlock);
                    theta = atan(mean(v)/mean(u))*180/pi;
                    [u, v, w] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
                    u_bar = mean(u);
                    w_bar = mean(w);
                    tauRe = -rho_a*mean((u-u_bar).*(w-w_bar));
                    u_cal = (u-InstrumentCalibrationIntercept)./InstrumentCalibrationFactor; %for cups, include both intercept and calibration factor
                    u_bar_cal = mean(u_cal);
                end
                u_profile = [u_profile; u_bar]; %add u_bar to profile
                u_profile_cal = [u_profile_cal; u_bar_cal]; %add u_bar_cal to profile
                theta_profile = [theta_profile; theta]; %add theta to profile
                tauRe_profile = [tauRe_profile; tauRe]; %add tauRe to profile
                z_profile = [z_profile; InstrumentInterval.StartHeight_m+0.25]; %add z to profile -- right now this is just a kluge because I haven't performed calibration with distance sensor
                InstrumentName_profile = [InstrumentName_profile; Instrument]; %add instrument name to profile
            end
        end
    end
    
    %sort profiles by height
    [z_profile, ind_sort] = sort(z_profile);
    u_profile = u_profile(ind_sort);
    u_profile_cal = u_profile_cal(ind_sort);
    theta_profile = theta_profile(ind_sort);
    InstrumentName_profile = InstrumentName_profile(ind_sort,:);
    
    %get mean wind angle
    theta_bar = mean(theta_profile(~isnan(theta_profile)));
    
    %get shear velocity at bottom
    ust_S1 = sqrt(tauRe_profile(strcmp(cellstr(InstrumentName_profile),'S1'))/rho_a);
    
    %plot profile
    figure(1); clf;
    set(gcf,'Visible','off');
    semilogy(u_profile,z_profile,'x',u_profile_cal,z_profile,'o','MarkerSize',10); hold on;
    text(u_profile_cal+0.05,z_profile,InstrumentName_profile,'FontSize',16);
    title([InstrumentInterval.Site,' ',...
        datestr(BlockStartTime,'yyyy-mm-dd HH:MM'),'-',...
        datestr(BlockEndTime,'HH:MM'),', ',...
        'u_{*} = ',num2str(ust_S1,2),' m/s, ',...
        '\theta = ',num2str(theta_bar,2),'\circ'],'FontSize',16);
    ylim([0.5 10]);
    xlabel('u (m/s)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend('Raw','Calibrated','Location','Northwest');
    set(h_legend,'Fontsize',16);
    set(gca,'FontSize',16);
    print([folder_Plots,'VelocityProfile_WithCups_',datestr(BlockStartTime,'yyyymmdd_HHMM'),'.png'],'-dpng');
end