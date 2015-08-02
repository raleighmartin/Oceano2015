%% SCRIPT TO COMPARE BSNE VALUES TO STRESSES
% Dependencies: 'CombineIntervals', 'IntervalsWithin',
% 'MaxOverlappingInterval', 'reorient_anemometers_vanboxel2004'

%%clear existing data and load processed data and metadata
% clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/BSNE/'; %folder for plots
saveMetadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
saveProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
load(saveMetadata_Path); %load metadata
% load(saveProcessedData_Path); %load processed data

%set physical parameters
rho_a = 1.23; %kg/m^3
kappa = 0.39; %von Karman parameter
z_max_fit = 2.5; %maximum instrument height for profile fit 

%set types of instruments for calculations
InstrumentTypes = {'Sonic'};
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

%create time blocks based on BSNE times
BlockStartTimes = [ProcessedData.BSNE.StartTime]';
BlockEndTimes = [ProcessedData.BSNE.EndTime]';
N_Blocks = length(BlockStartTimes);

%initialize cell arrays of profiles
u_profiles = cell(N_Blocks,1); %initialize list of velocity profiles
u_profiles_cal = cell(N_Blocks,1); %initialize list of calibrated velocity profiles
tauRe_profiles = cell(N_Blocks,1); %initialize list of Reynolds stresses
tauRe_profiles_cal = cell(N_Blocks,1); %initialize list of calibrated Reynolds stresses
z_profiles = cell(N_Blocks,1); %initialize list of instrument heights
theta_profiles = cell(N_Blocks,1); %initialize list of velocity angle profiles
InstrumentName_profiles = cell(N_Blocks,1); %initialize lists of instrument names

%initialize list of shear velocities
ust_log_list = zeros(N_Blocks,1); %raw u* for log law
z0_log_list = zeros(N_Blocks,1); %raw z0 for log law
ust_log_cal_list = zeros(N_Blocks,1); %calibrated u* for log law
z0_log_cal_list = zeros(N_Blocks,1); %calibrated z0 for log law
ust_Re_list = zeros(N_Blocks,1); %u* for Reynolds stress
ust_Re_cal_list = zeros(N_Blocks,1); %calibrated u* for Reynolds stress
theta_list = zeros(N_Blocks,1); %wind angle list

%initialize BSNE values
Q_list = zeros(N_Blocks,1);
zbar_list = zeros(N_Blocks,1);
zbar_max_list = zeros(N_Blocks,1);
zbar_min_list = zeros(N_Blocks,1);

%go through each block
for i = 1:N_Blocks
    BlockStartTime = BlockStartTimes(i);
    BlockEndTime = BlockEndTimes(i);
    
    %initialize profiles
    u_profile = []; %initialize profile of velocities
    u_profile_cal = []; %initialize profile of calibrated velocities
    tauRe_profile = []; %initialize list of Reynolds stresses
    tauRe_profile_cal = []; %initialize list of Reynolds stresses
    z_profile = []; %initialize profile of heights
    theta_profile = []; %initialize profile of velocity angles
    InstrumentName_profile = []; %initialize profile of instrument names
    for j = 1:N_InstrumentTypes
        InstrumentType = InstrumentTypes{j}; %set specific instrument type
        Instruments = fieldnames(ProcessedData.(InstrumentType));
        N_Instruments = length(Instruments);
        for k = 1:N_Instruments
            Instrument = Instruments{k}; %set specific instrument
            InstrumentCalibrationFactor = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'u')); %get calibration factor for instrument
            ind_InstrumentIntervalStartsBefore = find([ProcessedData.(InstrumentType).(Instrument).StartTime]<=BlockStartTime); %indices of instrument time intervals starting before block times
            ind_InstrumentIntervalEndsAfter = find([ProcessedData.(InstrumentType).(Instrument).EndTime]>=BlockEndTime); %indices of instrument time intervals ending after block times
            ind_InstrumentIntervalsContainingBlock = intersect(ind_InstrumentIntervalStartsBefore,ind_InstrumentIntervalEndsAfter); %indices of instrument time intervals containing block times
            if length(ind_InstrumentIntervalsContainingBlock)==1
                InstrumentInterval = ProcessedData.(InstrumentType).(Instrument)(ind_InstrumentIntervalsContainingBlock); %extract data for instrument time interval
                ind_TimesInsideBlock = find((InstrumentInterval.t.int>=BlockStartTime)&(InstrumentInterval.t.int<=BlockEndTime)); %get indices of times within block
                if strcmp(InstrumentType,'Cup')
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    u_bar = mean(InstrumentInterval.u.int(ind_TimesInsideBlock));
                    u_cal = u/InstrumentCalibrationFactor;
                    u_bar_cal = mean(u_cal);
                    theta = NaN;
                    tauRe = NaN;
                    tauRe_cal = NaN;
                elseif strcmp(InstrumentType,'Sonic')
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    v = InstrumentInterval.v.int(ind_TimesInsideBlock);
                    w = InstrumentInterval.w.int(ind_TimesInsideBlock);
                    theta = atan(mean(v)/mean(u))*180;
                    [u, v, w] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
                    u_bar = mean(u);
                    w_bar = mean(w);
                    tauRe = -rho_a*mean((u-u_bar).*(w-w_bar));
                    u_cal = u/InstrumentCalibrationFactor;
                    w_cal = w/InstrumentCalibrationFactor;
                    u_bar_cal = mean(u_cal);
                    w_bar_cal = mean(w_cal);
                    tauRe_cal = -rho_a*mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
                end
                u_profile = [u_profile; u_bar]; %add u_bar to profile
                u_profile_cal = [u_profile_cal; u_bar_cal]; %add u_bar_cal to profile
                tauRe_profile = [tauRe_profile; tauRe]; %add tauRe to profile
                tauRe_profile_cal = [tauRe_profile_cal; tauRe_cal]; %add tauRe to profile
                theta_profile = [theta_profile; theta]; %add theta to profile
                z_profile = [z_profile; InstrumentInterval.StartHeight_m+0.25]; %add z to profile -- right now this is just a kluge because I haven't performed calibration with distance sensor
                InstrumentName_profile = [InstrumentName_profile; Instrument]; %add instrument name to profile
            end
        end
    end
    
    %sort profiles by height
    [z_profile, ind_sort] = sort(z_profile);
    u_profile = u_profile(ind_sort);
    u_profile_cal = u_profile_cal(ind_sort);
    tauRe_profile = tauRe_profile(ind_sort);
    tauRe_profile_cal = tauRe_profile_cal(ind_sort);
    theta_profile = theta_profile(ind_sort);
    InstrumentName_profile = InstrumentName_profile(ind_sort);
    
    %compute shear velocities - law of wall
    profile_fit_ind = find(z_profile<=z_max_fit); %use only anemometers below height given above
    P = polyfit(log(z_profile(profile_fit_ind)),u_profile(profile_fit_ind),1);
    ust_log = P(1)*kappa; %raw u* for log law
    z0_log = exp(-P(2)/P(1)); %z0 for log law
    P = polyfit(log(z_profile(profile_fit_ind)),u_profile_cal(profile_fit_ind),1);
    ust_log_cal = P(1)*kappa; %calibrated u* for log law
    z0_log_cal = exp(-P(2)/P(1)); %calibrated z0 for log law
    ust_Re = sqrt(mean(tauRe_profile(profile_fit_ind))/rho_a);
    ust_Re_cal = sqrt(mean(tauRe_profile_cal(profile_fit_ind))/rho_a);
      
    %add individual profiles to cell arrays of profiles
    u_profiles{i} = u_profile;
    u_profiles_cal{i} = u_profile_cal;
    tauRe_profiles{i} = tauRe_profile;
    tauRe_profiles_cal{i} = tauRe_profile_cal;
    theta_profiles{i} = theta_profile;
    z_profiles{i} = z_profile;
    InstrumentName_profiles{i} = InstrumentName_profile;
    
    %add shear velocities to lists
    ust_log_list(i) = ust_log; %raw u* for log law
    ust_log_cal_list(i) = ust_log_cal; %calibrated u* for log law
    z0_log_list(i) = z0_log; %raw z0 for log law
    z0_log_cal_list(i) = z0_log_cal; %calibrated z0 for log law
    ust_Re_list(i) = ust_Re; %u* for Reynolds stress
    ust_Re_cal_list(i) = ust_Re_cal; %calibrated u* for Reynolds stress
    
    %get mean wind angle, add to list
    theta_list(i) = mean(theta_profile(~isnan(theta_profile)));
    
    %get BSNE values
    BSNE_i = ProcessedData.BSNE(i);  
    z = BSNE_i.z;
    qz = BSNE_i.qz;
    name = BSNE_i.name;
    q0 = BSNE_i.q0;
    Q = BSNE_i.Q;
    zbar = BSNE_i.zbar.mid;
    zbar_min = BSNE_i.zbar.min;
    zbar_max = BSNE_i.zbar.max;
    Q_time = mean([BSNE_i.StartTime; BSNE_i.EndTime]);
    
    %get fit values for mean, min, and max saltation heights
    z_fit = linspace(0,max(z),50);
    qz_fit = q0*exp(-z_fit/zbar);
    q0_min = mean(qz./(exp(-z/zbar_min)));
    qz_fit_min = q0_min*exp(-z_fit/zbar_min);
    q0_max = mean(qz./(exp(-z/zbar_max)));
    qz_fit_max = q0_max*exp(-z_fit/zbar_max);
    
    %plot profile
    figure(1); clf;
    
    %plot actual values
    semilogx(qz,z,'kx','MarkerSize',5); hold on;
    text(qz,z,name);
   
    %plot fit profiles
    plot(qz_fit,z_fit,'k')
    plot(qz_fit_min,z_fit,'r');
    plot(qz_fit_max,z_fit,'b');
    
    %plot mean heights
    plot(q0*exp(-1),zbar,'k+','MarkerSize',10)
    plot(q0_min*exp(-1),zbar_min,'r+','MarkerSize',10);
    plot(q0_max*exp(-1),zbar_max,'b+','MarkerSize',10);    
    
    %label plot
    title(datestr(Q_time_list(i)),'FontSize',16);
    xlabel('flux (g/m/s^2)','FontSize',16);
    ylabel('height (m)','FontSize',16);
    print([folder_Plots,'/Profiles/FluxProfile_',datestr(Q_time,'yyyy-mm-dd_HHMM'),'.png'],'-dpng');
    
    %add values to lists
    Q_list(i) = Q;
    zbar_list(i) = zbar;
    zbar_max_list(i) = zbar_max;
    zbar_min_list(i) = zbar_min;
    Q_time_list(i,1) = Q_time;
end


%get indices only of intervals in +/- 45 degree angle
onshorewind_ind = find(abs(theta_list)<=45);

%plot BSNE fluxes versus shear velocity - calibrated Reynolds
figure(1); clf;
plot(ust_Re_list(onshorewind_ind),Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,Re}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_Re_list(onshorewind_ind),Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsUstBSNE_Uncalibrated.png'],'-dpng');

%plot BSNE fluxes versus shear stress - calibrated Reynolds
figure(1); clf;
plot(rho_a*ust_Re_list(onshorewind_ind).^2,Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('\tau_{Re}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_Re_list(onshorewind_ind).^2,Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsStressBSNE_Uncalibrated.png'],'-dpng');

%plot BSNE fluxes versus shear velocity - calibrated Reynolds
figure(2); clf;
plot(ust_Re_cal_list(onshorewind_ind),Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,Re}','FontSize',14);
ylabel('Q_{BSNE} (g/m/s)','FontSize',14);
corr_matrix = corrcoef(ust_Re_cal_list(onshorewind_ind),Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsUstBSNE_Calibrated.png'],'-dpng');

%plot BSNE fluxes versus shear stress - calibrated Reynolds
figure(2); clf;
plot(rho_a*ust_Re_cal_list(onshorewind_ind).^2,Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('\tau_{Re}','FontSize',14);
ylabel('Q_{BSNE} (g/m/s)','FontSize',14);
corr_matrix = corrcoef(ust_Re_cal_list(onshorewind_ind).^2,Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsStressBSNE_Calibrated.png'],'-dpng');

%plot BSNE fluxes versus shear velocity - calibrated log
figure(3); clf;
plot(ust_log_list(onshorewind_ind),Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,log}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_log_list(onshorewind_ind),Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'LogUstBSNE_Uncalibrated.png'],'-dpng');

%plot BSNE fluxes versus shear stress - calibrated log
figure(3); clf;
plot(rho_a*ust_log_list(onshorewind_ind).^2,Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('\tau_{log}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_log_list(onshorewind_ind).^2,Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'LogStressBSNE_Uncalibrated.png'],'-dpng');

%plot BSNE fluxes versus shear velocity - uncalibrated log
figure(4); clf;
plot(ust_log_cal_list(onshorewind_ind),Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,log,cal}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_log_cal_list(onshorewind_ind),Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'LogUstBSNE_Calibrated.png'],'-dpng');

%plot BSNE fluxes versus shear stress - uncalibrated log
figure(4); clf;
plot(rho_a*ust_log_cal_list(onshorewind_ind).^2,Q_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('\tau_{log,cal}');
ylabel('Q_{BSNE} (g/m/s)');
corr_matrix = corrcoef(ust_log_cal_list(onshorewind_ind).^2,Q_list(onshorewind_ind));
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
set(gca,'FontSize',16);
print([folder_Plots,'LogStressBSNE_Calibrated.png'],'-dpng');


%plot flux heights versus shear velocity - calibrated Reynolds
figure(5); clf;
errorbar(ust_Re_list(onshorewind_ind),zbar_list(onshorewind_ind),...
    zbar_min_list(onshorewind_ind),zbar_max_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,Re}');
ylabel('z_{bar} (cm)');
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsZbar_Uncalibrated.png'],'-dpng');

%plot flux heights versus shear velocity - calibrated Reynolds
figure(6); clf;
errorbar(ust_Re_cal_list(onshorewind_ind),zbar_list(onshorewind_ind),...
    zbar_min_list(onshorewind_ind),zbar_max_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,Re}');
ylabel('z_{bar} (cm)');
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsZbar_Calibrated.png'],'-dpng');

%plot flux heights versus shear velocity - calibrated log
figure(7); clf;
errorbar(ust_log_list(onshorewind_ind),zbar_list(onshorewind_ind),...
    zbar_min_list(onshorewind_ind),zbar_max_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,log}');
ylabel('z_{bar} (cm)');
set(gca,'FontSize',16);
print([folder_Plots,'LogZbar_Uncalibrated.png'],'-dpng');

%plot flux heights versus shear velocity - uncalibrated log
figure(8); clf;
errorbar(ust_log_cal_list(onshorewind_ind),zbar_list(onshorewind_ind),...
    zbar_min_list(onshorewind_ind),zbar_max_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('u_{*,log,cal}');
ylabel('z_{bar} (cm)');
set(gca,'FontSize',16);
print([folder_Plots,'LogZbar_Calibrated.png'],'-dpng');

%plot BSNE height versus BSNE time
figure(9); clf;
plot(Q_time_list(onshorewind_ind),zbar_list(onshorewind_ind),'o','MarkerSize',5);
xlabel('time');
ylabel('z_{bar} (cm)');
set(gca,'FontSize',16);
print([folder_Plots,'TimeZbar.png'],'-dpng');
% 
% 
% %plot uncalibrated values of shear velocities
% figure(2); clf;
% plot(ust_Re_cal_list(onshorewind_ind),ust_log_cal_list(onshorewind_ind),'o','MarkerSize',5);
% hold on; plot([0 0.55],[0 0.55]);
% xlim([0 0.55]); ylim([0 0.55]);
% xlabel('u_{*,Re,cal}');
% ylabel('u_{*,log,cal}');
% set(gca,'FontSize',16);
% print([folder_Plots,'LogReynoldsComparison_calibrated.png'],'-dpng');
% 
% 
% %hybrid plot of shear velocities
% figure(3); clf;
% plot(ust_Re_list(onshorewind_ind),ust_log_cal_list(onshorewind_ind),'o','MarkerSize',5);
% hold on; plot([0 0.55],[0 0.55]);
% xlim([0 0.55]); ylim([0 0.55]);
% xlabel('u_{*,Re}');
% ylabel('u_{*,log,cal}');
% set(gca,'FontSize',16);
% print([folder_Plots,'LogReynoldsComparison_hybrid.png'],'-dpng');
% 
% 
% %sample velocity profile - cups
% profile_ind = find(abs(theta_list)==min(abs(theta_list))); %choose profile that has least angle deviation
% figure(4); clf;
% semilogy(u_profiles{profile_ind},z_profiles{profile_ind},'x',...
%     u_profiles_cal{profile_ind},z_profiles{profile_ind},'o','MarkerSize',10);
% hold on;
% plot(ust_log_cal_list(profile_ind)/kappa*log(z_profiles{profile_ind}/z0_log_cal_list(profile_ind)),z_profiles{profile_ind});
% ylim([0.5 10]);
% h_legend = legend('uncalibrated','calibrated','fit','Location','NorthWest');
% xlabel('u (m/s)');
% ylabel('z (m)');
% set(h_legend,'FontSize',16);
% set(gca,'FontSize',16);
% print([folder_Plots,'VelocityProfile_SonicsCups.png'],'-dpng');
% 
% 
% %sample velocity profile - no cups
% profile_ind = find(abs(theta_list)==min(abs(theta_list))); %choose profile that has least angle deviation
% sonic_ind = find(InstrumentName_profiles{profile_ind}=='S');
% figure(5); clf;
% semilogy(u_profiles{profile_ind}(sonic_ind),z_profiles{profile_ind}(sonic_ind),'x',...
%     u_profiles_cal{profile_ind}(sonic_ind),z_profiles{profile_ind}(sonic_ind),'o','MarkerSize',10);
% hold on;
% plot(ust_log_cal_list(profile_ind)/kappa*log(z_profiles{profile_ind}/z0_log_cal_list(profile_ind)),z_profiles{profile_ind});
% ylim([0.5 10]);
% h_legend = legend('uncalibrated','calibrated','fit','Location','NorthWest');
% xlabel('u (m/s)');
% ylabel('z (m)');
% set(h_legend,'FontSize',16);
% set(gca,'FontSize',16);
% print([folder_Plots,'VelocityProfile_SonicsOnly.png'],'-dpng');
% 
% 
% %sample momentum flux profile
% profile_ind = find(abs(theta_list)==min(abs(theta_list))); %choose profile that has least angle deviation
% figure(6); clf;
% semilogy(tauRe_profiles{profile_ind},z_profiles{profile_ind},'x',...
%     tauRe_profiles_cal{profile_ind},z_profiles{profile_ind},'o','MarkerSize',10);
% hold on;
% ylim([0.5 10]);
% h_legend = legend('uncalibrated','calibrated','Location','NorthWest');
% xlabel('tau_{Re} (Pa)');
% ylabel('z (m)');
% set(h_legend,'FontSize',16);
% set(gca,'FontSize',16);
% print([folder_Plots,'MomentumFluxProfile.png'],'-dpng');