%% SCRIPT TO GENERATE PROFILES OF ANEMOMETER DATA (NOT INCLUDING CUPS) AND PERFORM ANALYSES ASSOCIATED WITH THESE PROFILES

% Dependencies: 'CombineIntervals','CreateTimeBlocks','IntervalsWithin',
% 'MaxOverlappingInterval,'reorient_anemometers_vanboxel2004'

%%clear existing data and load processed data and metadata
clear all;
folder_ProcessedData = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/AnemometerProfiles/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
load(Metadata_Path); %load meta data
load(ProcessedData_Path); %load processed data

%set physical parameters
rho_a = 1.23; %air density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2
P_psi = [-0.2473 -1.2570 -2.3943 -2.4641 0.0312]; %fourth-order polynomial coefficients to determine psi for stability correction from z/L, from Kaimal and Finnigan (1994) Table 1.1

%set fitting parameters
z_max_fit = 1.5; %maximum instrument height for profile fit
theta_max_onshore = 45; %maximum angle (deg) to be considered onshore wind

%set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,30,0); %30 minutes

%set types of instruments for profiles
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

%create time blocks based on combined time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(CombinedStartTimes, CombinedEndTimes, ProfileTimeInterval);
N_Blocks = length(BlockStartTimes);

%initialize cell arrays of profiles
u_profiles = cell(N_Blocks,1); %initialize list of velocity profiles
u_profiles_cal = cell(N_Blocks,1); %initialize list of calibrated velocity profiles
tauRe_profiles = cell(N_Blocks,1); %initialize list of Reynolds stresses
tauRe_profiles_cal = cell(N_Blocks,1); %initialize list of calibrated Reynolds stresses
heatflux_profiles = cell(N_Blocks,1); %initialize list of heat fluxes
heatflux_profiles_cal = cell(N_Blocks,1); %initialize list of calibrated heat fluxes
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

%go through each block
for i = 1;
%for i = 1:N_Blocks
    BlockStartTime = BlockStartTimes(i);
    BlockEndTime = BlockEndTimes(i);
    
    %initialize profiles
    u_profile = []; %initialize profile of velocities
    u_profile_cal = []; %initialize profile of calibrated velocities
    T_profile = []; %initialize profile of temperatures
    T_profile_cal = []; %initialize profile of calibrated temperatures
    tauRe_profile = []; %initialize list of Reynolds stresses
    tauRe_profile_cal = []; %initialize list of calibrated Reynolds stresses
    heatflux_profile = []; %initialize list of heat fluxes
    heatflux_profile_cal = []; %initialize list of calibrated heat fluxes
    z_profile = []; %initialize profile of heights
    theta_profile = []; %initialize profile of velocity angles
    InstrumentName_profile = []; %initialize profile of instrument names
    
    %generate profiles
    for j = 1:N_InstrumentTypes
        InstrumentType = InstrumentTypes{j}; %set specific instrument type
        Instruments = fieldnames(ProcessedData.(InstrumentType));
        N_Instruments = length(Instruments);
        
        %go through each instrument
        for k = 1:N_Instruments
            %set specific instrument, determine indices of instrument in processed data
            Instrument = Instruments{k};
            ind_InstrumentIntervalStartsBefore = find([ProcessedData.(InstrumentType).(Instrument).StartTime]<=BlockStartTime); %indices of instrument time intervals starting before block times
            ind_InstrumentIntervalEndsAfter = find([ProcessedData.(InstrumentType).(Instrument).EndTime]>=BlockEndTime); %indices of instrument time intervals ending after block times
            ind_InstrumentIntervalsContainingBlock = intersect(ind_InstrumentIntervalStartsBefore,ind_InstrumentIntervalEndsAfter); %indices of instrument time intervals containing block times
            
            %get calibration factors and intercepts
            CalibrationFactor_u = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'u'));
            CalibrationIntercept_u = InstrumentVariables.CalibrationIntercept(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'u'));
            CalibrationFactor_v = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'v'));
            CalibrationIntercept_v = InstrumentVariables.CalibrationIntercept(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'w'));
            CalibrationFactor_w = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'w'));
            CalibrationIntercept_w = InstrumentVariables.CalibrationIntercept(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'w'));
            CalibrationFactor_T = InstrumentVariables.CalibrationFactor(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'T'));
            CalibrationIntercept_T = InstrumentVariables.CalibrationIntercept(strcmp(InstrumentVariables.Instrument,Instrument)&strcmp(InstrumentVariables.VarNameGeneral,'T'));
            
            %check to make sure block actually contains the instrument of interest, if so, perform analysis
            if length(ind_InstrumentIntervalsContainingBlock)==1
                InstrumentInterval = ProcessedData.(InstrumentType).(Instrument)(ind_InstrumentIntervalsContainingBlock); %extract data for instrument time interval
                ind_TimesInsideBlock = find((InstrumentInterval.t.int>=BlockStartTime)&(InstrumentInterval.t.int<=BlockEndTime)); %get indices of times within block
                
                %get raw (interpolated) velocities
                u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                v = InstrumentInterval.v.int(ind_TimesInsideBlock);
                w = InstrumentInterval.w.int(ind_TimesInsideBlock);
                
                %get other values
                T = InstrumentInterval.T.int(ind_TimesInsideBlock); %temperature
                z = InstrumentInterval.InstrumentHeight.z; %instrument height
                
                theta = atan(mean(v)/mean(u))*180/pi; %calculate angle of unrotated instrument
                [u, v, w] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
                
                %make computations
                u_bar = mean(u);
                w_bar = mean(w);
                T_bar = mean(T);
                tauRe = -rho_a*mean((u-u_bar).*(w-w_bar));
                heatflux = mean((w-w_bar).*(T-T_bar));                
                
                %get calibrated computations
                u_cal = (u-CalibrationIntercept_u)/CalibrationFactor_u;
                v_cal = (v-CalibrationIntercept_v)/CalibrationFactor_v;
                w_cal = (w-CalibrationIntercept_w)/CalibrationFactor_w;
                T_cal = (T-CalibrationIntercept_T)/CalibrationFactor_T;
                u_bar_cal = mean(u_cal);
                w_bar_cal = mean(w_cal);
                T_bar_cal = mean(T_cal);
                tauRe_cal = -rho_a*mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
                heatflux_cal = mean((T_cal-T_bar_cal).*(w_cal-w_bar_cal));
                
                %add values to profiles
                u_profile = [u_profile; u_bar]; %add u_bar to profile
                u_profile_cal = [u_profile_cal; u_bar_cal]; %add u_bar_cal to profile
                T_profile = [T_profile; T_bar]; %add T_bar to profile
                T_profile_cal = [T_profile_cal; T_bar_cal]; %add T_bar_cal to profile
                tauRe_profile = [tauRe_profile; tauRe]; %add tauRe to profile
                tauRe_profile_cal = [tauRe_profile_cal; tauRe_cal]; %add tauRe_cal to profile
                heatflux_profile = [heatflux_profile; heatflux]; %add heatflux to profile
                heatflux_profile_cal = [heatflux_profile_cal; heatflux_cal]; %add heatflux to profile
                theta_profile = [theta_profile; theta]; %add theta to profile
                z_profile = [z_profile; z]; %add z to profile
                InstrumentName_profile = [InstrumentName_profile; Instrument]; %add instrument name to profile
            end
        end
    end
    
    %sort profiles by height
    [z_profile, ind_sort] = sort(z_profile);
    u_profile = u_profile(ind_sort);
    u_profile_cal = u_profile_cal(ind_sort);
    T_profile = T_profile(ind_sort);
    T_profile_cal = T_profile_cal(ind_sort);
    tauRe_profile = tauRe_profile(ind_sort);
    tauRe_profile_cal = tauRe_profile_cal(ind_sort);
    heatflux_profile = heatflux_profile(ind_sort);
    heatflux_profile_cal = heatflux_profile_cal(ind_sort);
    theta_profile = theta_profile(ind_sort);
    InstrumentName_profile = InstrumentName_profile(ind_sort,:);
    
    %compute shear velocities from law of wall
    profile_fit_ind = find(z_profile<=z_max_fit); %use only anemometers below height given by z_max_fit for u* calc
    P = polyfit(log(z_profile(profile_fit_ind)),u_profile(profile_fit_ind),1); %fit to get parameters for log law
    ust_log = P(1)*kappa; %raw u* for log law
    z0_log = exp(-P(2)/P(1)); %z0 for log law
    P = polyfit(log(z_profile(profile_fit_ind)),u_profile_cal(profile_fit_ind),1); %fit to get parameters for calibrated log law
    ust_log_cal = P(1)*kappa; %calibrated u* for log law
    z0_log_cal = exp(-P(2)/P(1)); %calibrated z0 for log law
    
    %compute u* from Reynolds stress
    ust_Re_profile = sqrt(mean(tauRe_profile(1))/rho_a); %compute Reynolds stress u* profile
    ust_Re_profile_cal = sqrt(mean(tauRe_profile(1))/rho_a); %compute calibrated Reynolds stress u* profile
    ust_Re = ust_Re_profile(1); %use lowest Reynolds stress for u* calc
    ust_Re_cal = ust_Re_profile_cal(1); %calibrated Reynolds stress for u* calc
    
    %compute mean wind angle from profile
    theta_bar = mean(theta_profile);
    
    %compute z/L stability parameter profile
    zL_profile = (-(g./(T_profile+273.15))*heatflux_profile(1))./(ust_Re^3./(kappa*z_profile));
    zL_profile_cal = (-(g./(T_profile_cal+273.15))*heatflux_profile_cal(1))./(ust_Re_cal^3./(kappa*z_profile));
    
    %compute psi profile, law of wall stability correction
    psi_profile = polyval(P_psi,zL_profile);
    psi_profile_cal = polyval(P_psi,zL_profile_cal);

    %compute stability-corrected predicted velocity profile
    u_profile_pred = (ust_log/kappa)*(log(z_profile/z0_log)-psi_profile+psi_profile(1));
    u_profile_pred_cal = (ust_log_cal/kappa)*(log(z_profile/z0_log_cal)-psi_profile_cal+psi_profile_cal(1));

    %add individual profiles to cell arrays of profiles
    u_profiles{i} = u_profile;
    u_profiles_cal{i} = u_profile_cal;
    tauRe_profiles{i} = tauRe_profile;
    tauRe_profiles_cal{i} = tauRe_profile_cal;
    heatflux_profiles{i} = heatflux_profile;
    heatflux_profiles_cal{i} = heatflux_profile_cal;
    theta_profiles{i} = theta_profile;
    z_profiles{i} = z_profile;
    InstrumentName_profiles{i} = InstrumentName_profile;
    
    %add values to lists
    ust_log_list(i) = ust_log; %raw u* for log law
    ust_log_cal_list(i) = ust_log_cal; %calibrated u* for log law
    z0_log_list(i) = z0_log; %raw z0 for log law
    z0_log_cal_list(i) = z0_log_cal; %calibrated z0 for log law
    ust_Re_list(i) = ust_Re; %u* for Reynolds stress
    ust_Re_cal_list(i) = ust_Re_cal; %calibrated u* for Reynolds stress
    theta_list(i) = mean(theta_profile); %mean wind angle
    
    %plot velocity profile
    figure(1); clf;
    set(gcf,'Visible','off');
    semilogy(u_profile,z_profile,'bx',u_profile_cal,z_profile,'ro','MarkerSize',10); hold on;
    plot(u_profile_pred,z_profile,'b',u_profile_pred_cal,z_profile,'r');
    title([InstrumentInterval.Site,' ',...
        datestr(BlockStartTime,'yyyy-mm-dd HH:MM'),'-',...
        datestr(BlockEndTime,'HH:MM'),', ',...
        'u_{*} = ',num2str(ust_Re,2),' m/s, ',...
        '\theta = ',num2str(theta_bar,2),'\circ'],'FontSize',16);
    ylim([0.5 10]);
    xlabel('u (m/s)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend('Raw','Calibrated','Pred (raw)','Pred (cal)','Location','Northwest');
    set(h_legend,'Fontsize',16);
    set(gca,'FontSize',16);
    print([folder_Plots,'Profiles/VelocityProfile_',datestr(BlockStartTime,'yyyymmdd_HHMM'),'.png'],'-dpng');

    %plot temperature profile
    figure(2); clf;
    set(gcf,'Visible','off');
    semilogy(T_profile,z_profile,'x',T_profile_cal,z_profile,'o','MarkerSize',10); hold on;
    title([InstrumentInterval.Site,' ',...
        datestr(BlockStartTime,'yyyy-mm-dd HH:MM'),'-',...
        datestr(BlockEndTime,'HH:MM'),', ',...
        'u_{*} = ',num2str(ust_Re,2),' m/s, ',...
        '\theta = ',num2str(theta_bar,2),'\circ'],'FontSize',16);
    ylim([0.5 10]);
    xlabel('T ({\circ}C)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend('Raw','Calibrated','Location','Northeast');
    set(h_legend,'Fontsize',16);
    set(gca,'FontSize',16);
    print([folder_Plots,'Profiles/TemperatureProfile_',datestr(BlockStartTime,'yyyymmdd_HHMM'),'.png'],'-dpng');
    
    %plot momentum flux profile
    figure(3); clf;
    set(gcf,'Visible','off');
    semilogy(tauRe_profile,z_profile,'x',tauRe_profile_cal,z_profile,'o','MarkerSize',10); hold on;
    title([InstrumentInterval.Site,' ',...
        datestr(BlockStartTime,'yyyy-mm-dd HH:MM'),'-',...
        datestr(BlockEndTime,'HH:MM'),', ',...
        'u_{*} = ',num2str(ust_Re,2),' m/s, ',...
        '\theta = ',num2str(theta_bar,2),'\circ'],'FontSize',16);
    ylim([0.5 10]);
    xlabel('\tau_{Re} (Pa)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend('uncalibrated','calibrated','Location','NorthWest');
    set(h_legend,'FontSize',16);
    set(gca,'FontSize',16);
    print([folder_Plots,'Profiles/MomentumFluxProfile_',datestr(BlockStartTime,'yyyymmdd_HHMM'),'.png'],'-dpng');
    
    %plot heat flux profile
    figure(4); clf;
    set(gcf,'Visible','off');
    semilogy(heatflux_profile,z_profile,'x',heatflux_profile_cal,z_profile,'o','MarkerSize',10); hold on;
    title([InstrumentInterval.Site,' ',...
        datestr(BlockStartTime,'yyyy-mm-dd HH:MM'),'-',...
        datestr(BlockEndTime,'HH:MM'),', ',...
        'u_{*} = ',num2str(ust_Re,2),' m/s, ',...
        '\theta = ',num2str(theta_bar,2),'\circ'],'FontSize',16);
    ylim([0.5 10]);
    xlabel('\langle{w`T`}\rangle (K m/s)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend('uncalibrated','calibrated','Location','NorthWest');
    set(h_legend,'FontSize',16);
    set(gca,'FontSize',16);
    print([folder_Plots,'Profiles/HeatFluxProfile_',datestr(BlockStartTime,'yyyymmdd_HHMM'),'.png'],'-dpng');
end

%get indices only of intervals in +/- theta_max_onshore angle (onshore winds)
onshorewind_ind = find(abs(theta_list)<=theta_max_onshore);

%plot uncalibrated values of shear velocities
figure(1); clf;
plot(ust_Re_list(onshorewind_ind),ust_log_list(onshorewind_ind),'o','MarkerSize',5);
hold on; plot([0 0.55],[0 0.55]);
xlim([0 0.55]); ylim([0 0.55]);
xlabel('u_{*,Re}');
ylabel('u_{*,log}');
set(gca,'FontSize',16);
print([folder_Plots,'LogReynoldsComparison_uncalibrated.png'],'-dpng');

%plot calibrated values of shear velocities
figure(2); clf;
plot(ust_Re_cal_list(onshorewind_ind),ust_log_cal_list(onshorewind_ind),'o','MarkerSize',5);
hold on; plot([0 0.55],[0 0.55]);
xlim([0 0.55]); ylim([0 0.55]);
xlabel('u_{*,Re,cal}');
ylabel('u_{*,log,cal}');
set(gca,'FontSize',16);
print([folder_Plots,'LogReynoldsComparison_calibrated.png'],'-dpng');

%hybrid plot of shear velocities
figure(3); clf;
plot(ust_Re_list(onshorewind_ind),ust_log_cal_list(onshorewind_ind),'o','MarkerSize',5);
hold on; plot([0 0.55],[0 0.55]);
xlim([0 0.55]); ylim([0 0.55]);
xlabel('u_{*,Re}');
ylabel('u_{*,log,cal}');
set(gca,'FontSize',16);
print([folder_Plots,'LogReynoldsComparison_hybrid.png'],'-dpng');