%% SCRIPT TO LOOK AT REYNOLDS VERSUS LOG LAW SHEAR STRESS

%% data loading -- can comment out if re-running script
%initialize
clear all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/vonKarman/'; %folder for plots

%information about sites for analysis
Sites = {'Jericoacoara','RanchoGuadalupe','Oceano'};
Markers = {'bx','ro','gv'};
N_Sites = length(Sites);

for i = 1:N_Sites
    %load processed data
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Sites{i});
    Data{i} = load(ProcessedData_Path); %load processed data
    Metadata{i} = load(Metadata_Path); %load metadata
end

%% parameter values
%set physical parameters
rho_a = 1.23; %air density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%set fitting parameters
z_max_fit = 2.5; %maximum instrument height (m) for profile fit

%set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,15,0); %15 minutes

%% INITIALIZE CELL ARRAYS
%initialize cell arrays for all sites
z_profiles = cell(N_Sites,1); %initialize profiles of instrument heights
z_lowest_list = cell(N_Sites,1); %height for lowest anemometer
z_dudz_profiles = cell(N_Sites,1); %initialize profiles of heights for du/dz
AnemometerName_profiles = cell(N_Sites,1); %initialize profiles of instrument names

%initialize cell arrays of raw profiles
u_profiles_raw = cell(N_Sites,1); %initialize list of velocity profiles
ustRe_profiles_raw = cell(N_Sites,1); %initialize list of Reynolds shear velocities
ustRe_dudz_profiles_raw = cell(N_Sites,1); %initialize list of Reynolds shear velocities
heatflux_profiles_raw = cell(N_Sites,1); %initialize list of heat fluxes
dudz_profiles_raw = cell(N_Sites, 1); %inialize list of du/dz profiles
ustdudz_profiles_raw = cell(N_Sites,1); %initialize list of expected shear velocities from du/dz - raw
ustdudzStab_profiles_raw = cell(N_Sites,1); %initialize list of expected shear velocities from du/dz - raw, stability corrected

%initialize cell arrays of calibrated profiles
u_profiles_cal = cell(N_Sites,1); %initialize list of calibrated velocity profiles
ustRe_profiles_cal = cell(N_Sites,1); %initialize list of calibrated Reynolds shear velocities
ustRe_dudz_profiles_cal = cell(N_Sites,1); %initialize list of calibrated Reynolds shear velocities
heatflux_profiles_cal = cell(N_Sites,1); %initialize list of calibrated heat fluxes
dudz_profiles_cal = cell(N_Sites, 1); %inialize list of du/dz profiles
ustdudz_profiles_cal = cell(N_Sites,1); %initialize list of expected shear velocities from du/dz - calibrated
ustdudzStab_profiles_cal = cell(N_Sites,1); %initialize list of expected shear velocities from du/dz - calibrated, stability corrected

%initialize lists of raw profile aggregate values
ustLog_list_raw = cell(N_Sites,1); %raw u* for log law
z0_list_raw = cell(N_Sites,1); %raw z0 for log law
ustLogStab_list_raw = cell(N_Sites,1); %stability-corrected u* for log law
z0Stab_list_raw = cell(N_Sites,1); %stability-corrected z0 for log law
ustRe_lowest_raw = cell(N_Sites,1); %raw u* for Reynolds stress of lowest anemometer

%initialize lists of calibrated profile aggregate values
ustLog_list_cal = cell(N_Sites,1); %calibrated u* for log law
z0_list_cal = cell(N_Sites,1); %calibrated z0 for log law
ustLogStab_list_cal = cell(N_Sites,1); %stability-corrected calibrated u* for log law
z0Stab_list_cal = cell(N_Sites,1); %stability-corrected calibrated z0 for log law
ustRe_lowest_cal = cell(N_Sites,1); %calibrated u* for Reynolds stress of lowest anemometer

%initialize lists of flux values
Q_list = cell(N_Sites,1); %total flux
zbar_list = cell(N_Sites,1); %flux height

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites

    %choose anemometer type based on site of interest
    if strcmp(Sites{i},'Oceano')
        AnemometerType = 'Sonic';
    elseif strcmp(Sites{i},'RanchoGuadalupe')
        AnemometerType = 'Ultrasonic';
    elseif strcmp(Sites{i},'Jericoacoara')
        AnemometerType = 'Ultrasonic';
    end
    
    %extract wind data from overall processed data file
    WindData = Data{i}.ProcessedData.(AnemometerType);
    Anemometers = fieldnames(WindData);
    N_Anemometers = length(Anemometers);
    
    %get start times and end times for cup anemometers and sonics
    ind_Wind = find(strcmp(Metadata{i}.InstrumentMetadata.InstrumentType,AnemometerType));
    InstrumentStartTimes = Metadata{i}.InstrumentMetadata.StartTime(ind_Wind);
    InstrumentEndTimes = Metadata{i}.InstrumentMetadata.EndTime(ind_Wind);

    %combine time intervals based on start and end times for individual instruments
    [CombinedStartTimes, CombinedEndTimes] = CombineIntervals(InstrumentStartTimes, InstrumentEndTimes);
    
    %create time blocks based on combined time intervals
    [BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(CombinedStartTimes, CombinedEndTimes, ProfileTimeInterval);
    N_Blocks = length(BlockStartTimes);
    
    %initialize cell arrays of basic profiles
    z_profiles{i} = cell(N_Blocks,1); %initialize profiles of instrument heights
    z_lowest_list{i} = zeros(N_Blocks,1); %height for lowest anemometer
    z_dudz_profiles{i} = cell(N_Blocks,1); %initialize heights for du/dz profiles
    AnemometerName_profiles{i} = cell(N_Blocks,1); %initialize profiles of instrument names

    %initialize cell arrays of raw profiles
    u_profiles_raw{i} = cell(N_Blocks,1); %initialize list of velocity profiles
    ustRe_profiles_raw{i} = cell(N_Blocks,1); %initialize list of Reynolds shear velocities
    ustRe_dudz_profiles_raw{i} = cell(N_Blocks,1); %initialize list of Reynolds shear velocities
    heatflux_profiles_raw{i} = cell(N_Blocks,1); %initialize list of heat fluxes
    dudz_profiles_raw{i} = cell(N_Blocks, 1); %inialize list of du/dz profiles
    ustdudz_profiles_raw{i} = cell(N_Blocks,1); %initialize list of expected shear velocities from du/dz - raw
    ustdudzStab_profiles_raw{i} = cell(N_Blocks,1); %initialize list of expected shear velocities from du/dz - raw, stability corrected
    
    %initialize cell arrays of calibrated profiles
    u_profiles_cal{i} = cell(N_Blocks,1); %initialize list of calibrated velocity profiles
    ustRe_profiles_cal{i} = cell(N_Blocks,1); %initialize list of calibrated Reynolds shear velocities
    ustRe_dudz_profiles_cal{i} = cell(N_Blocks,1); %initialize list of calibrated Reynolds shear velocities
    heatflux_profiles_cal{i} = cell(N_Blocks,1); %initialize list of calibrated heat fluxes
    dudz_profiles_cal{i} = cell(N_Blocks, 1); %inialize list of du/dz profiles
    ustdudz_profiles_cal{i} = cell(N_Blocks,1); %initialize list of expected shear velocities from du/dz - calibrated
    ustdudzStab_profiles_cal{i} = cell(N_Blocks,1); %initialize list of expected shear velocities from du/dz - calibrated, stability corrected
    
    %initialize lists of raw profile aggregate values
    ustLog_list_raw{i} = zeros(N_Blocks,1); %raw u* for log law
    z0_list_raw{i} = zeros(N_Blocks,1); %raw z0 for log law
    ustLogStab_list_raw{i} = zeros(N_Blocks,1); %stability-corrected u* for log law
    z0Stab_list_raw{i} = zeros(N_Blocks,1); %stability-corrected z0 for log law
    ustRe_lowest_raw{i} = zeros(N_Blocks,1); %u* for Reynolds stress of lowest anemometer
    
    %initialize lists of calibrated profile aggregate values
    ustLog_list_cal{i} = zeros(N_Blocks,1); %calibrated u* for log law
    z0_list_cal{i} = zeros(N_Blocks,1); %calibrated z0 for log law
    ustLogStab_list_cal{i} = zeros(N_Blocks,1); %stability-corrected calibrated u* for log law
    z0Stab_list_cal{i} = zeros(N_Blocks,1); %stability-corrected calibrated z0 for log law
    ustRe_lowest_cal{i} = zeros(N_Blocks,1); %calibrated u* for Reynolds stress

    %initialize lists of flux values
    Q_list{i} = zeros(N_Blocks,1)*NaN; %total flux
    zbar_list{i} = zeros(N_Blocks,1)*NaN; %flux height
    
    %go through each time block
    for j = 1:N_Blocks
        
        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks)] 
        
        %get start and end times for block
        StartTime = BlockStartTimes(j);
        EndTime = BlockEndTimes(j);
        
        %initialize profiles - basic values
        z_profile = [];
        AnemometerName_profile = [];

        %initialize profiles - raw
        u_profile_raw = [];
        T_profile_raw = [];
        ustRe_profile_raw = [];
        heatflux_profile_raw = [];

        %initialize profiles - calibrated
        u_profile_cal = [];
        T_profile_cal = [];
        ustRe_profile_cal = [];
        heatflux_profile_cal = [];
        
        %go through each anemometer
        for k = 1:N_Anemometers
            
            %get times and indices for variable in time interval
            [~, t, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData.(Anemometers{k}),StartTime,EndTime,'u','int','int');
            
            %make computations if data cover full interval
            if (length(IntervalN)==1)&&((min(t)==StartTime&&max(t)==EndTime))
                
                %get velocity components
                u_raw = WindData.(Anemometers{k})(IntervalN).u.int(IntervalInd{:});
                v_raw = WindData.(Anemometers{k})(IntervalN).v.int(IntervalInd{:});
                w_raw = WindData.(Anemometers{k})(IntervalN).w.int(IntervalInd{:});

                %rotate instrument
                [u_raw, v_raw, w_raw] = reorient_anemometers_vanboxel2004(u_raw, v_raw, w_raw); %rotate instrument
                
                %get temperatures
                T_raw = WindData.(Anemometers{k})(IntervalN).T.int(IntervalInd{:});
               
                %get instrument height
                z = WindData.(Anemometers{k})(IntervalN).InstrumentHeight.z;
                
                %compute calibrated values
                CalibrationFactor_u = WindData.(Anemometers{k})(IntervalN).u.CalibrationFactor;
                CalibrationFactor_v = WindData.(Anemometers{k})(IntervalN).v.CalibrationFactor;
                CalibrationFactor_w = WindData.(Anemometers{k})(IntervalN).w.CalibrationFactor;
                CalibrationFactor_T = WindData.(Anemometers{k})(IntervalN).T.CalibrationFactor;
                CalibrationIntercept_u = WindData.(Anemometers{k})(IntervalN).u.CalibrationIntercept;
                CalibrationIntercept_v = WindData.(Anemometers{k})(IntervalN).v.CalibrationIntercept;
                CalibrationIntercept_w = WindData.(Anemometers{k})(IntervalN).w.CalibrationIntercept;
                CalibrationIntercept_T = WindData.(Anemometers{k})(IntervalN).T.CalibrationIntercept;
                u_cal = (u_raw-CalibrationIntercept_u)/CalibrationFactor_u;
                v_cal = (v_raw-CalibrationIntercept_v)/CalibrationFactor_v;
                w_cal = (w_raw-CalibrationIntercept_w)/CalibrationFactor_w;
                T_cal = (T_raw-CalibrationIntercept_T)/CalibrationFactor_T;
                
                %make computations - raw values
                u_bar_raw = mean(u_raw);
                w_bar_raw = mean(w_raw);
                T_bar_raw = mean(T_raw);
                ustRe_kernal = mean((u_raw-u_bar_raw).*(w_raw-w_bar_raw));
                if ustRe_kernal<=0
                    ustRe_raw = sqrt(-ustRe_kernal);
                else
                    ustRe_raw = NaN;
                end
                heatflux_raw = mean((w_raw-w_bar_raw).*(T_raw-T_bar_raw));                
                
                %make computations - calibrated values
                u_bar_cal = mean(u_cal);
                w_bar_cal = mean(w_cal);
                T_bar_cal = mean(T_cal);
                ustRe_kernal = mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
                if ustRe_kernal<=0
                    ustRe_cal = sqrt(-ustRe_kernal);
                else
                    ustRe_cal = NaN;
                end
                heatflux_cal = mean((T_cal-T_bar_cal).*(w_cal-w_bar_cal));
                                
                %add basic values to profiles
                z_profile = [z_profile; z]; %add z to profile
                AnemometerName_profile = [AnemometerName_profile; Anemometers{k}]; %add instrument name to profile
                
                %add raw values to profiles
                u_profile_raw = [u_profile_raw; u_bar_raw]; %add u_bar to profile
                T_profile_raw = [T_profile_raw; T_bar_raw]; %add T_bar to profile
                ustRe_profile_raw = [ustRe_profile_raw; ustRe_raw]; %add tauRe to profile
                heatflux_profile_raw = [heatflux_profile_raw; heatflux_raw]; %add heatflux to profile
                               
                %add calibrated values to profiles
                u_profile_cal = [u_profile_cal; u_bar_cal]; %add u_bar_cal to profile
                T_profile_cal = [T_profile_cal; T_bar_cal]; %add T_bar_cal to profile
                ustRe_profile_cal = [ustRe_profile_cal; ustRe_cal]; %add tauRe_cal to profile
                heatflux_profile_cal = [heatflux_profile_cal; heatflux_cal]; %add heatflux to profile
            end
        end
        
        %sort profiles by height
        [z_profile, ind_sort] = sort(z_profile);
        AnemometerName_profile = AnemometerName_profile(ind_sort,:);

        %sort raw values
        u_profile_raw = u_profile_raw(ind_sort);
        T_profile_raw = T_profile_raw(ind_sort);
        ustRe_profile_raw = ustRe_profile_raw(ind_sort);
        heatflux_profile_raw = heatflux_profile_raw(ind_sort);

        %sort calibrated values
        u_profile_cal = u_profile_cal(ind_sort);
        T_profile_cal = T_profile_cal(ind_sort);
        ustRe_profile_cal = ustRe_profile_cal(ind_sort);
        heatflux_profile_cal = heatflux_profile_cal(ind_sort);

        %compute raw and calibrated stability parameter profiles
        zL_profile_raw = (-(g./(T_profile_raw+273.15)).*heatflux_profile_raw(1))./(ustRe_profile_raw(1).^3./(kappa*z_profile));
        zL_profile_cal = (-(g./(T_profile_cal+273.15)).*heatflux_profile_cal(1))./(ustRe_profile_cal(1).^3./(kappa*z_profile));
         
        %compute shear velocities from law of wall
        ind_fit = find(z_profile<=z_max_fit); %use only anemometers below height given by z_max_fit for u* calc
        z_fit = z_profile(ind_fit);
        u_fit_raw = u_profile_raw(ind_fit);
        u_fit_cal = u_profile_cal(ind_fit);
        zL_fit_raw = zL_profile_raw(ind_fit);
        zL_fit_cal = zL_profile_cal(ind_fit);
        
        %fit log-law using raw and calibrated values
        [ustLog_raw, z0_raw, ustLogStab_raw, z0Stab_raw] = FitLogLaw(z_fit,u_fit_raw,zL_fit_raw,kappa);
        [ustLog_cal, z0_cal, ustLogStab_cal, z0Stab_cal] = FitLogLaw(z_fit,u_fit_cal,zL_fit_cal,kappa);
    
        %compute du/dz within profiles
        dz_profile = diff(z_profile);
        du_profile_raw = diff(u_profile_raw); %raw du (m/s)
        du_profile_cal = diff(u_profile_cal); %calibrated du (m/s)
        dudz_profile_raw = du_profile_raw./dz_profile; %raw du/dz (1/s)
        dudz_profile_cal = du_profile_cal./dz_profile; %calibrated du/dz (1/s)
        z_dudz_profile = sqrt(z_profile(1:(end-1)).*z_profile(2:end)); %height for du/dz (m)

        %compute raw and calibrated stability parameter profiles for z of
        %du/dz profile by geometric average, then compute phi for similarity
        zL_dudz_profile_raw = -sqrt(zL_profile_raw(1:end-1).*zL_profile_raw(2:end));
        zL_dudz_profile_cal = -sqrt(zL_profile_cal(1:end-1).*zL_profile_cal(2:end));
        phi_dudz_profile_raw = (1-15*zL_dudz_profile_raw).^(-1/4);
        phi_dudz_profile_cal = (1-15*zL_dudz_profile_cal).^(-1/4);
        
        %compute expected ust based on du/dz and similarity parameter
        ustdudz_profile_raw = (kappa*z_dudz_profile.*dudz_profile_raw); %no stability correction
        ustdudz_profile_cal = (kappa*z_dudz_profile.*dudz_profile_cal); %no stability correction
        ustdudzStab_profile_raw = (kappa*z_dudz_profile.*dudz_profile_raw./phi_dudz_profile_raw);
        ustdudzStab_profile_cal = (kappa*z_dudz_profile.*dudz_profile_cal./phi_dudz_profile_cal);
        
        %compute ustRe at individual heights of du/dz by averaging
        ustRe_dudz_profile_raw = mean([ustRe_profile_raw(1:end-1),ustRe_profile_raw(2:end)]')';
        ustRe_dudz_profile_cal = mean([ustRe_profile_cal(1:end-1),ustRe_profile_cal(2:end)]')';

        %add individual profiles to cell arrays of profiles - general profiles
        z_profiles{i}{j} = z_profile;
        z_lowest_list{i}(j) = z_profile(1);
        z_dudz_profiles{i}{j} = z_dudz_profile;
        AnemometerName_profiles{i}{j} = AnemometerName_profile;
        
        %add individual profiles to cell arrays of profiles - raw profiles
        u_profiles_raw{i}{j} = u_profile_raw;
        ustRe_profiles_raw{i}{j} = ustRe_profile_raw;
        ustRe_dudz_profiles_raw{i}{j} = ustRe_dudz_profile_raw;
        heatflux_profiles_raw{i}{j} = heatflux_profile_raw;
        dudz_profiles_raw{i}{j} = dudz_profile_raw;
        ustdudz_profiles_raw{i}{j} = ustdudz_profile_raw;
        ustdudzStab_profiles_raw{i}{j} = ustdudzStab_profile_raw;
        ustdudz_profiles_cal{i}{j} = ustdudz_profile_cal;
        ustdudzStab_profiles_cal{i}{j} = ustdudzStab_profile_cal;
        
        %add individual profiles to cell arrays of profiles - calibrated profiles
        u_profiles_cal{i}{j} = u_profile_cal;
        ustRe_profiles_cal{i}{j} = ustRe_profile_cal;
        ustRe_dudz_profiles_cal{i}{j} = ustRe_dudz_profile_cal;
        heatflux_profiles_cal{i}{j} = heatflux_profile_cal;
        dudz_profiles_cal{i}{j} = dudz_profile_cal;

        %add values to lists - raw values
        ustLog_list_raw{i}(j) = ustLog_raw; %raw u* for log law
        z0_list_raw{i}(j) = z0_raw; %raw z0 for log law
        ustLogStab_list_raw{i}(j) = ustLogStab_raw; %raw stability corrected u* for log law
        z0Stab_list_raw{i}(j) = z0Stab_raw; %raw stability corrected z0 for log law
        ustRe_lowest_raw{i}(j) = ustRe_profile_raw(1); %u* for Reynolds stress
        
        %add values to lists - calibrated values
        ustLog_list_cal{i}(j) = ustLog_cal; %calibrated u* for log law
        z0_list_cal{i}(j) = z0_cal; %calibrated z0 for log law
        ustLogStab_list_cal{i}(j) = ustLogStab_cal; %calibrated stability corrected u* for log law
        z0Stab_list_cal{i}(j) = z0Stab_cal; %calibrated stability corrected z0 for log law
        ustRe_lowest_cal{i}(j) = ustRe_profile_cal(1); %calibrated u* for Reynolds stress
        
        %FLUX CALCULATIONS FOR INTERVAL
        %extract flux data from overall processed data file
        Flux = Data{i}.ProcessedData.TotalFlux;
        N_Flux = length(Flux);
        
        %extract time interval
        [~, t, IntervalN, IntervalInd] = ExtractVariableTimeInterval(Flux,StartTime,EndTime,'Q','calc','calc');
        
        %make computations if data cover full interval
        if (length(IntervalN)==1)&&((min(t)==StartTime&&max(t)==EndTime))
            qz_profile = mean(Flux(IntervalN).qz.calc([IntervalInd{:}],:));
            zq_profile = Flux(IntervalN).z.calc;
            [zbar,~,Q] = qz_profilefit(zq_profile,qz_profile);
            
            %convert to 0 if NaN
            if isnan(Q)
                Q=0;
            end
            if isnan(zbar)
                zbar=0;
            end
            
            %add to list
            Q_list{i}(j) = Q; %g/m/s
            zbar_list{i}(j) = zbar; %m
        end
    end
end

%plot u*log versus u*Reynolds (lowest anemometer) using calibrated values without stability correction
figure(1); clf; hold on;
for i=1:N_Sites
    plot(ustRe_lowest_cal{i},ustLog_list_cal{i},Markers{i});
end
plot([0 0.6],[0 0.6],'k');
h_legend = legend(Sites,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlabel('u_{*,Re,lowest} (m/s)','FontSize',16);
ylabel('u_{*,log} (m/s)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsLog_Uncorrected_lowest.png'],'-dpng');

%plot u*log versus u*Reynolds (1 m) using calibrated values without stability correction
figure(2); clf; hold on;
z_calc = 1; %1 meter height of anemometer for plot
for i=1:N_Sites
    ustRe_list = zeros(length(z_profiles{i}),1);
    for j=1:length(z_profiles{i})
        ind_1m = find(abs(z_profiles{i}{j}-z_calc)==min(abs(z_profiles{i}{j}-z_calc)));
        ustRe_list(j) = mean(ustRe_profiles_cal{i}{j}(ind_1m));
    end
    plot(ustRe_list,ustLog_list_cal{i},Markers{i});
end
plot([0 0.6],[0 0.6],'k');
h_legend = legend(Sites,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlabel('u_{*,Re,1m} (m/s)','FontSize',16);
ylabel('u_{*,log} (m/s)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsLog_Uncorrected_1m.png'],'-dpng');

%plot u*log versus u*Reynolds (2 m) using calibrated values without stability correction
figure(3); clf; hold on;
z_calc = 2; %1 meter height of anemometer for plot
for i=1:N_Sites
    ustRe_list = zeros(length(z_profiles{i}),1);
    for j=1:length(z_profiles{i})
        ind_1m = find(abs(z_profiles{i}{j}-z_calc)==min(abs(z_profiles{i}{j}-z_calc)));
        ustRe_list(j) = mean(ustRe_profiles_cal{i}{j}(ind_1m));
    end
    plot(ustRe_list,ustLog_list_cal{i},Markers{i});
end
plot([0 0.6],[0 0.6],'k');
h_legend = legend(Sites,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlabel('u_{*,Re,2m} (m/s)','FontSize',16);
ylabel('u_{*,log} (m/s)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsLog_Uncorrected_2m.png'],'-dpng');

%plot u*log versus u*Reynolds (lowest anemometer) using calibrated values with stability correction
figure(4); clf; hold on;
for i=1:N_Sites
    plot(ustRe_lowest_cal{i},ustLogStab_list_cal{i},Markers{i});
end
plot([0 0.6],[0 0.6],'k');
h_legend = legend(Sites,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlabel('u_{*,Re,lowest} (m/s)','FontSize',16);
ylabel('u_{*,log} (m/s) - Stability Corrected','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ReynoldsLog_Corrected.png'],'-dpng');

%plot ratio of u*log (stability corrected) over u*Reynolds (lowest anemometer) versus sediment flux
figure(5); clf; hold on;
for i=1:N_Sites
    plot(Q_list{i},ustLogStab_list_cal{i}./ustRe_lowest_cal{i},Markers{i});
end
plot([0 60],[1 1]);
xlim([0 60]);
ylim([0.75 1.5]);
h_legend = legend(Sites,'Location','SouthEast');
set(h_legend,'FontSize',16);
xlabel('Q (g/m/s)','FontSize',16);
ylabel('u_{*,log}/u_{*,Re,lowest}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'StressRatioFlux.png'],'-dpng');

%plot ratio of u*log (stability corrected) over u*Reynolds (lowest anemometer) versus saltation height
figure(6); clf; hold on;
for i=1:N_Sites
    plot(zbar_list{i},ustLogStab_list_cal{i}./ustRe_lowest_cal{i},Markers{i});
end
plot([0 0.2],[1 1]);
xlim([0 0.2]);
ylim([0.75 1.5]);
h_legend = legend(Sites,'Location','SouthEast');
set(h_legend,'FontSize',16);
xlabel('z_{salt} (m)','FontSize',16);
ylabel('u_{*,log}/u_{*,Re,lowest}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'StressRatioZsalt.png'],'-dpng');

%plot ratio of u*log (not stability corrected) over u*Reynolds (lowest anemometer) versus product of sediment flux and saltation height
figure(7); clf; hold on;
for i=1:N_Sites
    plot(Q_list{i}.*zbar_list{i},ustLog_list_cal{i}./ustRe_lowest_cal{i},Markers{i});
end
plot([0 10],[1 1]);
ylim([0.75 1.5]);
h_legend = legend(Sites,'Location','SouthEast');
set(h_legend,'FontSize',16);
xlabel('Q * z_{salt} (g/s)','FontSize',16);
ylabel('u_{*,log}/u_{*,Re,lowest}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'StressRatioQZ_notstabilitycorrected.png'],'-dpng');

%plot ratio of u*log (stability corrected) over u*Reynolds (lowest anemometer) versus product of sediment flux and saltation height
figure(8); clf; hold on;
for i=1:N_Sites
    plot(Q_list{i}.*zbar_list{i},ustLogStab_list_cal{i}./ustRe_lowest_cal{i},Markers{i});
end
plot([0 10],[1 1]);
ylim([0.75 1.5]);
h_legend = legend(Sites,'Location','SouthEast');
set(h_legend,'FontSize',16);
xlabel('Q * z_{salt} (g/s)','FontSize',16);
ylabel('u_{*,log}/u_{*,Re,lowest}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'StressRatioQZ_stabilitycorrected.png'],'-dpng');

%plot ratio of expected u* from local velocity profile (stability corrected) versus u*Reynolds (average of closest anemometers)
figure(9); clf; hold on;
ustRe_binedges = 0:0.05:0.6;
ustRe_bincenters = 0.025:0.05:0.575;
N_bins = length(ustRe_bincenters);
for i = 1:N_Sites
    ustRatio_binavgs = zeros(1,N_bins);
    ustdudzStab_list = vertcat(ustdudzStab_profiles_cal{i}{:});
    ustRe_list = vertcat(ustRe_dudz_profiles_cal{i}{:});
    zU_list = vertcat(z_dudz_profiles{i}{:});
    ustRatio_sublist = ustdudzStab_list(zU_list<=z_max_fit)./ustRe_list(zU_list<=z_max_fit); %sublist of u*expected/u*Re ratio filtered by z<z_max
    ustRe_sublist = ustRe_list(zU_list<=z_max_fit); %sublist of u*Re filtered by z<z_max
    for j = 1:N_bins
        ustRatio_bin = ustRatio_sublist(ustRe_sublist>=ustRe_binedges(j)&ustRe_sublist<=ustRe_binedges(j+1));
        ustRatio_binavgs(j) = median(ustRatio_bin(~isnan(ustRatio_bin)));
    end
    plot(ustRe_bincenters, ustRatio_binavgs,['-',Markers{i}]);
end
plot([0 0.6],[1 1]);
ylim([0 2]);
h_legend = legend(Sites,'Location','SouthEast');
set(h_legend,'FontSize',16);
xlabel('u_{*,Re,localavg} (m/s)','FontSize',16);
ylabel('(\kappa z)/(u_{*,Re,localavg})(du/dz)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ustdudzReRatio_ustRe.png'],'-dpng');

%plot ratio of expected u* from local velocity profile (stability corrected) versus Q (average of closest anemometers), conditioned on z_U/z_salt
figure(10); clf; hold on;
zuNorm_binedges = 0:15:60;
zuNorm_bincenters = 7.5:15:52.5;
LineStyles = {'-','--',':','-.'};
N_zubins = length(zuNorm_bincenters);
flux_binedges = 0:10:60;
flux_bincenters = 5:10:55;
N_fluxbins = length(flux_bincenters);
legend_entries = cell(N_Sites*N_zubins,1);
for i = 1:N_Sites
    ustdudzStab_list = vertcat(ustdudzStab_profiles_cal{i}{:}); %get combined list of u*dudz
    ustRe_list = vertcat(ustRe_dudz_profiles_cal{i}{:}); %get combined list of u*Re
    ustRatio_list = ustdudzStab_list./ustRe_list; %list of u*dudz/u*Re ratio
    N_profiles = length(zbar_list{i});
    zuNorm_list = []; %initialize list of normalized heights
    flux_list = []; %initialize list of fluxes
    for j = 1:N_profiles
        zuNorm_list = [zuNorm_list; z_dudz_profiles{i}{j}./zbar_list{i}(j)]; %add to list of normalized anemometer heights
        flux_list = [flux_list; Q_list{i}(j)*ones(size(z_dudz_profiles{i}{j}))]; %add to list of fluxes
    end
    
    %go through each normalized anemometer height bin
    ustRatio_binavgs = zeros(N_zubins,N_fluxbins); %initialize bin averages of ustdudz/ustRe
    for j = 1:N_zubins
        flux_sublist = flux_list(zuNorm_list>=zuNorm_binedges(j)&zuNorm_list<=zuNorm_binedges(j+1));
        ustRatio_sublist = ustRatio_list(zuNorm_list>=zuNorm_binedges(j)&zuNorm_list<=zuNorm_binedges(j+1));
        
        %go through each flux bin
        for k = 1:N_fluxbins
            ustRatio_bin = ustRatio_sublist(flux_sublist>=flux_binedges(k)&flux_sublist<=flux_binedges(k+1));
            ustRatio_binavgs(j,k) = median(ustRatio_bin);
        end
        plot(flux_bincenters, ustRatio_binavgs(j,:),[LineStyles{j},Markers{i}]);
        legend_entries{(N_zubins+1)*(i-1)+j}=[Sites{i},', z_{U}/z_{salt} = ',int2str(zuNorm_binedges(j)),' - ',int2str(zuNorm_binedges(j+1))];
    end
    
    %also, plot the case where Q = 0 and zuNorm is undefined
    ustRatio_0 = ustRatio_list(flux_list==0);
    plot(0,median(ustRatio_0(~isnan(ustRatio_0))),Markers{i},'MarkerSize',12)
    median(ustRatio_0(~isnan(ustRatio_0)))
    legend_entries{(N_zubins+1)*i}=[Sites{i},', Q = 0 g/m/s'];
end
plot([0 0.6],[1 1]);
h_legend = legend(legend_entries,'Location','EastOutside');
set(h_legend,'FontSize',10);
xlabel('Q (g/m/s)','FontSize',16);
ylabel('(\kappa z)/(u_{*,Re,localavg})(du/dz)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ustdudzReRatio_Q.png'],'-dpng');

%clear excess files and save remaining for future analysis
clear('Data','Metadata','WindData','Flux');
SaveData_Path = strcat(folder_ProcessedData,'vonKarmanData');
save(SaveData_Path);