%% SCRIPT TO LOOK AT THRESHOLD CROSSINGS

%% data loading -- can comment out if re-running script
%initialize
clear all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/Thresholds/'; %folder for plots

%information about sites for analysis
Sites = {'Jericoacoara','RanchoGuadalupe','Oceano'};
Markers = {'x','o','v'};
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

%parameters for threshold detection and plotting windows
T_min = duration(0,0,1); %minimum duration for an event (s)
T_window = duration(0,0,5); %duration of window for plotting (+/- s)
dt_Wenglor = duration(0,0,0.04); %Wenglor timestep
t_Wenglor_window = seconds(-T_window):seconds(dt_Wenglor):seconds(T_window);

%% initialize cell arrays of values for all sites
t_wind_windows_all = cell(N_Sites,1);
q_init_windows_all = cell(N_Sites,1);
q_cess_windows_all = cell(N_Sites,1);
u_init_windows_all = cell(N_Sites,1);
u_cess_windows_all = cell(N_Sites,1);
z_init_windows_all = cell(N_Sites,1);
z_cess_windows_all = cell(N_Sites,1);
d_init_windows_all = cell(N_Sites,1);
d_cess_windows_all = cell(N_Sites,1);
zd_init_windows_all = cell(N_Sites,1); %z/d ratio
zd_cess_windows_all = cell(N_Sites,1); %z/d ratio

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    %% initialize matrices of values for site
    q_init_windows_all{i} = [];
    q_cess_windows_all{i} = [];
    u_init_windows_all{i} = [];
    u_cess_windows_all{i} = [];
    z_init_windows_all{i} = [];
    z_cess_windows_all{i} = [];
    d_init_windows_all{i} = [];
    d_cess_windows_all{i} = [];
    
    %% perform analysis
    %choose anemometer type and grain-size location based on site of interest
    if strcmp(Sites{i},'Oceano')
        AnemometerType = 'Sonic';
        GrainSizeLocation = 'Wenglor';
        dt_wind = duration(0,0,0.02);
    elseif strcmp(Sites{i},'RanchoGuadalupe')
        AnemometerType = 'Ultrasonic';
        GrainSizeLocation = 'Wenglor';
        dt_wind = duration(0,0,0.04);
    elseif strcmp(Sites{i},'Jericoacoara')
        AnemometerType = 'Ultrasonic';
        GrainSizeLocation = 'BSNE_A';
        dt_wind = duration(0,0,0.04);
    end
    
    %extract wind data from overall processed data file
    WindData = Data{i}.ProcessedData.(AnemometerType);
    Anemometers = fieldnames(WindData);
    N_Anemometers = length(Anemometers);
    t_wind_windows_all{i} = seconds(-T_window):seconds(dt_wind):seconds(T_window);
    
    %extract Wenglor surface grain-size data from overall processed data file
    GrainSizeData = Data{i}.ProcessedData.GrainSize.Surface(strcmp({Data{i}.ProcessedData.GrainSize.Surface.Location},GrainSizeLocation));
    GrainSizeTimes = [GrainSizeData.CollectionTime];
    
    %extract flux data from overall processed data file
    FluxData = Data{i}.ProcessedData.TotalFlux;
    N_FluxIntervals = length(FluxData);
    
    %go through each flux interval
    for j=1:N_FluxIntervals
        
        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_FluxIntervals),', ',datestr(now)]
                
        %get surface grain size for this interval
        IntervalTime = mean([FluxData(1).StartTime, FluxData(1).EndTime]); %time of midpoint of flux interval
        GrainSizeTimeDifference = abs(GrainSizeTimes-IntervalTime); %difference between this time and time of grain size collections
        d50_Interval = GrainSizeData(GrainSizeTimeDifference==min(GrainSizeTimeDifference)).d_50_mm; %get surface grain size for time minimizing this difference
        
        %Wenglor info
        z_Wenglor = FluxData(j).z.calc; %get Wenglor heights for interval
        N_Wenglor = length(z_Wenglor); %get number of Wenglors
        t_Wenglor = FluxData(j).t.calc; %get times of Wenglor timeseries
        StartTime_Wenglor = min(t_Wenglor); %get first time for Wenglor timeseries
        EndTime_Wenglor = max(t_Wenglor); %get last time for Wenglor timeseries
        
        %get wind data for Wenglor interval
        [~, t_wind, IntervalN, IntervalInd] = ExtractVariableTimeInterval...
            (WindData.(Anemometers{1}),StartTime_Wenglor,EndTime_Wenglor,'u','int','int');
        u_lowest = WindData.(Anemometers{1})(IntervalN).u.int([IntervalInd{:}]);
        v_lowest = WindData.(Anemometers{1})(IntervalN).v.int([IntervalInd{:}]);
        w_lowest = WindData.(Anemometers{1})(IntervalN).w.int([IntervalInd{:}]);
        [u_lowest, v_lowest, w_lowest] = reorient_anemometers_vanboxel2004(u_lowest, v_lowest, w_lowest); %rotate instrument
        
        %get times and indices of Wenglor timeseries that overlap with wind timeseries
        [t_Wenglor,ind_Wenglor_wind,~] = intersect(t_Wenglor,t_wind);
                
        %go through each Wenglor
        for k=1:N_Wenglor
            
            %get flux data for Wenglor
            q_Wenglor = FluxData(j).qz.calc(ind_Wenglor_wind,k);
            
            %find indices of time windows associated with initiation and cessation
            [ind_init_windows, ind_cess_windows, t_init, t_cess] = fluxevent_detection(q_Wenglor,t_Wenglor,T_min,T_window);
            N_init_windows = length(t_init);
            N_cess_windows = length(t_cess);
            
            %Wenglor information about initiation and cessation events
            q_init_windows = q_Wenglor(ind_init_windows); %get fluxes for initiation windows
            q_cess_windows = q_Wenglor(ind_cess_windows); %get fluxes for cessation windows
            t_Wenglor_init_windows = t_Wenglor(ind_init_windows); %Wenglor times for initiation windows
            t_Wenglor_cess_windows = t_Wenglor(ind_cess_windows); %Wenglor times for cessation windows
            
            %rotate matrices if they contain only 1 window
            if N_init_windows == 1
                q_init_windows = q_init_windows';
                t_Wenglor_init_windows = t_Wenglor_init_windows';
            end
            if N_cess_windows == 1
                q_cess_windows = q_cess_windows';
                t_Wenglor_cess_windows = t_Wenglor_cess_windows';
            end
            
            %get wind velocities for initiation windows
            StartTime_init_windows = min(t_Wenglor_init_windows'); %start times for windows
            EndTime_init_windows = max(t_Wenglor_init_windows'); %end times for windows            
            u_init_windows = zeros(N_init_windows,2*round(T_window/dt_wind)+1);
            for l = 1:N_init_windows
                u_init_windows(l,:) = u_lowest((t_wind>=StartTime_init_windows(l))&(t_wind<=EndTime_init_windows(l)));
            end
        
            %get wind velocities for cessation windows
            StartTime_cess_windows = min(t_Wenglor_cess_windows'); %start times for windows
            EndTime_cess_windows = max(t_Wenglor_cess_windows'); %end times for windows
            u_cess_windows = zeros(N_cess_windows,2*round(T_window/dt_wind)+1);
            for l = 1:N_cess_windows
                u_cess_windows(l,:) = u_lowest((t_wind>=StartTime_cess_windows(l))&(t_wind<=EndTime_cess_windows(l)));
            end
            
            %% add to matrices of values for site if windows are not empty
            if N_init_windows~=0
                q_init_windows_all{i} = [q_init_windows_all{i}; q_init_windows];
                u_init_windows_all{i} = [u_init_windows_all{i}; u_init_windows];
                z_init_windows_all{i} = [z_init_windows_all{i}; ones(N_init_windows,1)*z_Wenglor(k)];
                d_init_windows_all{i} = [d_init_windows_all{i}; ones(N_init_windows,1)*d50_Interval];
            end
            if N_cess_windows~=0
                q_cess_windows_all{i} = [q_cess_windows_all{i}; q_cess_windows];
                u_cess_windows_all{i} = [u_cess_windows_all{i}; u_cess_windows];
                z_cess_windows_all{i} = [z_cess_windows_all{i}; ones(N_cess_windows,1)*z_Wenglor(k)];
                d_cess_windows_all{i} = [d_cess_windows_all{i}; ones(N_cess_windows,1)*d50_Interval];
            end
        end
    end
    %get ratios of z/d
    zd_init_windows_all{i} = 1000*z_init_windows_all{i}./d_init_windows_all{i}; %z/d ratio
    zd_cess_windows_all{i} = 1000*z_cess_windows_all{i}./d_cess_windows_all{i}; %z/d ratio
end

%% Create plots conditioned on z/d
zd_min = [0, 250];
zd_max = [250, 1000];

for i=1:length(zd_min)
    %initiation
    figure(i*2-1); clf;
    subplot(2,1,1); hold on;
    subplot(2,1,2); hold on;
    N_init_zd = cell(N_Sites,1);
    combined_legend = cell(N_Sites,1);
    for j=1:N_Sites
        zd_ind = find((zd_init_windows_all{j}>=zd_min(i))&(zd_init_windows_all{j}<=zd_max(i)));
        N_init_zd{j} = strcat('N = ',int2str(length(zd_ind)));
        combined_legend{j} = [Sites{j},', ',N_init_zd{j}];
        subplot(2,1,1);
        plot(t_wind_windows_all{j},mean(u_init_windows_all{j}(zd_ind,:)));
        subplot(2,1,2);
        plot(t_Wenglor_window,mean(q_init_windows_all{j}(zd_ind,:)));
    end
    subplot(2,1,1);
    xlabel('t_{init} (s)','FontSize',16);
    ylabel('u (m/s)','FontSize',16);
    set(gca,'FontSize',16);
    title(['z/d = ',int2str(zd_min(i)),' - ',int2str(zd_max(i)),'; T_{min} = ',int2str(seconds(T_min)*1000),' ms']);
    subplot(2,1,2);
    xlabel('t_{init} (s)','FontSize',16);
    ylabel('q (g/m^2/s)','FontSize',16);
    set(gca,'FontSize',12);
    h_legend = legend(combined_legend,'Location','Northwest');
    set(h_legend,'FontSize',12);
    print([folder_Plots,'init_zd_',int2str(zd_min(i)),'_',int2str(zd_max(i)),'_Tmin_',int2str(seconds(T_min)*1000),'.png'],'-dpng');

    %cessation
    figure(i*2); clf;
    subplot(2,1,1); hold on;
    subplot(2,1,2); hold on;
    N_cess_zd = cell(N_Sites,1);
    combined_legend = cell(N_Sites,1);
    for j=1:N_Sites
        zd_ind = find((zd_cess_windows_all{j}>=zd_min(i))&(zd_cess_windows_all{j}<=zd_max(i)));
        N_cess_zd{j} = strcat('N = ',int2str(length(zd_ind)));
        combined_legend{j} = [Sites{j},', ',N_cess_zd{j}];
        subplot(2,1,1);
        plot(t_wind_windows_all{j},mean(u_cess_windows_all{j}(zd_ind,:)));
        subplot(2,1,2);
        plot(t_Wenglor_window,mean(q_cess_windows_all{j}(zd_ind,:)));
    end
    subplot(2,1,1);
    xlabel('t_{cess} (s)','FontSize',16);
    ylabel('u (m/s)','FontSize',16);
    set(gca,'FontSize',16);
    title(['z/d = ',int2str(zd_min(i)),' - ',int2str(zd_max(i)),'; T_{min} = ',int2str(seconds(T_min)*1000),' ms']);
    subplot(2,1,2);
    xlabel('t_{cess} (s)','FontSize',16);
    ylabel('q (g/m^2/s)','FontSize',16);
    set(gca,'FontSize',16);
    h_legend = legend(combined_legend,'Location','Northeast');
    set(h_legend,'FontSize',12);
    print([folder_Plots,'cess_zd_',int2str(zd_min(i)),'_',int2str(zd_max(i)),'_Tmin_',int2str(seconds(T_min)*1000),'.png'],'-dpng');
end