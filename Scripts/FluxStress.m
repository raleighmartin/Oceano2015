%% SCRIPT TO GENERATE PROFILES OF WENGLOR FLUX DATA AND COMPARE TO SHEAR VELOCITIES

%initialize
clear all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots
SaveData_Path = strcat(folder_ProcessedData,'FluxStressCalcs_5minuteInterval');

%information about sites for analysis
Sites = {'Jericoacoara','RanchoGuadalupe','Oceano'};
Markers = {'bx','ro','gv'};
N_Sites = length(Sites);

% %data loading -- can comment out if re-running script
% %load processed and metadata for each site, add to structured arrays of all data and metadata
% for i = 1:N_Sites
%     ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
%     Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Sites{i});
%     Data{i} = load(ProcessedData_Path); %load processed data
%     Metadata{i} = load(Metadata_Path); %load metadata
% end

%load literature data
LiteratureData_Path = [folder_ProcessedData,'GreeleyNamikasData'];
load(LiteratureData_Path);

%% parameter values
%set physical parameters
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

% %set time interval for computing velocity profiles
% ProfileTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
% N_RunningPerProfile = floor(ProfileTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times
% 
% %initialize lists of flux values
% Q_list = cell(N_Sites,1); %total flux
% zbar_list = cell(N_Sites,1); %flux height
% 
% %initialize lists of wind values
% ustRe_raw_list = cell(N_Sites,1); %raw u* for Reynolds stress
% tauRe_raw_list = cell(N_Sites,1); %raw tau for Reynolds stress
% ustRe_cal_list = cell(N_Sites,1); %calibrated u* for Reynolds stress
% tauRe_cal_list = cell(N_Sites,1); %calibrated tau for Reynolds stress
% 
% %initialize lists of grain sizes and dates
% date_list = cell(N_Sites,1); %lists of dates corresponding to calculated values
% d50_list = cell(N_Sites,1); %lists of surface grain sizes corresponding to calculated values
% 
% % subset lists excluding Q = negative, inf, or NaN entries
% zbar_list_subset = cell(N_Sites,1);
% Q_list_subset = cell(N_Sites,1);
% ustRe_raw_list_subset = cell(N_Sites,1);
% tauRe_raw_list_subset = cell(N_Sites,1);
% ustRe_cal_list_subset = cell(N_Sites,1);
% tauRe_cal_list_subset = cell(N_Sites,1);
% date_list_subset = cell(N_Sites,1);
% d50_list_subset = cell(N_Sites,1);
% 
% %% PERFORM ANALYSIS FOR EACH SITE
% for i = 1:N_Sites
%     
%     %choose anemometer type based on site of interest
%     if strcmp(Sites{i},'Oceano')
%         AnemometerType = 'Sonic';
%         Anemometer = 'S1';
%     elseif strcmp(Sites{i},'RanchoGuadalupe')
%         AnemometerType = 'Ultrasonic';
%         Anemometer = 'U1';
%     elseif strcmp(Sites{i},'Jericoacoara')
%         AnemometerType = 'Ultrasonic';
%         Anemometer = 'U1';
%     end
%     
%     %extract wind and flux data from overall processed data file
%     WindData = Data{i}.ProcessedData.(AnemometerType).(Anemometer);
%     FluxData = Data{i}.ProcessedData.TotalFlux;
%     
%     %get start times and end times for flux and wind observations
%     WindStartTimes = [WindData.StartTime]';
%     WindEndTimes = [WindData.EndTime]';
%     FluxStartTimes = [FluxData.StartTime]';
%     FluxEndTimes = [FluxData.EndTime]';
%     
%     %get start and end times for intersecting flux and wind intervals
%     [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes,WindEndTimes,FluxStartTimes,FluxEndTimes);
%     
%     %create time blocks based on intersecting time intervals
%     [BlockStartTimes, BlockEndTimes] = ...
%             CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, ProfileTimeInterval);
%     
%     %Add in additional blocks for running average offsets
%     for j = 1:N_RunningPerProfile
%         [BlockStartTimes_j, BlockEndTimes_j] = ...
%             CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, ProfileTimeInterval);
%         BlockStartTimes = [BlockStartTimes; BlockStartTimes_j];
%         BlockEndTimes = [BlockEndTimes; BlockEndTimes_j];
%     end
%     BlockStartTimes = sort(BlockStartTimes);
%     BlockEndTimes = sort(BlockEndTimes);
%     N_Blocks = length(BlockStartTimes);
%     
%     %initialize lists of flux values
%     Q_list{i} = zeros(N_Blocks,1); %total flux
%     zbar_list{i} = zeros(N_Blocks,1); %flux height
% 
%     %initialize lists of wind values
%     ustRe_raw_list{i} = zeros(N_Blocks,1); %raw u* for Reynolds stress
%     tauRe_raw_list{i} = zeros(N_Blocks,1); %raw tau for Reynolds stress
%     ustRe_cal_list{i} = zeros(N_Blocks,1); %calibrated u* for Reynolds stress
%     tauRe_cal_list{i} = zeros(N_Blocks,1); %calibrated tau for Reynolds stress
%     
%     %create list of dates
%     date_list{i} = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
%     
%     %initialize lists of grain sizes
%     d50_list{i} = zeros(N_Blocks,1); %lists of surface grain sizes corresponding to calculated values
% 
%     %go through time blocks
%     for j = 1:N_Blocks
% 
%         %display processing status
%         processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks),', ',datestr(now)]
%         
%         %get specific start and end time
%         StartTime = BlockStartTimes(j);
%         EndTime = BlockEndTimes(j);
%         
%         %FLUX CALCULATIONS FOR INTERVAL
%         %extract time interval
%         [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData,StartTime,EndTime,'Q','calc','calc');
%         
%         %restrict to time intervals with data -- THIS SHOULDN'T BE
%         %NECESSARY -- NEED TO GO BACK AND LOOK AT PROCESSING SCRIPT.  But
%         %for now just put NaNs where data are missing
%         if ~isempty(IntervalN)
%             %restrict to interval with most data (happens if interval is on edge of two datasets)
%             IntervalN = IntervalN(cellfun('length',IntervalInd)==max(cellfun('length',IntervalInd)));
%             IntervalInd = IntervalInd{cellfun('length',IntervalInd)==max(cellfun('length',IntervalInd))};
%         
%             %compute mean qz profile
%             qz_profile = mean(FluxData(IntervalN).qz.calc(IntervalInd,:));
%             zq_profile = FluxData(IntervalN).z.calc;
%             [zbar,q0,Q] = qz_profilefit(zq_profile,qz_profile);
%             
%             %convert to 0 if NaN
%             if isnan(Q)
%                 Q=0;
%             end
%             if isnan(zbar)
%                 zbar=0;
%             end
%             
%             %add to list
%             Q_list{i}(j) = Q; %g/m/s
%             zbar_list{i}(j) = zbar; %m
%         else
%             Q_list{i}(j) = NaN;
%             zbar_list{i}(j) = NaN;
%         end
%         
%         %WIND CALCULATIONS FOR INTERVAL
%         %extract time interval
%         [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'u','int','int');
% 
%         %get velocity values
%         u = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'u','int','int');
%         v = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'v','int','int');
%         w = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'w','int','int');
%         
%         %rotate instrument, call these 'raw' values
%         [u_raw, ~, w_raw] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
%         
%         %get calibration factors
%         CalibrationFactor_u = WindData(IntervalN).u.CalibrationFactor;
%         CalibrationFactor_w = WindData(IntervalN).w.CalibrationFactor;
%         CalibrationIntercept_u = WindData(IntervalN).u.CalibrationIntercept;
%         CalibrationIntercept_w = WindData(IntervalN).w.CalibrationIntercept;
%         
%         %apply calibration, call these 'cal' values
%         u_cal = (u_raw-CalibrationIntercept_u)/CalibrationFactor_u;
%         w_cal = (w_raw-CalibrationIntercept_w)/CalibrationFactor_w;
% 
%         %make computations - raw values
%         u_bar_raw = mean(u_raw);
%         w_bar_raw = mean(w_raw);
%         ustRe_kernal = mean((u_raw-u_bar_raw).*(w_raw-w_bar_raw));
%         if ustRe_kernal<=0
%             ustRe_raw = sqrt(-ustRe_kernal);
%         else
%             ustRe_raw = NaN;
%         end
%         
%         %make computations - calibrated values
%         u_bar_cal = mean(u_cal);
%         w_bar_cal = mean(w_cal);
%         ustRe_kernal = mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
%         if ustRe_kernal<=0
%             ustRe_cal = sqrt(-ustRe_kernal);
%         else
%             ustRe_cal = NaN;
%         end
%         
%         %add to list
%         ustRe_raw_list{i}(j) = ustRe_raw;
%         tauRe_raw_list{i}(j) = rho_a*ustRe_raw.^2;
%         ustRe_cal_list{i}(j) = ustRe_cal;
%         tauRe_cal_list{i}(j) = rho_a*ustRe_cal.^2;
%         
%         %GRAIN SIZE INFO
%         ind_date = find([Data{i}.ProcessedData.GrainSize.Surface.Date]==WindData(IntervalN).Date);
%         d50_list{i}(j) = median([Data{i}.ProcessedData.GrainSize.Surface(ind_date).d_50_mm]);
%     end
%     
%     % get subset lists excluding Q = negative, inf, or NaN entries
%     ind_good = find((~isnan(Q_list{i})&(Q_list{i}>=0))&Q_list{i}~=Inf);
%     zbar_list_subset{i} = zbar_list{i}(ind_good);
%     Q_list_subset{i} = Q_list{i}(ind_good);
%     ustRe_raw_list_subset{i} = ustRe_raw_list{i}(ind_good);
%     ustRe_cal_list_subset{i} = ustRe_cal_list{i}(ind_good);
%     tauRe_raw_list_subset{i} = tauRe_raw_list{i}(ind_good);
%     tauRe_cal_list_subset{i} = tauRe_cal_list{i}(ind_good);
%     date_list_subset{i} = date_list{i}(ind_good);
%     d50_list_subset{i} = d50_list{i}(ind_good);
% end
% 
% %% SAVE DATA
% clear Data;
% clear WindData;
% clear FluxData;
% save(SaveData_Path);

%% LOAD DATA IF NECESSARY
load(SaveData_Path);

%% SEPARATE AND PLOT RESULTS BY DAY
threshold_analysis_dates = [...
    datetime(2014,11,13);...
    datetime(2015,3,23);...
    datetime(2015,5,15);...
    datetime(2015,5,16);...
    datetime(2015,5,18);...
    datetime(2015,5,19);... %data on this day do not work with Jasper's script
    datetime(2015,5,23);...
    datetime(2015,5,24);...
    datetime(2015,5,27);...
    datetime(2015,6,2);...
    datetime(2015,6,3)];  %data on this day do not work with Jasper's script
tauThr_day = [];
ustThr_day = [];
ustThr_day_Jasper = [];
d50_threshold_day = [];
date_threshold_day = [];

for i = 1:N_Sites
    Dates = unique(date_list_subset{i});
    N_Dates = length(Dates);
    
    %go through each day
    for j = 1:N_Dates
        
        %get values for this date
        ind_date = find(date_list_subset{i}==Dates(j));
        tauRe_date = tauRe_cal_list_subset{i}(ind_date);
        ustRe_date = ustRe_cal_list_subset{i}(ind_date);
        Q_date = Q_list_subset{i}(ind_date);
        d50_date = d50_list_subset{i}(ind_date);
        
        %if among analysis days, fit to determine threshold, using only times with positive flux
        if ~isempty(find(threshold_analysis_dates==Dates(j)))
            %perform fit to get Threshold
            ind_Qpositive = find(Q_date>0);
            P_linear = polyfit(tauRe_date(ind_Qpositive),Q_date(ind_Qpositive),1);
            tauThr = -P_linear(2)/P_linear(1);
            ustThr = sqrt(tauThr/rho_a);

            %add to lists for each day
            tauThr_day = [tauThr_day; tauThr];
            ustThr_day = [ustThr_day; ustThr];
            d50_threshold_day = [d50_threshold_day; d50_date(1)];
            date_threshold_day = [date_threshold_day; Dates(j)];

            %get values for plotting fit
            tauRe_fit = linspace(tauThr,max(tauRe_date),50);
            ustRe_fit = sqrt(tauRe_fit/rho_a);
            Q_fit = polyval(P_linear,tauRe_fit);
            
            %use Jasper's script to get threshold, add to list
            try
                [ustThr_Jasper, ~] = determine_uit(ustRe_date', Q_date', rho_a);
            catch
                ustThr_Jasper = NaN;
            end
            ustThr_day_Jasper = [ustThr_day_Jasper; ustThr_Jasper];
        end
            
        %plot flux versus shear stress
        figure(1); clf;
        plot(tauRe_date,Q_date, Markers{i}); hold on;
        if ~isempty(find(threshold_analysis_dates==Dates(j)))
            plot(tauRe_fit,Q_fit);
        end
        title(datestr(Dates(j)),'FontSize',16);
        xlabel('\tau_{Re} (Pa)','FontSize',16);
        ylabel('Q (g/m/s)','FontSize',16);
        ylim([0 35]);
        set(gca,'FontSize',16);
        print([folder_Plots,'FluxTau_Wenglor_',datestr(Dates(j)),'.png'],'-dpng');
        
        %plot flux versus shear velocity
        figure(2); clf;
        plot(ustRe_date,Q_date, Markers{i}); hold on;
        if ~isempty(find(threshold_analysis_dates==Dates(j)))
            plot(ustRe_fit,Q_fit);
        end
        title(datestr(Dates(j)),'FontSize',16);
        xlabel('u_{*,Re} (m/s)','FontSize',16);
        ylabel('Q (g/m/s)','FontSize',16);
        ylim([0 35]);
        set(gca,'FontSize',16);
        print([folder_Plots,'FluxUst_Wenglor_',datestr(Dates(j)),'.png'],'-dpng');
        
        %for days with sufficient observations, plot threshold-normalized
        %fluxes
        if ~isempty(find(threshold_analysis_dates==Dates(j)))
            
            %get normalized quantities
            ind_highflux = find(Q_date>5); %flux must be at least 5 g/m/s for non-dimensionalization
            ust_nondim_date = ustRe_date(ind_highflux)/ustThr;
            tau_ex_date = tauRe_date(ind_highflux)-tauThr;
            Q_nondim2_date = Q_date(ind_highflux)*g./(1000*rho_a*tau_ex_date.*ustThr);
            Q_nondim3_date = Q_date(ind_highflux)*g./(1000*rho_a*tau_ex_date.*ustRe_date(ind_highflux));
            
            %plot square normalized flux versus excess shear velocity
            figure(3); clf;
            plot(ust_nondim_date,Q_nondim2_date, Markers{i});
            xlabel('u_{*}/u_{*,th}','FontSize',16);
            ylabel('gQ/{(\rho_{a}\tau_{ex}u_{*th})}','FontSize',16);
            title([datestr(Dates(j)),': Squared law']);
            set(gca,'FontSize',16);
            set(gca,'xscale','log');
            print([folder_Plots,'FluxNondim2_UstNonDim_',datestr(Dates(j)),'.png'],'-dpng');

            %plot cube normalized flux versus excess shear velocity
            figure(4); clf;
            plot(ust_nondim_date,Q_nondim3_date, Markers{i});
            xlabel('u_{*}/u_{*,th}','FontSize',16);
            ylabel('gQ/{(\rho_{a}\tau_{ex}u_{*})}','FontSize',16);
            title([datestr(Dates(j)),': Cubed law']);
            set(gca,'FontSize',16);
            set(gca,'xscale','log');
            print([folder_Plots,'FluxNondim3_UstNonDim_',datestr(Dates(j)),'.png'],'-dpng');
        end
    end
end

%get fit of u*thr versus shear velocity
P_linear = polyfit(d50_threshold_day,ustThr_day,1);
P_linear_Jasper = polyfit(d50_threshold_day,ustThr_day_Jasper,1);
P_sqrt = polyfit(sqrt(d50_threshold_day),ustThr_day,1);
P_sqrt_Jasper = polyfit(sqrt(d50_threshold_day),ustThr_day_Jasper,1);
C_sqrt = mean(ustThr_day./sqrt(d50_threshold_day));

%plot u*thr versus grain diameter
figure(3); clf;
for i = 1:N_Sites
    Dates = unique(date_list_subset{i});
    [~, ind_Site, ~] = intersect(date_threshold_day, Dates);

    d50_plot = d50_threshold_day(ind_Site);
    ustThr_plot = ustThr_day(ind_Site);
    ustThr_Jasper_plot = ustThr_day_Jasper(ind_Site);
    
    subplot(1,2,1); hold on;
    plot(d50_plot,ustThr_plot,Markers{i},'MarkerSize',10);
    title('original','FontSize',16);
    subplot(1,2,2); hold on;
    plot(d50_plot,ustThr_Jasper_plot,Markers{i},'MarkerSize',10);
    title('Jasper','FontSize',16);
end

%set up predictions
d50_range = linspace(min(d50_threshold_day),max(d50_threshold_day),50);

%set up legend
legend_values = Sites;
legend_values{length(Sites)+1} = 'linear fit - original';
legend_values{length(Sites)+2} = 'linear fit - Jasper';
legend_values{length(Sites)+3} = 'sqrt fit - original';
legend_values{length(Sites)+4} = 'sqrt fit - Jasper';
legend_values{length(Sites)+5} = 'sqrt fit - 0 intercept';

%add predictions and legend to both subplots
for i = 1:2
    subplot(1,2,i);
    plot(d50_range,polyval(P_linear,d50_range),'k');
    plot(d50_range,polyval(P_linear_Jasper,d50_range),'k-.');
    plot(d50_range,polyval(P_sqrt,sqrt(d50_range)),'c');
    plot(d50_range,polyval(P_sqrt_Jasper,sqrt(d50_range)),'c-.');
    plot(d50_range,C_sqrt*sqrt(d50_range),'b-','LineWidth',3);
    h_legend = legend(legend_values,'Location','NorthWest');
    set(h_legend,'FontSize',16);
    xlabel('d_{50} (mm)');
    ylabel('u_{*,thr} (m/s)');
    xlim([0.3 0.6]);
    ylim([0.23 0.4]);
    set(gca,'FontSize',16);
end
print([folder_Plots,'UstThr_D50_days.png'],'-dpng');


%% SEPARATE AND PLOT RESULTS BY GRAIN SIZE BINS

%set bin edges
d50_bins_min = [0.33, 0.39, 0.41, 0.51];
d50_bins_max = [0.36, 0.41, 0.44, 0.54];
N_bins = length(d50_bins_min);

%initialize list of calculated values for bins
tauThr_bins = zeros(N_bins,1);
ustThr_bins = zeros(N_bins,1);
ustThr_bins_Jasper = zeros(N_bins,1);
d50_threshold_bins = zeros(N_bins,1);

%go through each bin
for i = 1:N_bins

    %initialize lists for bins
    tauRe_bin = [];
    ustRe_bin = [];
    Q_bin = [];
    d50_bin = [];
    
    %go through each site
    for j = 1:N_Sites
        
        %get values for this bin for this site
        ind_bin_Site = find((d50_list_subset{j}>=d50_bins_min(i))&(d50_list_subset{j}<=d50_bins_max(i)));

        %if bin contains values for this Site, get values and add to overall list
        if ~isempty(ind_bin_Site)
            tauRe_bin_Site = tauRe_cal_list_subset{j}(ind_bin_Site);
            ustRe_bin_Site = ustRe_cal_list_subset{j}(ind_bin_Site);
            Q_bin_Site = Q_list_subset{j}(ind_bin_Site);
            d50_bin_Site = d50_list_subset{j}(ind_bin_Site);
            
            %add to lists for bins
            tauRe_bin = [tauRe_bin; tauRe_bin_Site];
            ustRe_bin = [ustRe_bin; ustRe_bin_Site];
            Q_bin = [Q_bin; Q_bin_Site];
            d50_bin = [d50_bin; d50_bin_Site];
        end
    end
    
    %for each bin, fit to determine threshold, using only times with positive flux
    ind_Qpositive = find(Q_bin>0);
    P_linear = polyfit(tauRe_bin(ind_Qpositive),Q_bin(ind_Qpositive),1);
    tauThr = -P_linear(2)/P_linear(1);
    ustThr = sqrt(tauThr/rho_a);
    
    %add to lists for each bin
    tauThr_bins(i) = tauThr;
    ustThr_bins(i) = ustThr;
    d50_threshold_bins(i) = mean(d50_bin);
    d50_print = [int2str(round(1000*mean(d50_bin))),'um'];
    
    %get values for plotting fit
    tauRe_fit = linspace(tauThr,max(tauRe_bin),50);
    ustRe_fit = sqrt(tauRe_fit/rho_a);
    Q_fit = polyval(P_linear,tauRe_fit);
    
    %use Jasper's script to get threshold, add to list
    [ustThr_Jasper, ~] = determine_uit(ustRe_bin', Q_bin', rho_a);
    ustThr_bins_Jasper(i) = ustThr_Jasper;
    
    %plot flux versus shear stress
    figure(4); clf;
    plot(tauRe_bin,Q_bin, 'o'); hold on;
    plot(tauRe_fit,Q_fit);
    title(['d = ',d50_print],'FontSize',16);
    xlabel('\tau_{Re} (Pa)','FontSize',16);
    ylabel('Q (g/m/s)','FontSize',16);
    ylim([0 35]);
    set(gca,'FontSize',16);
    print([folder_Plots,'FluxTau_Wenglor_d',d50_print,'.png'],'-dpng');

    %plot flux versus shear velocity
    figure(5); clf;
    plot(ustRe_bin,Q_bin, 'o'); hold on;
    plot(ustRe_fit,Q_fit);
    title(['d = ',d50_print],'FontSize',16);
    xlabel('u_{*,Re} (m/s)','FontSize',16);
    ylabel('Q (g/m/s)','FontSize',16);
    ylim([0 35]);
    set(gca,'FontSize',16);
    print([folder_Plots,'FluxUst_Wenglor_d',d50_print,'.png'],'-dpng');

end

%get fit of u*thr versus shear velocity
Pfit_ustThr_d50_linear = polyfit(d50_threshold_bins,ustThr_bins,1);
Pfit_ustThr_d50_linear_Jasper = polyfit(d50_threshold_bins, ustThr_bins_Jasper,1);
Pfit_ustThr_d50_sqrt = polyfit(sqrt(d50_threshold_bins),ustThr_bins,1);
Pfit_ustThr_d50_sqrt_Jasper = polyfit(sqrt(d50_threshold_bins), ustThr_bins_Jasper,1);
d50_range = linspace(min(d50_threshold_bins),max(d50_threshold_bins),50);
Cfit_ustThr_d50_sqrt = mean(ustThr_bins./sqrt(d50_threshold_bins));

%plot u*thr versus shear velocity
figure(6); clf; hold on;
subplot(1,2,1); hold on;
plot(d50_threshold_bins,ustThr_bins,'o','MarkerSize',15);
plot(d50_range,polyval(Pfit_ustThr_d50_linear,d50_range),'k');
plot(d50_range,polyval(Pfit_ustThr_d50_linear_Jasper,d50_range),'k-.');
plot(d50_range,polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_range)),'c');
plot(d50_range,polyval(Pfit_ustThr_d50_sqrt_Jasper,sqrt(d50_range)),'c-.');
plot(d50_range,Cfit_ustThr_d50_sqrt*sqrt(d50_range),'b-','LineWidth',3);

title('original','FontSize',16);
subplot(1,2,2); hold on;
plot(d50_threshold_bins,ustThr_bins_Jasper,'o','MarkerSize',15);
plot(d50_range,polyval(Pfit_ustThr_d50_linear,d50_range),'k');
plot(d50_range,polyval(Pfit_ustThr_d50_linear_Jasper,d50_range),'k-.');
plot(d50_range,polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_range)),'c');
plot(d50_range,polyval(Pfit_ustThr_d50_sqrt_Jasper,sqrt(d50_range)),'c-.');
plot(d50_range,Cfit_ustThr_d50_sqrt*sqrt(d50_range),'b-','LineWidth',3);

title('Jasper','FontSize',16);
for i = 1:2
    subplot(1,2,i); hold on;
    xlabel('d_{50} (mm)');
    ylabel('u_{*,thr} (m/s)');
    set(gca,'FontSize',16);
    h_legend = legend({'data','original fit','Jasper fit','original sqrt','Jasper sqrt','sqrt 0 intercept'},'Location','NorthWest');
    set(h_legend,'FontSize',16);
    xlim([0.3 0.55]);
    ylim([0.24 0.36]);
end
print([folder_Plots,'UstThr_D50_bins.png'],'-dpng');


%% PLOT FLUX VS U* and NON-DIMENSIONALIZED Q VS N-D U*, BASED ON FIT ABOVE.  ALSO MAKE PLOTS FOR TAU

%initialize plots
figure(7); clf; hold on;
figure(8); clf; hold on;
figure(9); clf; hold on;
figure(10); clf; hold on;
figure(11); clf; hold on;
figure(12); clf; hold on;
figure(13); clf; hold on;
figure(14); clf; hold on;
figure(15); clf; hold on;
figure(16); clf; hold on;

%initialize full lists
Q_all = [];
Q_nondim_all = [];
Q_nondim2_all = [];
Q_nondim3_all = [];
ust_nondim_all = [];
tau_nondim_all = [];
ust_ex_all = [];
tau_ex_all = [];

for i = 1:N_Sites
    %get values from lists
    Q_Site = Q_list_subset{i};
    d50_Site = d50_list_subset{i};
    ust_Site = ustRe_cal_list_subset{i};
    tau_Site = tauRe_cal_list_subset{i};
    
    %perform non-dimensionalization and excess calculations
    Q_nondim_Site = Q_Site./sqrt((rho_s/rho_a)*g*d50_Site.^3);
    ust_nondim_Site = ust_Site./polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_Site));
    tau_nondim_Site = ust_nondim_Site.^2;
    ust_ex_Site = ust_Site - polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_Site));
    tau_ex_Site = tau_Site - rho_a*polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_Site)).^2;
    Q_nondim2_Site = Q_Site*g./(1000*rho_a*tau_ex_Site.*polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_Site)));
    Q_nondim3_Site = Q_Site*g./(1000*rho_a*tau_ex_Site.*ust_Site);
    
    %add to lists
    Q_all = [Q_all; Q_Site];
    Q_nondim_all = [Q_nondim_all; Q_nondim_Site];
    Q_nondim2_all = [Q_nondim2_all; Q_nondim2_Site];
    Q_nondim3_all = [Q_nondim3_all; Q_nondim3_Site];
    ust_nondim_all = [ust_nondim_all; ust_nondim_Site];
    tau_nondim_all = [tau_nondim_all; tau_nondim_Site];
    ust_ex_all = [ust_ex_all; ust_ex_Site];
    tau_ex_all = [tau_ex_all; tau_ex_Site];
    
    %add to plots
    figure(7);
    plot(ust_Site,Q_Site,Markers{i},'MarkerSize',2);
    figure(8);
    plot(ust_nondim_Site,Q_nondim_Site,Markers{i},'MarkerSize',2);
    figure(9);
    plot(ust_nondim_Site,Q_nondim2_Site,Markers{i},'MarkerSize',2);
    figure(10);
    plot(ust_nondim_Site,Q_nondim3_Site,Markers{i},'MarkerSize',2);
    figure(11);
    plot(ust_ex_Site,Q_nondim_Site,Markers{i},'MarkerSize',2);
    figure(12);
    plot(ust_ex_Site,Q_Site,Markers{i},'MarkerSize',2);
    figure(13);
    plot(tau_Site,Q_Site,Markers{i},'MarkerSize',2);
    figure(14);
    plot(tau_nondim_Site,Q_nondim_Site,Markers{i},'MarkerSize',2);
    figure(15);
    plot(tau_ex_Site,Q_nondim_Site,Markers{i},'MarkerSize',2);
    figure(16);
    plot(tau_ex_Site,Q_Site,Markers{i},'MarkerSize',2);
end

%get Namikas / Greeley values, add to plots
figure(7);
plot(ust_mpers_Greeley96,Q_gperm2s_Greeley96,'k^');
plot(ust_mpers_Namikas03,Q_gperm2s_Namikas03,'kd');

figure(8);
ust_nondim_Greeley96 = ust_mpers_Greeley96/polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Greeley96));
ust_nondim_Namikas03 = ust_mpers_Namikas03/polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Namikas03));
Q_nondim_Greeley96 = Q_gperm2s_Greeley96./sqrt((rho_s/rho_a)*g*d50_mm_Greeley96.^3);
Q_nondim_Namikas03 = Q_gperm2s_Namikas03./sqrt((rho_s/rho_a)*g*d50_mm_Namikas03.^3);
plot(ust_nondim_Greeley96,Q_nondim_Greeley96,'k^');
plot(ust_nondim_Namikas03,Q_nondim_Namikas03,'kd');

figure(11);
ust_ex_Greeley96 = ust_mpers_Greeley96-polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Greeley96));
ust_ex_Namikas03 = ust_mpers_Namikas03-polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Namikas03));
plot(ust_ex_Greeley96,Q_nondim_Greeley96,'k^');
plot(ust_ex_Namikas03,Q_nondim_Namikas03,'kd');

figure(12);
plot(ust_ex_Greeley96,Q_gperm2s_Greeley96,'k^');
plot(ust_ex_Namikas03,Q_gperm2s_Namikas03,'kd');

figure(13);
tau_Greeley96 = rho_a*ust_mpers_Greeley96.^2;
tau_Namikas03 = rho_a*ust_mpers_Namikas03.^2;
plot(tau_Greeley96,Q_gperm2s_Greeley96,'k^');
plot(tau_Namikas03,Q_gperm2s_Namikas03,'kd');

figure(14);
tau_nondim_Greeley96 = ust_nondim_Greeley96.^2;
tau_nondim_Namikas03 = ust_nondim_Namikas03.^2;
plot(tau_nondim_Greeley96,Q_nondim_Greeley96,'k^');
plot(tau_nondim_Namikas03,Q_nondim_Namikas03,'kd');

figure(15);
tau_ex_Greeley96 = tau_Greeley96-rho_a.*polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Greeley96))^2;
tau_ex_Namikas03 = tau_Namikas03-rho_a.*polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Namikas03))^2;
plot(tau_ex_Greeley96,Q_nondim_Greeley96,'k^');
plot(tau_ex_Namikas03,Q_nondim_Namikas03,'kd');

figure(16);
plot(tau_ex_Greeley96,Q_gperm2s_Greeley96,'k^');
plot(tau_ex_Namikas03,Q_gperm2s_Namikas03,'kd');

%add fits to plots
ind_positive = find(ust_nondim_all>1);

%Generate fits and determine correlation coefficients
%dimensionless Q versus dimensionless tau
C_Q_nondim_tau_nondim_fit = mean(Q_nondim_all(ind_positive)./(tau_nondim_all(ind_positive)-1));
R_Q_nondim_tau_nondim = corrcoef(Q_nondim_all(ind_positive),tau_nondim_all(ind_positive)-1);
R_Q_nondim_tau_nondim = R_Q_nondim_tau_nondim(2);

%dimensionless Q versus dimensionless u* (squared relation)
C_Q_nondim2_ust_nondim_fit = median(Q_nondim2_all(ind_positive));

%dimensionless Q versus dimensionless u* (cubed relation)
C_Q_nondim3_ust_nondim_fit = median(Q_nondim3_all(ind_positive));

%bin values for fits
ust_nondim_bins_min = 1:.05:1.85;
ust_nondim_bins_max = 1.05:.05:1.9;
ust_nondim_bins_mid = mean([ust_nondim_bins_min;ust_nondim_bins_max]);
N_bins = length(ust_nondim_bins_mid);

C_Q_nondim2_binned = zeros(N_bins,1);
C_Q_nondim3_binned = zeros(N_bins,1);
for i=1:N_bins
    bin_ind = find(ust_nondim_all>=ust_nondim_bins_min(i)&ust_nondim_all<=ust_nondim_bins_max(i));
    C_Q_nondim2_binned(i) = median(Q_nondim2_all(bin_ind));
    C_Q_nondim3_binned(i) = median(Q_nondim3_all(bin_ind));
end

%dimensionless Q versus dimensionless tau
C_Q_nondim_tau_nondim_fit = mean(Q_nondim_all(ind_positive)./(tau_nondim_all(ind_positive)-1));
R_Q_nondim_tau_nondim = corrcoef(Q_nondim_all(ind_positive),tau_nondim_all(ind_positive)-1);
R_Q_nondim_tau_nondim = R_Q_nondim_tau_nondim(2);

%dimensionless Q versus tau excess
C_Q_nondim_tau_ex_fit = mean(Q_nondim_all(ind_positive)./(tau_ex_all(ind_positive)));
R_Q_nondim_tau_ex = corrcoef(Q_nondim_all(ind_positive),tau_ex_all(ind_positive));
R_Q_nondim_tau_ex = R_Q_nondim_tau_ex(2);

%dimensionless Q versus ust excess
C_Q_nondim_ust_ex_fit = mean(Q_nondim_all(ind_positive)./(ust_ex_all(ind_positive)));
R_Q_nondim_ust_ex = corrcoef(Q_nondim_all(ind_positive),ust_ex_all(ind_positive));
R_Q_nondim_ust_ex = R_Q_nondim_ust_ex(2);

%dimensional Q versus tau excess
C_Q_ust_ex_fit = mean(Q_all(ind_positive)./(ust_ex_all(ind_positive)));
R_Q_ust_ex = corrcoef(Q_all(ind_positive),(ust_ex_all(ind_positive)));
R_Q_ust_ex = R_Q_ust_ex(2);
C_Q_tau_ex_fit = mean(Q_all(ind_positive)./(tau_ex_all(ind_positive)));
R_Q_tau_ex = corrcoef(Q_all(ind_positive),(tau_ex_all(ind_positive)));
R_Q_tau_ex = R_Q_tau_ex(2);

%get values for plotting fit
tau_nondim_fit = linspace(1, max(tau_nondim_all),50);
tau_ex_fit = linspace(0, max(tau_ex_all),50);
ust_nondim_fit = sqrt(tau_nondim_fit);
ust_ex_fit = linspace(0, max(ust_ex_all),50);
Q_nondim_tau_nondim_fit = C_Q_nondim_tau_nondim_fit*(tau_nondim_fit-1);
Q_nondim_tau_ex_fit = C_Q_nondim_tau_ex_fit*tau_ex_fit;
Q_nondim_ust_ex_fit = C_Q_nondim_ust_ex_fit*ust_ex_fit;
Q_ust_ex_fit = C_Q_ust_ex_fit*ust_ex_fit;
Q_tau_ex_fit = C_Q_tau_ex_fit*tau_ex_fit;

%add fits to plots
figure(8);
plot(ust_nondim_fit,Q_nondim_tau_nondim_fit,'k','LineWidth',2);
figure(9);
plot(xlim,C_Q_nondim2_ust_nondim_fit*[1 1],'k','LineWidth',2);
plot(ust_nondim_bins_mid,C_Q_nondim2_binned,'k--','LineWidth',2);
figure(10);
plot(xlim,C_Q_nondim3_ust_nondim_fit*[1 1],'k','LineWidth',2);
plot(ust_nondim_bins_mid,C_Q_nondim3_binned,'k--','LineWidth',2);
figure(11);
plot(ust_ex_fit,Q_nondim_ust_ex_fit,'k','LineWidth',2);
figure(12);
plot(ust_ex_fit,Q_ust_ex_fit,'k','LineWidth',2);
figure(14);
plot(tau_nondim_fit,Q_nondim_tau_nondim_fit,'k','LineWidth',2);
figure(15);
plot(tau_ex_fit,Q_nondim_tau_ex_fit,'k','LineWidth',2);
figure(16);
plot(tau_ex_fit,Q_tau_ex_fit,'k','LineWidth',2);

%label plots
figure(7);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'Flux_Ust_all.png'],'-dpng');

figure(8);
xlabel('u_{*}/u_{*,th}','FontSize',16);
ylabel('Q(\rho_{s}gD/\rho_{a})^{-1/2}D^{-1}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_nondim_tau_nondim)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlim([0 2.5]);
ylim([0 1]);
print([folder_Plots,'FluxNondim_UstNonDim_all.png'],'-dpng');

figure(9);
xlabel('u_{*}/u_{*,th}','FontSize',16);
ylabel('gQ/{(\rho_{a}\tau_{ex}u_{*th})}','FontSize',16);
xlim([1 max(ust_nondim_all)]);
ylim([0 10]);
title('Squared law');
set(gca,'FontSize',16);
set(gca,'xscale','log');
legend_values = Sites;
legend_values{length(Sites)+1} = 'fit';
legend_values{length(Sites)+2} = 'binned fit';
h_legend = legend(legend_values,'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxNondim2_UstNonDim_all.png'],'-dpng');

figure(10);
xlabel('u_{*}/u_{*,th}','FontSize',16);
ylabel('gQ/{(\rho_{a}\tau_{ex}u_{*})}','FontSize',16);
xlim([1 max(ust_nondim_all)]);
ylim([0 10]);
title('Cubed law');
set(gca,'FontSize',16);
set(gca,'xscale','log');
legend_values = Sites;
legend_values{length(Sites)+1} = 'fit';
legend_values{length(Sites)+2} = 'binned fit';
h_legend = legend(legend_values,'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxNondim3_UstNonDim_all.png'],'-dpng');

figure(11);
xlabel('u_{*} - u_{*,th} (m/s)','FontSize',16);
ylabel('Q(\rho_{s}gD/\rho_{a})^{-1/2}D^{-1}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_nondim_ust_ex)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlim([-0.1 0.3]);
ylim([0 1]);
print([folder_Plots,'FluxNondim_UstEx_all.png'],'-dpng');

figure(12);
xlabel('u_{*} - u_{*,th} (m/s)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_ust_ex)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'Flux_UstEx_all.png'],'-dpng');

figure(13);
xlabel('\tau (Pa)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'Flux_Tau_all.png'],'-dpng');

figure(14);
xlabel('\tau/\tau_{th}','FontSize',16);
ylabel('Q(\rho_{s}gD/\rho_{a})^{-1/2}D^{-1}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_nondim_tau_nondim)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlim([0 6]);
ylim([0 1]);
print([folder_Plots,'FluxNondim_TauNonDim_all.png'],'-dpng');

figure(15);
xlabel('\tau - \tau_{th} (Pa)','FontSize',16);
ylabel('Q(\rho_{s}gD/\rho_{a})^{-1/2}D^{-1}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_nondim_tau_ex)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
xlim([-0.05 0.3]);
ylim([0 1]);
print([folder_Plots,'FluxNondim_TauEx_all.png'],'-dpng');

figure(16);
xlabel('\tau - \tau_{th} (Pa)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
legend_values{length(Sites)+3} = ['fit, R= ',num2str(R_Q_tau_ex)];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'Flux_TauEx_all.png'],'-dpng');


%% PLOT RESULTS - SALTATION HEIGHT

%% z_salt versus u*
figure(17); clf; hold on;
ust_bins_min = 0.15:0.05:0.55;
ust_bins_max = 0.2:0.05:0.6;
ust_bins_mid = mean([ust_bins_min; ust_bins_max]);
N_bins = length(ust_bins_mid);

for i=1:N_Sites
    zbar_bins_mean = zeros(N_bins,1);
    zbar_bins_SE = zeros(N_bins,1);
    Q_positive_ind = find(Q_list_subset{i}>0);
    for j = 1:N_bins
        ust_bin_ind = find(ustRe_cal_list_subset{i}>=ust_bins_min(j)&ustRe_cal_list_subset{i}<=ust_bins_max(j));
        zbar_bin_all = zbar_list_subset{i}(intersect(Q_positive_ind,ust_bin_ind));
        not_outlier_ind = find(zbar_bin_all>0&zbar_bin_all<0.30);
        zbar_bins_mean(j) = mean(zbar_bin_all(not_outlier_ind));
        zbar_bins_SE(j) = std(zbar_bin_all(not_outlier_ind))./sqrt(length(not_outlier_ind));
    end
    errorbar(ust_bins_mid,zbar_bins_mean,zbar_bins_SE,Markers{i});
end
%add Greeley and Namikas data
plot(ust_mpers_Greeley96,zbar_m_Greeley96,'^k');
plot(ust_mpers_Namikas03,zbar_m_Namikas03,'dk');

%label plot
xlabel('u_{*,Re} (m/s)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
xlim([0 0.65]);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightUst_Wenglor.png'],'-dpng');


%% saltation height versus d_surface
figure(18); clf; hold on;
for i=1:N_Sites
    d50_bins = unique(d50_list_subset{i});
    N_bins = length(d50_bins);
    zbar_bins_mean = zeros(N_bins,1);
    zbar_bins_SE = zeros(N_bins,1);
    Q_positive_ind = find(Q_list_subset{i}>0);
    for j = 1:N_bins
        d50_bin_ind = find(d50_list_subset{i}==d50_bins(j));
        zbar_bin_all = zbar_list_subset{i}(intersect(Q_positive_ind,d50_bin_ind));
        not_outlier_ind = find(zbar_bin_all>0&zbar_bin_all<0.30);
        zbar_bins_mean(j) = mean(zbar_bin_all(not_outlier_ind));
        zbar_bins_SE(j) = std(zbar_bin_all(not_outlier_ind))./sqrt(length(not_outlier_ind));
    end
    errorbar(d50_bins,zbar_bins_mean,zbar_bins_SE,Markers{i});
end
%add Greeley and Namikas data
zbar_mean_Greeley96 = mean(zbar_m_Greeley96);
zbar_SE_Greeley96 = std(zbar_m_Greeley96)/sqrt(N_Greeley96);
errorbar(d50_mm_Greeley96,zbar_mean_Greeley96,zbar_SE_Greeley96,'^k');
zbar_mean_Namikas03 = mean(zbar_m_Namikas03);
zbar_SE_Namikas03 = std(zbar_m_Namikas03)/sqrt(N_Namikas03);
errorbar(d50_mm_Namikas03,zbar_mean_Namikas03,zbar_SE_Namikas03,'dk');

%label plot
xlabel('d_{50} (mm)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightD50_Wenglor.png'],'-dpng');


%% saltation height / d_surface versus u*/u*th
figure(19); clf; hold on;
ust_nondim_bins_min = 0.5:0.25:1.75;
ust_nondim_bins_max = 0.75:0.25:2;
ust_nondim_bins_mid = mean([ust_nondim_bins_min; ust_nondim_bins_max]);
N_bins = length(ust_nondim_bins_mid);
for i = 1:N_Sites;
    zbar_nondim_bins_mean = zeros(N_bins,1);
    zbar_nondim_bins_SE = zeros(N_bins,1);
    Q_positive_ind = find(Q_list_subset{i}>0);
    ust_nondim = ustRe_cal_list_subset{i}./polyval(Pfit_ustThr_d50_linear,d50_list_subset{i});
    zbar_nondim = 1000*zbar_list_subset{i}./d50_list_subset{i};
    for j = 1:N_bins
        ust_bin_ind = find(ust_nondim>=ust_nondim_bins_min(j)&ust_nondim<=ust_nondim_bins_max(j));
        zbar_bin_all = zbar_nondim(intersect(Q_positive_ind,ust_bin_ind));
        not_outlier_ind = find(zbar_bin_all>0&zbar_bin_all<500);
        zbar_nondim_bins_mean(j) = mean(zbar_bin_all(not_outlier_ind));
        zbar_nondim_bins_SE(j) = std(zbar_bin_all(not_outlier_ind))/sqrt(length(not_outlier_ind));
    end
    errorbar(ust_nondim_bins_mid,zbar_nondim_bins_mean,zbar_nondim_bins_SE,Markers{i});
end
%add Greeley and Namikas data
ust_nondim_Greeley96 = ust_mpers_Greeley96/polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Greeley96));
zbar_nondim_Greeley96 = 1e3*zbar_m_Greeley96/d50_mm_Greeley96;
plot(ust_nondim_Greeley96,zbar_nondim_Greeley96,'^k');
ust_nondim_Namikas03 = ust_mpers_Namikas03/polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Namikas03));
zbar_nondim_Namikas03 = 1e3*zbar_m_Namikas03/d50_mm_Namikas03;
plot(ust_nondim_Namikas03,zbar_nondim_Namikas03,'dk');

%label plot
xlabel('u_{*,Re}/u_{*,thr}','FontSize',16);
ylabel('z_{salt}/d_{50}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightzOverD_Wenglor.png'],'-dpng');


%% saltation height / d_surface versus u*_ex
figure(20); clf; hold on;
ust_ex_bins_min = -0.05:0.025:0.225;
ust_ex_bins_max = -0.025:0.025:0.25;
ust_ex_bins_mid = mean([ust_ex_bins_min; ust_ex_bins_max]);
N_bins = length(ust_ex_bins_mid);
for i = 1:N_Sites;
    zbar_ex_bins_mean = zeros(N_bins,1);
    zbar_ex_bins_error = zeros(N_bins,1);
    Q_positive_ind = find(Q_list_subset{i}>0);
    ust_ex = ustRe_cal_list_subset{i}-polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_list_subset{i}));
    zbar_ex = 1000*zbar_list_subset{i}./d50_list_subset{i};
    for j = 1:N_bins
        ust_bin_ind = find(ust_ex>=ust_ex_bins_min(j)&ust_ex<=ust_ex_bins_max(j));
        zbar_bin_all = zbar_ex(intersect(Q_positive_ind,ust_bin_ind));
        not_outlier_ind = find(zbar_bin_all>0&zbar_bin_all<500);
        zbar_ex_bins_mean(j) = mean(zbar_bin_all(not_outlier_ind));
        zbar_ex_bins_SE(j) = std(zbar_bin_all(not_outlier_ind))/sqrt(length(not_outlier_ind));
    end
    errorbar(ust_ex_bins_mid,zbar_ex_bins_mean,zbar_ex_bins_SE,Markers{i},'MarkerSize',10);
end

%add Greeley and Namikas data
ust_ex_Greeley96 = ust_mpers_Greeley96-polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Greeley96));
zbar_nondim_Greeley96 = 1e3*zbar_m_Greeley96/d50_mm_Greeley96;
plot(ust_ex_Greeley96,zbar_nondim_Greeley96,'^k','MarkerSize',10);
ust_ex_Namikas03 = ust_mpers_Namikas03-polyval(Pfit_ustThr_d50_sqrt,sqrt(d50_mm_Namikas03));
zbar_nondim_Namikas03 = 1e3*zbar_m_Namikas03/d50_mm_Namikas03;
plot(ust_ex_Namikas03,zbar_nondim_Namikas03,'dk','MarkerSize',10);

%label plot
xlabel('u_{*,Re}-u_{*,thr} (m/s)','FontSize',16);
ylabel('z_{salt}/d_{50}','FontSize',16);
ylim([130 280]);
xlim([-0.05 0.35]);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 5])
print([folder_Plots,'FluxHeightzOverD_ustEx_Wenglor.png'],'-dpng');