%% SCRIPT TO GENERATE PROFILES OF FLUX DATA AND COMPARE TO SHEAR VELOCITIES
% Dependencies: 'ExtractVariableTimeInterval', 'CreateTimeBlocks',
% 'reorient_anemometers_vanboxel2004', 'xcorrelation'

%%clear existing data and load processed data and metadata
clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationLags/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to processed data
load(Metadata_Path); %load meta data
load(ProcessedData_Path); %load processed data

%set physical parameters
rho_a = 1.23; %air density kg/m^3
g = 9.8; %gravity m/s^2

%set time interval for computing velocity profiles
AnalysisTimeInterval = duration(0,30,0); %30 minutes

%set minimum number of Wenglor counts
MinimumCountsForInterval = seconds(AnalysisTimeInterval)/10; %must average at least one per 10 seconds

%set reference instrument for lag calculations
ReferenceInstrument = 'S1';
ReferenceType = 'Sonic';

%set dates for comparison (some had errors)
DatesForComparison = [...
    datetime(2015,5,24);...
    datetime(2015,5,27);...
    datetime(2015,5,28);...
    datetime(2015,6,1);...
    datetime(2015,6,2);...
    datetime(2015,6,3);...
    datetime(2015,6,4)...
    ];

%get indices of time intervals for reference instrument
ReferenceInstrumentDates = [ProcessedData.(ReferenceType).(ReferenceInstrument).Date];
[~,ind_ReferenceInstrument, ~] = intersect(ReferenceInstrumentDates,DatesForComparison);

%get start times and end times for reference instrument
ReferenceStartTimes = [ProcessedData.(ReferenceType).(ReferenceInstrument)(ind_ReferenceInstrument).StartTime];
ReferenceEndTimes = [ProcessedData.(ReferenceType).(ReferenceInstrument)(ind_ReferenceInstrument).EndTime];

%create time blocks based on time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(ReferenceStartTimes, ReferenceEndTimes, AnalysisTimeInterval);
N_Blocks = length(BlockStartTimes);

% %intialize quantities
% Wenglors = fieldnames(ProcessedData.Wenglor); %get list of all Wenglors for analysis
% Sonics = fieldnames(ProcessedData.Sonic); %get list of all Sonics for analysis
% N_Wenglors = length(Wenglors);
% N_Sonics = length(Sonics);
% lags_Wenglors = zeros(N_Blocks,N_Wenglors)*NaN;
% lags_Sonics = zeros(N_Blocks,N_Sonics)*NaN;
% z_Wenglors = zeros(N_Blocks,N_Wenglors)*NaN;
% z_Sonics = zeros(N_Blocks,N_Sonics)*NaN;
% qz_Wenglors = zeros(N_Blocks,N_Wenglors)*NaN;
% Q_Wenglors = zeros(N_Blocks)*NaN;
% theta_Reference = zeros(N_Blocks)*NaN;
% ustRe_Reference = zeros(N_Blocks)*NaN;
% 
% 
% %go through blocks to get lags
% for i = 1:N_Blocks
%     
%     %print out progress of processing
%     Progress = [int2str(i),' of ',int2str(N_Blocks)]
%     
%     %get values for reference instrument
%     [u_Reference, t_Reference] = ExtractVariableTimeInterval(ProcessedData.(ReferenceType).(ReferenceInstrument),BlockStartTimes(i),BlockEndTimes(i),'u','int','int');
%     v_Reference = ExtractVariableTimeInterval(ProcessedData.(ReferenceType).(ReferenceInstrument),BlockStartTimes(i),BlockEndTimes(i),'v','int','int');
%     w_Reference = ExtractVariableTimeInterval(ProcessedData.(ReferenceType).(ReferenceInstrument),BlockStartTimes(i),BlockEndTimes(i),'w','int','int');
%     theta_Reference(i) = atan(mean(v_Reference)/mean(u_Reference))*180/pi; %calculate angle of unrotated instrument, add to list
%     [u_Reference, v_Reference, w_Reference] = reorient_anemometers_vanboxel2004(u_Reference, v_Reference, w_Reference); %rotate instrument
%     u_bar_Reference = mean(u_Reference);
%     w_bar_Reference = mean(w_Reference);
%     tauRe_Reference = -rho_a*mean((u_Reference-u_bar_Reference).*(w_Reference-w_bar_Reference)); %compute Reynolds stress
%     ustRe_Reference(i) = sqrt(tauRe_Reference/rho_a); %compute u*, add to list                
%     
%     %go through Wenglors
%     for j = 1:N_Wenglors
%         
%         %get information for Wenglor during interval
%         [n_Wenglor, t_Wenglor, interval_Wenglor, ~] = ExtractVariableTimeInterval(ProcessedData.Wenglor.(Wenglors{j}),BlockStartTimes(i),BlockEndTimes(i),'n','int','int');
%         
%         %check to make sure times fill whole interval window and that
%         %Wenglor does not have error code and that significant saltation
%         %was detected
%         if ((~isempty(t_Wenglor)&&length(interval_Wenglor)==1)&&...
%                 ((min(t_Wenglor)==BlockStartTimes(i))&&(max(t_Wenglor)==BlockEndTimes(i))))&&...
%                 ((ProcessedData.Wenglor.(Wenglors{j})(interval_Wenglor).ErrorCode==0)&&...
%                 sum(n_Wenglor)>=MinimumCountsForInterval)
%             
%             %now, make computations for cross-correlation
%             [~,ind_Reference,ind_Wenglor] = intersect(t_Reference,t_Wenglor); %get indices of intersecting times
%             dt = mode(diff(t_Reference(ind_Reference))); %get modal dt
%             [~, ~, T_max, ~] = xcorrelation(u_Reference(ind_Reference),n_Wenglor(ind_Wenglor),dt); %compute cross-correlation
%             lags_Wenglors(i,j) = seconds(T_max); %convert time into seconds
%             z_Wenglors(i,j) = ProcessedData.Wenglor.(Wenglors{j})(interval_Wenglor).InstrumentHeight.z; %get instrument height
%         
%             %get flux information
%             qz_Wenglors(i,j) = mean(ExtractVariableTimeInterval(ProcessedData.Wenglor.(Wenglors{j}),BlockStartTimes(i),BlockEndTimes(i),'flux','qz','int'));
%             
%         end
%     end
%     
%     %compute total flux
%     z_fit = z_Wenglors(i,~isnan(z_Wenglors(i,:)))
%     qz_fit = qz_Wenglors(i,~isnan(z_Wenglors(i,:)))
%     fit_params = polyfit(z_fit,log(qz_fit),1); %fit only to nonzero values of qz        
%     zbar = -1/fit_params(1); %m
%     q0 = exp(fit_params(2)); %g/m^2/s
%     Q = q0*zbar; %g/m/s
% 
%     if (~isinf(Q)&&~isnan(Q))&&(length(z_fit)>=2)
%         Q_Wenglors(i) = Q; %set to Q if polyfit is well-defined
%     else
%         Q_Wenglors(i) = 0; %set to zero is polyfit is ill-defined
%     end
%     
%     %go through Sonics
%     for j = 1:N_Sonics
%         
%         %get information for Sonic during interval
%         [u_Sonic, t_Sonic, interval_Sonic, ~] = ExtractVariableTimeInterval(ProcessedData.Sonic.(Sonics{j}),BlockStartTimes(i),BlockEndTimes(i),'u','int','int');
%         
%         %check to make sure times fill whole interval window
%         if (~isempty(t_Sonic)&&length(interval_Sonic)==1)&&((min(t_Sonic)==BlockStartTimes(i))&&(max(t_Sonic)==BlockEndTimes(i)))
%             
%             %get v, w, and perform rotation
%             v_Sonic = ExtractVariableTimeInterval(ProcessedData.Sonic.(Sonics{j}),BlockStartTimes(i),BlockEndTimes(i),'v','int','int');
%             w_Sonic = ExtractVariableTimeInterval(ProcessedData.Sonic.(Sonics{j}),BlockStartTimes(i),BlockEndTimes(i),'w','int','int');
%             [u_Sonic, v_Sonic, w_Sonic] = reorient_anemometers_vanboxel2004(u_Sonic, v_Sonic, w_Sonic); %rotate instrument
% 
%             %now, make computations for cross-correlation
%             [~,ind_Reference,ind_Sonic] = intersect(t_Reference,t_Sonic); %get indices of intersecting times
%             dt = mode(diff(t_Reference(ind_Reference))); %get modal dt
%             [~, ~, T_max, ~] = xcorrelation(u_Reference(ind_Reference),u_Sonic(ind_Sonic),dt); %compute cross-correlation
%             lags_Sonics(i,j) = seconds(T_max); %convert time into seconds
%             z_Sonics(i,j) = ProcessedData.Sonic.(Sonics{j})(interval_Sonic).InstrumentHeight.z; %get instrument height
%         
%         end
%         
%     end
%     
%     %plot lags
%     figure(1); clf;
%     plot(lags_Sonics(i,:), z_Sonics(i,:), 'o','MarkerSize',10); hold on;
%     plot(lags_Wenglors(i,:), z_Wenglors(i,:), 'x','MarkerSize',10);
%     xlabel('lag (s)');
%     ylabel('z (m)');
%     set(gca,'FontSize',16);
%     legend('Sonics','Wenglors','Location','NorthEast');
%     
%     title([datestr(BlockStartTimes(i),'yyyy-mm-dd HH:MM'),'-',...
%         datestr(BlockEndTimes(i),'HH:MM'),', ',...
%         'u_{*} = ',num2str(ustRe_Reference(i),2),' m/s, ',...
%         '\theta = ',num2str(theta_Reference(i),2),'\circ, ',...
%         'Q = ',num2str(Q_Wenglors(i),2),' g/m/s'],'FontSize',14);
%     
%     print([folder_Plots,'lags_',datestr(BlockStartTimes(i),'yyyy-mm-dd_HHMM'),'_linearscale.png'],'-dpng');
%     set(gca,'yscale','log');
%     print([folder_Plots,'lags_',datestr(BlockStartTimes(i),'yyyy-mm-dd_HHMM'),'_logscale.png'],'-dpng');
% end

load([folder_ProcessedData,'lag_computations_Oceano']);

median_lag_Wenglors = zeros(N_Blocks,1);
expected_lag_Sonics = zeros(N_Blocks,1);
for i = 1:N_Blocks
    median_lag_Wenglors(i) = median(lags_Wenglors(i,~isnan(lags_Wenglors(i,:))));
    fit_params = polyfit(z_Sonics(i,1:3),lags_Sonics(i,1:3),1);
    expected_lag_Sonics(i) = fit_params(2);
end

figure(2); clf;
semilogx(Q_Wenglors(~isnan(median_lag_Wenglors)),median_lag_Wenglors(~isnan(median_lag_Wenglors)),'o','MarkerSize',10)
xlabel('flux (g/m/s)','FontSize',16);
ylabel('median saltation lag (s)','FontSize',16);
ylim([0 5]);
set(gca,'FontSize',16);
title('Not Inclination Corrected');
print([folder_Plots,'AllLagsVersusFlux_NotInclinationCorrected.png'],'-dpng');

figure(3); clf;
semilogx(Q_Wenglors(~isnan(median_lag_Wenglors)),median_lag_Wenglors(~isnan(median_lag_Wenglors))-expected_lag_Sonics(~isnan(median_lag_Wenglors)),'o','MarkerSize',10)
xlabel('flux (g/m/s)','FontSize',16);
ylabel('median saltation lag (s)','FontSize',16);
ylim([-0.1 5]);
set(gca,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 4])
%title('Inclination Corrected');
print([folder_Plots,'AllLagsVersusFlux_InclinationCorrected.png'],'-dpng');

figure(4); clf;
plot(ustRe_Reference, expected_lag_Sonics,'o','MarkerSize',10)
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('inclination lag (s)','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'AllInclinationLags.png'],'-dpng');

% save([folder_ProcessedData,'lag_computations'],...
%     'BlockStartTimes','BlockEndTimes','ReferenceInstrument',...
%     'Wenglors','Sonics','lags_Wenglors','lags_Sonics',...
%     'z_Wenglors','z_Sonics','qz_Wenglors','Q_Wenglors',...
%     'theta_Reference','ustRe_Reference');