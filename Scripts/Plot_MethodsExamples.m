%% SCRIPT TO PLOT SAMPLE OF WENGLOR CALIBRATION VALUES

%%clear existing data and load processed data and metadata
%clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/WenglorCalibration/'; %folder for plots
saveMetadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
saveProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
load(saveMetadata_Path); %load metadata
%load(saveProcessedData_Path); %load processed data

PlotDate = datetime(2015,5,23);

%Wenglors = fieldnames(ProcessedData.Wenglor);
Wenglors = fieldnames(Data{3}.ProcessedData.Wenglor);
N_Wenglors = length(Wenglors);

%analysis of calibration factors
figure(1); clf; hold on;

for i = 1:N_Wenglors
%    [C, t, ~, ~] = ExtractVariableTimeInterval(ProcessedData.Wenglor.(Wenglors{i}),PlotDate,PlotDate+duration(24,0,0),'flux','qzPerCount','int');
    [C, t, ~, ~] = ExtractVariableTimeInterval(Data{3}.ProcessedData.Wenglor.(Wenglors{i}),PlotDate,PlotDate+duration(24,0,0),'flux','qzPerCount','int');
    plot(t,C,'LineWidth',3);
end

title(datestr(PlotDate));
ylabel('C_{qn} (g/m^2/s per count)');
legend(Wenglors,'Location','EastOutside');
set(gca,'FontSize',16);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 3])
print([folder_Plots,'Wenglor_Calibration_',datestr(PlotDate,'yyyy-mm-dd'),'.png'],'-dpng');

%analysis of thresholds
Wenglor_ThresholdAnalysis = 'W1';
Sonic_ThresholdAnalysis = 'S1';
window_average_steps = round(logspace(log10(1),log10(15000),20));
window_average_T = duration(0,0,0.04)*window_average_steps;
N_windows = length(window_average_steps);

[~, ~, int_N, ~] = ExtractVariableTimeInterval(ProcessedData.Wenglor.(Wenglor_ThresholdAnalysis),PlotDate,PlotDate+duration(24,0,0),'flux','qzPerCount','int');
Wenglor_StartTime = ProcessedData.Wenglor.(Wenglor_ThresholdAnalysis)(int_N).StartTime;
Wenglor_EndTime = ProcessedData.Wenglor.(Wenglor_ThresholdAnalysis)(int_N).EndTime;
t = ProcessedData.Wenglor.(Wenglor_ThresholdAnalysis)(int_N).t.int;
n = ProcessedData.Wenglor.(Wenglor_ThresholdAnalysis)(int_N).n.int;
[u, ~] = ExtractVariableTimeInterval(ProcessedData.Sonic.(Sonic_ThresholdAnalysis),Wenglor_StartTime,Wenglor_EndTime,'u','int','int');

counts_time_fraction = zeros(N_windows,1);
for i = 1:N_windows
    [n_avg, ~] = window_average(n, t, window_average_T(i));
    counts_time_fraction(i) = length(find(n_avg>0))/length(n_avg);
end

figure(2); clf;
%plot frequencies
subplot(1,2,1)
plot(seconds(window_average_T),counts_time_fraction,'o','MarkerSize',10)
ylim([0 1]);
set(gca,'xscale','log');
xlim([min(seconds(window_average_T)), max(seconds(window_average_T))]);
xlabel('averaging time (s)');
ylabel('f_{transport}');
title(['Wenglor counts: ', datestr(PlotDate)]);
set(gca,'FontSize',16);

%plot distribution of wind velocities
subplot(1,2,2)
plot(sort(u),1-linspace(0,1,length(u)),'LineWidth',3);
ylim([0 1]);
xlabel('u (m/s)');
ylabel('1-CDF');
title(['Wind velocities: ', datestr(PlotDate)]);
set(gca,'FontSize',16);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 3])
print([folder_Plots,'Wenglor_Threshold_',datestr(PlotDate,'yyyy-mm-dd'),'.png'],'-dpng');
