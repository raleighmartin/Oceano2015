%% SCRIPT TO MAKE PLOTS FOR CALIBRATION OF DUST DATA

% DEPENDENT SCRIPTS USED HERE:
% CombineIntervals.m
% CreateTimeBlocks.m
% IntervalsWithin.m
% MaxOverlappingInterval.m

%% clear existing data and close all plots
clearvars;
close all;

%% add path for data analysis scripts -- this may need to be changed locally
folder_Scripts = '../../AeolianFieldworkAnalysis/Scripts/'; %folder with general scripts for data analysis
path(path,folder_Scripts); %add this path to locations for running functions

%% indicate where to save plots
folder_Plots = '../PlotOutput/Calibration/';

%% load processed and meta data
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/';
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_OceanoCalibration');
Metadata_Path = strcat(folder_ProcessedData,'Metadata_OceanoCalibration');
load(ProcessedData_Path); %load processed data 
load(Metadata_Path); %load metadata

%set time interval for computing calibration values
CalibrationTimeInterval = duration(0,15,0); %15 minutes

%set physical parameters
rho_a = 1.23; %kg/m^3

%set types of instruments for calibration
InstrumentTypes = {'Dust'};
N_InstrumentTypes = length(InstrumentTypes);

%get indices in Instrument Metadata of dust sensors
ind_Dust = find(strcmp(InstrumentMetadata.InstrumentType, 'Dust'));

%get list of all dust sensors
DustSensors = unique(InstrumentMetadata.Instrument(ind_Dust));
N_DustSensors = length(DustSensors);

%get start times and end times for all instruments
DustStartTimes = InstrumentMetadata.StartTime(ind_Dust);
DustEndTimes = InstrumentMetadata.EndTime(ind_Dust);

%combine time intervals based on start and end times for individual instruments
[CombinedStartTimes, CombinedEndTimes] = CombineIntervals(DustStartTimes, DustEndTimes);

%create time blocks based on combined time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(CombinedStartTimes, CombinedEndTimes, CalibrationTimeInterval);
N_Blocks = length(BlockStartTimes);

%initialize matrices of values - default value is NaN (so that it will be obvious if no data exists for this time interval)
n1_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n1 bin
n2_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n2 bin
n3_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n3 bin
n4_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n4 bin
n5_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n5 bin
n6_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n6 bin
n7_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n7 bin
n8_matrix = zeros(N_Blocks,N_DustSensors)*NaN; %n8 bin

%go through each block
for i = 1:N_Blocks
    BlockStartTime = BlockStartTimes(i);
    BlockEndTime = BlockEndTimes(i);
    
    for k = 1:N_DustSensors
        Sensor = DustSensors{k}; %set specific instrument

        %get indices for instrument values within block
        ind_DustIntervalStartsBefore = find([ProcessedData.Dust.(Sensor).StartTime]<=BlockStartTime); %indices of sensor time intervals starting before block times
        ind_DustIntervalEndsAfter = find([ProcessedData.Dust.(Sensor).EndTime]>=BlockEndTime); %indices of sensor time intervals ending after block times
        ind_DustIntervalsContainingBlock = intersect(ind_DustIntervalStartsBefore,ind_DustIntervalEndsAfter); %indices of sensor time intervals containing block times

        %get data only if interval is within block (otherwise it will remain NaN
        if ~isempty(ind_DustIntervalsContainingBlock)
            
            %extract data for this interval
            DustInterval = ProcessedData.Dust.(Sensor)(ind_DustIntervalsContainingBlock); %extract data for data interval containing time block
            ind_TimesInsideBlock = find((DustInterval.t.int>=BlockStartTime)&(DustInterval.t.int<=BlockEndTime)); %get indices within data interval of block

            %get mean dust concentration
            n1_bar = mean(DustInterval.n1.int(ind_TimesInsideBlock));
            n2_bar = mean(DustInterval.n2.int(ind_TimesInsideBlock));
            n3_bar = mean(DustInterval.n3.int(ind_TimesInsideBlock));
            n4_bar = mean(DustInterval.n4.int(ind_TimesInsideBlock));
            n5_bar = mean(DustInterval.n5.int(ind_TimesInsideBlock));
            n6_bar = mean(DustInterval.n6.int(ind_TimesInsideBlock));
            n7_bar = mean(DustInterval.n7.int(ind_TimesInsideBlock));
            n8_bar = mean(DustInterval.n8.int(ind_TimesInsideBlock));

            %add values to matrices
            n1_matrix(i,k) = n1_bar;
            n2_matrix(i,k) = n2_bar;
            n3_matrix(i,k) = n3_bar;
            n4_matrix(i,k) = n4_bar;
            n5_matrix(i,k) = n5_bar;
            n6_matrix(i,k) = n6_bar;
            n7_matrix(i,k) = n7_bar;
            n8_matrix(i,k) = n8_bar;
        end
    end
end

%get indices of D1 dust matrix for comparison
ind_D1 = find(strcmp(DustSensors,'D1')); %get column of matrix corresponding to D1

%go through each sensor for comparison to D1
for i = 1:N_DustSensors
    
    %get sensor
    Sensor = DustSensors{i};
   
    %get indices for comparison to D1
    ind_comparison = find(~isnan(n1_matrix(:,i)));

    %get dust concentrations for D1 and sensor of interest
    n1_D1 = n1_matrix(ind_comparison,ind_D1);
    n1_Sensor = n1_matrix(ind_comparison,i);
    n2_D1 = n2_matrix(ind_comparison,ind_D1);
    n2_Sensor = n2_matrix(ind_comparison,i);
    n3_D1 = n3_matrix(ind_comparison,ind_D1);
    n3_Sensor = n3_matrix(ind_comparison,i);
    n4_D1 = n4_matrix(ind_comparison,ind_D1);
    n4_Sensor = n4_matrix(ind_comparison,i);
    n5_D1 = n5_matrix(ind_comparison,ind_D1);
    n5_Sensor = n5_matrix(ind_comparison,i);
    n6_D1 = n6_matrix(ind_comparison,ind_D1);
    n6_Sensor = n6_matrix(ind_comparison,i);
    n7_D1 = n7_matrix(ind_comparison,ind_D1);
    n7_Sensor = n7_matrix(ind_comparison,i);
    n8_D1 = n8_matrix(ind_comparison,ind_D1);
    n8_Sensor = n8_matrix(ind_comparison,i);

    %compute simple fit
    cal_n1_simple = mean(n1_Sensor./n1_D1);
    cal_n2_simple = mean(n2_Sensor./n2_D1);
    cal_n3_simple = mean(n3_Sensor./n3_D1);
    cal_n4_simple = mean(n4_Sensor./n4_D1);
    cal_n5_simple = mean(n5_Sensor./n5_D1);
    cal_n6_simple = mean(n6_Sensor./n6_D1);
    cal_n7_simple = mean(n7_Sensor./n7_D1);
    cal_n8_simple = mean(n8_Sensor./n8_D1);

    %compute linear fit
    cal_n1_linearfit = polyfit(n1_D1,n1_Sensor,1);
    cal_n2_linearfit = polyfit(n2_D1,n2_Sensor,1);
    cal_n3_linearfit = polyfit(n3_D1,n3_Sensor,1);
    cal_n4_linearfit = polyfit(n4_D1,n4_Sensor,1);
    cal_n5_linearfit = polyfit(n5_D1,n5_Sensor,1);
    cal_n6_linearfit = polyfit(n6_D1,n6_Sensor,1);
    cal_n7_linearfit = polyfit(n7_D1,n7_Sensor,1);
    cal_n8_linearfit = polyfit(n8_D1,n8_Sensor,1);

    %get concentrations for simply calibrated Instrument
    n1_simplecal_Instrument = n1_Sensor/cal_n1_simple;
    n2_simplecal_Instrument = n2_Sensor/cal_n2_simple;
    n3_simplecal_Instrument = n3_Sensor/cal_n3_simple;
    n4_simplecal_Instrument = n4_Sensor/cal_n4_simple;
    n5_simplecal_Instrument = n5_Sensor/cal_n5_simple;
    n6_simplecal_Instrument = n6_Sensor/cal_n6_simple;
    n7_simplecal_Instrument = n7_Sensor/cal_n7_simple;
    n8_simplecal_Instrument = n8_Sensor/cal_n8_simple;

    %get concentration values for linear calibrated Instrument
    n1_linearcal_Instrument = (n1_Sensor-cal_n1_linearfit(2))/cal_n1_linearfit(1);
    n2_linearcal_Instrument = (n2_Sensor-cal_n2_linearfit(2))/cal_n2_linearfit(1);
    n3_linearcal_Instrument = (n3_Sensor-cal_n3_linearfit(2))/cal_n3_linearfit(1);
    n4_linearcal_Instrument = (n4_Sensor-cal_n4_linearfit(2))/cal_n4_linearfit(1);
    n5_linearcal_Instrument = (n5_Sensor-cal_n5_linearfit(2))/cal_n5_linearfit(1);
    n6_linearcal_Instrument = (n6_Sensor-cal_n6_linearfit(2))/cal_n6_linearfit(1);
    n7_linearcal_Instrument = (n7_Sensor-cal_n7_linearfit(2))/cal_n7_linearfit(1);
    n8_linearcal_Instrument = (n8_Sensor-cal_n8_linearfit(2))/cal_n8_linearfit(1);

    %compute relative error for dust concentration
    n1_nocal_relerr = mean(abs((n1_Sensor-n1_D1)./(n1_D1)));
    n1_simplecal_relerr = mean(abs((n1_simplecal_Instrument-n1_D1)./(n1_D1)));
    n1_linearcal_relerr = mean(abs((n1_linearcal_Instrument-n1_D1)./(n1_D1)));

    n2_nocal_relerr = mean(abs((n2_Sensor-n2_D1)./(n2_D1)));
    n2_simplecal_relerr = mean(abs((n2_simplecal_Instrument-n2_D1)./(n2_D1)));
    n2_linearcal_relerr = mean(abs((n2_linearcal_Instrument-n2_D1)./(n2_D1)));

    n3_nocal_relerr = mean(abs((n3_Sensor-n3_D1)./(n3_D1)));
    n3_simplecal_relerr = mean(abs((n3_simplecal_Instrument-n3_D1)./(n3_D1)));
    n3_linearcal_relerr = mean(abs((n3_linearcal_Instrument-n3_D1)./(n3_D1)));

    n4_nocal_relerr = mean(abs((n4_Sensor-n4_D1)./(n4_D1)));
    n4_simplecal_relerr = mean(abs((n4_simplecal_Instrument-n4_D1)./(n4_D1)));
    n4_linearcal_relerr = mean(abs((n4_linearcal_Instrument-n4_D1)./(n4_D1)));

    n5_nocal_relerr = mean(abs((n5_Sensor-n5_D1)./(n5_D1)));
    n5_simplecal_relerr = mean(abs((n5_simplecal_Instrument-n5_D1)./(n5_D1)));
    n5_linearcal_relerr = mean(abs((n5_linearcal_Instrument-n5_D1)./(n5_D1)));

    n6_nocal_relerr = mean(abs((n6_Sensor-n6_D1)./(n6_D1)));
    n6_simplecal_relerr = mean(abs((n6_simplecal_Instrument-n6_D1)./(n6_D1)));
    n6_linearcal_relerr = mean(abs((n6_linearcal_Instrument-n6_D1)./(n6_D1)));

    n7_nocal_relerr = mean(abs((n7_Sensor-n7_D1)./(n7_D1)));
    n7_simplecal_relerr = mean(abs((n7_simplecal_Instrument-n7_D1)./(n7_D1)));
    n7_linearcal_relerr = mean(abs((n7_linearcal_Instrument-n7_D1)./(n7_D1)));

    n8_nocal_relerr = mean(abs((n8_Sensor-n8_D1)./(n8_D1)));
    n8_simplecal_relerr = mean(abs((n8_simplecal_Instrument-n8_D1)./(n8_D1)));
    n8_linearcal_relerr = mean(abs((n8_linearcal_Instrument-n8_D1)./(n8_D1)));

    %plot dust concentrations - n1
    figure(10*i+1); clf;
    plot(n1_D1,n1_Sensor,'o'); hold on;
    plot([0 max(n1_D1)],cal_n1_simple*[0 max(n1_D1)]); %plot simple fit calibration
    plot([0 max(n1_D1)],polyval(cal_n1_linearfit,[0 max(n1_D1)])); %plot linear fit calibration
    xlabel('n1_{D1}','FontSize',14);
    ylabel(['n1_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n1_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n1_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n1_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n1_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n1_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n2
    figure(10*i+2); clf;
    plot(n2_D1,n2_Sensor,'o'); hold on;
    plot([0 max(n2_D1)],cal_n2_simple*[0 max(n2_D1)]); %plot simple fit calibration
    plot([0 max(n2_D1)],polyval(cal_n2_linearfit,[0 max(n2_D1)])); %plot linear fit calibration
    xlabel('n2_{D1}','FontSize',14);
    ylabel(['n2_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n2_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n2_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n2_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n2_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n2_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n3
    figure(10*i+3); clf;
    plot(n3_D1,n3_Sensor,'o'); hold on;
    plot([0 max(n3_D1)],cal_n3_simple*[0 max(n3_D1)]); %plot simple fit calibration
    plot([0 max(n3_D1)],polyval(cal_n3_linearfit,[0 max(n3_D1)])); %plot linear fit calibration
    xlabel('n3_{D1}','FontSize',14);
    ylabel(['n3_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n3_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n3_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n3_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n3_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n3_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n4
    figure(10*i+4); clf;
    plot(n4_D1,n4_Sensor,'o'); hold on;
    plot([0 max(n4_D1)],cal_n4_simple*[0 max(n4_D1)]); %plot simple fit calibration
    plot([0 max(n4_D1)],polyval(cal_n4_linearfit,[0 max(n4_D1)])); %plot linear fit calibration
    xlabel('n4_{D1}','FontSize',14);
    ylabel(['n4_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n4_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n4_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n4_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n4_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n4_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n5
    figure(10*i+5); clf;
    plot(n5_D1,n5_Sensor,'o'); hold on;
    plot([0 max(n5_D1)],cal_n5_simple*[0 max(n5_D1)]); %plot simple fit calibration
    plot([0 max(n5_D1)],polyval(cal_n5_linearfit,[0 max(n5_D1)])); %plot linear fit calibration
    xlabel('n5_{D1}','FontSize',14);
    ylabel(['n5_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n5_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n5_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n5_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n5_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n5_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n6
    figure(10*i+6); clf;
    plot(n6_D1,n6_Sensor,'o'); hold on;
    plot([0 max(n6_D1)],cal_n6_simple*[0 max(n6_D1)]); %plot simple fit calibration
    plot([0 max(n6_D1)],polyval(cal_n6_linearfit,[0 max(n6_D1)])); %plot linear fit calibration
    xlabel('n6_{D1}','FontSize',14);
    ylabel(['n6_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n6_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n6_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n6_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n6_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n6_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n7
    figure(10*i+7); clf;
    plot(n7_D1,n7_Sensor,'o'); hold on;
    plot([0 max(n7_D1)],cal_n7_simple*[0 max(n7_D1)]); %plot simple fit calibration
    plot([0 max(n7_D1)],polyval(cal_n7_linearfit,[0 max(n7_D1)])); %plot linear fit calibration
    xlabel('n7_{D1}','FontSize',14);
    ylabel(['n7_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n7_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n7_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n7_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n7_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n7_Calibration_D1_',Sensor,'.png'],'-dpng');

    %plot dust concentrations - n8
    figure(10*i+8); clf;
    plot(n8_D1,n8_Sensor,'o'); hold on;
    plot([0 max(n8_D1)],cal_n8_simple*[0 max(n8_D1)]); %plot simple fit calibration
    plot([0 max(n8_D1)],polyval(cal_n8_linearfit,[0 max(n8_D1)])); %plot linear fit calibration
    xlabel('n8_{D1}','FontSize',14);
    ylabel(['n8_{,',Sensor,'} (Pa)'],'FontSize',14);
    title(strcat('Cal factor = ',num2str(cal_n8_simple)),'FontSize',14);
    h_legend = legend(...
        ['data, \epsilon = ',num2str(n8_nocal_relerr)],...
        ['simple fit, \epsilon =',num2str(n8_simplecal_relerr)],...
        ['linear fit, \epsilon =',num2str(n8_linearcal_relerr)],...
        'Location','NorthWest');
    set(gca,'FontSize',14);
    set(h_legend,'FontSize',14);
    print([folder_Plots,'n8_Calibration_D1_',Sensor,'.png'],'-dpng');
end

%% Restore function path to default value
restoredefaultpath;