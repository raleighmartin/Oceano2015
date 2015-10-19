%% SCRIPT TO MAKE PLOTS FOR CALIBRATION OF DUST, WIND, AND WEATHER DATA

%% clear existing data and load processed data and metadata
clearvars;

%% add path for data analysis scripts -- this may need to be changed locally
folder_Scripts = '../../AeolianFieldworkAnalysis/Scripts/'; %folder with general scripts for data analysis
path(path,folder_Scripts); %add this path to locations for running functions
% DEPENDENT SCRIPTS USED HERE:
% CombineIntervals.m
% CreateTimeBlocks.m
% IntervalsWithin.m
% MaxOverlappingInterval.m

%% load processed and meta data
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/';
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_OceanoCalibration');
Metadata_Path = strcat(folder_ProcessedData,'Metadata_OceanoCalibration');
load(ProcessedData_Path); %load processed data 
load(Metadata_Path); %load metadata

%% indicate where to save plots
folder_Plots = '../PlotOutput/Calibration/';

%set time interval for computing calibration values
CalibrationTimeInterval = duration(0,15,0); %15 minutes

%set physical parameters
rho_a = 1.23; %kg/m^3

%set types of instruments for calibration
InstrumentTypes = {'Cup','Weather','Dust','Sonic','Ultrasonic'};
N_InstrumentTypes = length(InstrumentTypes);

%get indices in Instrument Metadata of these instrument types
ind_AllInstrumentTypes = [];
for i = 1:N_InstrumentTypes
    ind_AllInstrumentTypes = [ind_AllInstrumentTypes;...
        find(strcmp(InstrumentMetadata.InstrumentType, InstrumentTypes{i}))];
end

%get list of all instruments and associated types
AllInstruments = unique(InstrumentMetadata.Instrument(ind_AllInstrumentTypes));
N_AllInstruments = length(AllInstruments);
AllInstruments_Type = cell(N_AllInstruments,1);
for i = 1:N_AllInstruments
    AllInstruments_Type{i} = char(unique(InstrumentMetadata.InstrumentType(strcmp(InstrumentMetadata.Instrument,AllInstruments{i}))));
end

%get start times and end times for all instruments
InstrumentStartTimes = InstrumentMetadata.StartTime(ind_AllInstrumentTypes);
InstrumentEndTimes = InstrumentMetadata.EndTime(ind_AllInstrumentTypes);

%combine time intervals based on start and end times for individual instruments
[CombinedStartTimes, CombinedEndTimes] = CombineIntervals(InstrumentStartTimes, InstrumentEndTimes);

%create time blocks based on combined time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(CombinedStartTimes, CombinedEndTimes, CalibrationTimeInterval);
N_Blocks = length(BlockStartTimes);

%initialize matrices of values
u_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %streamwise velocity
theta_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %velocity angle
tauRe_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %Reynolds stress
heatflux_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %heat flux
T_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %temperatures

n1_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n1 bin
n2_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n2 bin
n3_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n3 bin
n4_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n4 bin
n5_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n5 bin
n6_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n6 bin
n7_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n7 bin
n8_matrix = zeros(N_Blocks,N_AllInstruments)*NaN; %n8 bin

%go through each block
for i = 1:N_Blocks
    BlockStartTime = BlockStartTimes(i);
    BlockEndTime = BlockEndTimes(i);
    
    for j = 1:N_InstrumentTypes
        InstrumentType = InstrumentTypes{j}; %set specific instrument type
        Instruments = fieldnames(ProcessedData.(InstrumentType));
        N_Instruments = length(Instruments);
        for k = 1:N_Instruments
            Instrument = Instruments{k}; %set specific instrument
            
            %get indices for instrument values within block
            ind_AllInstruments = find(strcmp(AllInstruments,Instrument));
            ind_InstrumentIntervalStartsBefore = find([ProcessedData.(InstrumentType).(Instrument).StartTime]<=BlockStartTime); %indices of instrument time intervals starting before block times
            ind_InstrumentIntervalEndsAfter = find([ProcessedData.(InstrumentType).(Instrument).EndTime]>=BlockEndTime); %indices of instrument time intervals ending after block times
            ind_InstrumentIntervalsContainingBlock = intersect(ind_InstrumentIntervalStartsBefore,ind_InstrumentIntervalEndsAfter); %indices of instrument time intervals containing block times
            
            %get data only if interval is within block
            if ~isempty(ind_InstrumentIntervalsContainingBlock)
                %extract data for this interval
                InstrumentInterval = ProcessedData.(InstrumentType).(Instrument)(ind_InstrumentIntervalsContainingBlock); %extract data for instrument time interval
                ind_TimesInsideBlock = find((InstrumentInterval.t.int>=BlockStartTime)&(InstrumentInterval.t.int<=BlockEndTime)); %get indices of times within block

                %calculations for cups
                if strcmp(InstrumentType,'Cup')
                    %get u
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    u_bar = mean(u);
                    
                    %add values to matrices
                    u_matrix(i,ind_AllInstruments) = u_bar;

                %calculations for sonics
                elseif strcmp(InstrumentType,'Ultrasonic')||strcmp(InstrumentType,'Sonic')

                    %get temperature
                    T = InstrumentInterval.T.int(ind_TimesInsideBlock);
                    T_bar = mean(T);

                    %get velocities
                    u = InstrumentInterval.u.int(ind_TimesInsideBlock);
                    v = InstrumentInterval.v.int(ind_TimesInsideBlock);
                    w = InstrumentInterval.w.int(ind_TimesInsideBlock);

                    %compute velocity angle and perform rotation
                    theta = atan(mean(v)/mean(u))*180/pi;
                    [u, v, w] = reorient_anemometers_vanboxel2004(u, v, w);

                    %compute u_bar and w_bar
                    u_bar = mean(u);
                    w_bar = mean(w);

                    %compute Reynolds stress and heat flux
                    tauRe = -rho_a*mean((u-u_bar).*(w-w_bar));
                    heatflux = mean((w-w_bar).*(T-T_bar));
                    
                    %add values to matrices
                    T_matrix(i,ind_AllInstruments) = T_bar;
                    u_matrix(i,ind_AllInstruments) = u_bar;
                    theta_matrix(i,ind_AllInstruments) = theta;
                    tauRe_matrix(i,ind_AllInstruments) = tauRe;
                    heatflux_matrix(i,ind_AllInstruments) = heatflux;
                    
                %calculations for weather station
                elseif strcmp(InstrumentType,'Weather')

                    %get temperature
                    T = InstrumentInterval.T.int(ind_TimesInsideBlock);
                    T_bar = mean(T);

                    %add values to matrices
                    T_matrix(i,ind_AllInstruments) = T_bar;

                %calculations for dust
                elseif strcmp(InstrumentType,'Dust')

                    %get dust concentration
                    n1_bar = mean(InstrumentInterval.n1.int(ind_TimesInsideBlock));
                    n2_bar = mean(InstrumentInterval.n2.int(ind_TimesInsideBlock));
                    n3_bar = mean(InstrumentInterval.n3.int(ind_TimesInsideBlock));
                    n4_bar = mean(InstrumentInterval.n4.int(ind_TimesInsideBlock));
                    n5_bar = mean(InstrumentInterval.n5.int(ind_TimesInsideBlock));
                    n6_bar = mean(InstrumentInterval.n6.int(ind_TimesInsideBlock));
                    n7_bar = mean(InstrumentInterval.n7.int(ind_TimesInsideBlock));
                    n8_bar = mean(InstrumentInterval.n8.int(ind_TimesInsideBlock));

                    %add values to matrices
                    n1_matrix(i,ind_AllInstruments) = n1_bar;
                    n2_matrix(i,ind_AllInstruments) = n2_bar;
                    n3_matrix(i,ind_AllInstruments) = n3_bar;
                    n4_matrix(i,ind_AllInstruments) = n4_bar;
                    n5_matrix(i,ind_AllInstruments) = n5_bar;
                    n6_matrix(i,ind_AllInstruments) = n6_bar;
                    n7_matrix(i,ind_AllInstruments) = n7_bar;
                    n8_matrix(i,ind_AllInstruments) = n8_bar;
                end
            end
        end
    end
end

%get indices of S1 sonic matrix for comparison
ind_S1 = find(strcmp(AllInstruments,'S1')); %get column of matrix corresponding to S1
ind_D1 = find(strcmp(AllInstruments,'D1')); %get column of matrix corresponding to D1
theta_S1 = theta_matrix(:,ind_S1); %get list of thetas for S1
ind_onshorewind = find(abs(theta_S1)<=45); %get indices only of those within +/- 45 degrees

%go through each instrument for comparison to S1 or D1
for i = 1:N_AllInstruments
    
    %get instrument
    Instrument = AllInstruments{i};
    
    %look at mean velocities
    if (strcmp(AllInstruments_Type(i),'Ultrasonic')||strcmp(AllInstruments_Type(i),'Sonic'))||strcmp(AllInstruments_Type(i),'Cup');

        %get indices for comparison to S1
        ind_comparison = intersect(find(~isnan(u_matrix(:,i))),ind_onshorewind);

        %get u-values for S1 and instrument
        u_S1 = u_matrix(ind_comparison,ind_S1);
        u_Instrument = u_matrix(ind_comparison,i); 

        %compute simple calibration factor and linear fit
        cal_u_simple = mean(u_Instrument./u_S1);
        cal_u_linearfit = polyfit(u_S1,u_Instrument,1);

        %get u-values for calibrated Instrument
        u_simplecal_Instrument = u_Instrument/cal_u_simple;
        u_linearcal_Instrument = (u_Instrument-cal_u_linearfit(2))/cal_u_linearfit(1);

        %compute relative error
        u_nocal_relerr = mean(abs((u_Instrument-u_S1)./(u_S1)));
        u_simplecal_relerr = mean(abs((u_simplecal_Instrument-u_S1)./(u_S1)));
        u_linearcal_relerr = mean(abs((u_linearcal_Instrument-u_S1)./(u_S1)));

        %plot comparison
        figure(i); clf;
        plot(u_S1,u_Instrument,'o'); hold on;
        plot([0 15],[0 cal_u_simple*15]); %plot simple calibration
        plot([0 15],polyval(cal_u_linearfit,[0 15])); %plot linear fit calibration
        xlabel('u_{S1} (m/s)','FontSize',14);
        ylabel(['u_{',Instrument,'} (m/s)'],'FontSize',14);
        if Instrument(1)=='C'
            title(strcat('Cal eq = ',num2str(cal_u_linearfit(1)),'u_{S1} + (',num2str(cal_u_linearfit(2)),')'),'FontSize',14);
        else
            title(strcat('Cal factor = ',num2str(cal_u_simple)),'FontSize',14);
        end
        set(gca,'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(u_nocal_relerr)],...
            ['simple calib., \epsilon = ',num2str(u_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(u_linearcal_relerr)],...
            'Location','NorthWest');
        set(h_legend,'FontSize',14);
        print([folder_Plots,'u_Calibration_S1_',Instrument,'.png'],'-dpng');
    end
    
    %look at mean temperatures
    if (strcmp(AllInstruments_Type(i),'Ultrasonic')||strcmp(AllInstruments_Type(i),'Sonic'))||strcmp(AllInstruments_Type(i),'Weather');
        %get indices for comparison to S1
        ind_comparison = intersect(find(~isnan(T_matrix(:,i))),ind_onshorewind);

        %get u-values for S1 and instrument
        T_S1 = T_matrix(ind_comparison,ind_S1);
        T_Instrument = T_matrix(ind_comparison,i); 

        %compute linear fit
        cal_T_linearfit = polyfit(T_S1,T_Instrument,1);

        %get T-values for calibrated Instrument
        T_linearcal_Instrument = (T_Instrument-cal_T_linearfit(2))/cal_T_linearfit(1);

        %compute relative error
        T_nocal_relerr = mean(abs((T_Instrument-T_S1)./(T_S1)));
        T_cal_relerr = mean(abs((T_linearcal_Instrument-T_S1)./(T_S1)));

        %plot comparison
        figure(i+N_AllInstruments); clf;
        plot(T_S1,T_Instrument,'o'); hold on;
        plot([14 16.5],polyval(cal_T_linearfit,[14 16.5])); %plot linear fit calibration
        xlabel('T_{S1} (m/s)','FontSize',14);
        ylabel(['T_{',Instrument,'} (m/s)'],'FontSize',14);
        title(strcat('Cal eq = ',num2str(cal_T_linearfit(1)),'T_{S1} + (',num2str(cal_T_linearfit(2)),')'),'FontSize',14);
        set(gca,'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(T_nocal_relerr)],...
            ['linear fit, \epsilon =',num2str(T_cal_relerr)],...
            'Location','NorthWest');
        set(h_legend,'FontSize',14);
        print([folder_Plots,'T_Calibration_S1_',Instrument,'.png'],'-dpng');
    end
    
    %look at Reynolds stresses and heat fluxes
    if strcmp(AllInstruments_Type(i),'Ultrasonic')||strcmp(AllInstruments_Type(i),'Sonic')
        
        %get indices for comparison to S1
        ind_comparison = intersect(find(~isnan(tauRe_matrix(:,i))),ind_onshorewind);
        
        %get Reynolds stresses for S1 and instrument
        tauRe_S1 = tauRe_matrix(ind_comparison,ind_S1);
        tauRe_Instrument = tauRe_matrix(ind_comparison,i);
        
        %get calibrated Reynolds stresses for instrument
        tauRe_cal_Instrument = tauRe_Instrument/cal_u_simple.^2;
        
        %compute relative error for Reynolds stress
        tauRe_nocal_relerr = mean(abs((tauRe_Instrument-tauRe_S1)./(tauRe_S1)));
        tauRe_cal_relerr = mean(abs((tauRe_cal_Instrument-tauRe_S1)./(tauRe_S1)));
        
        %get heat fluxes for S1 and instrument;
        heatflux_S1 = heatflux_matrix(ind_comparison,ind_S1);
        heatflux_Instrument = heatflux_matrix(ind_comparison,i);
        
        %get calibrated heat flux for instrument
        heatflux_cal_Instrument = heatflux_Instrument/(cal_u_simple*cal_T_linearfit(1));
        
        %compute relative error for heat flux
        heatflux_nocal_relerr = mean(abs((heatflux_Instrument-heatflux_S1)./(heatflux_S1)));
        heatflux_cal_relerr = mean(abs((heatflux_cal_Instrument-heatflux_S1)./(heatflux_S1)));
        
        %compute linear fits for Reynolds stress and heat flux
        tauRe_linearfit = polyfit(tauRe_S1,tauRe_Instrument,1);
        heatflux_linearfit = polyfit(heatflux_S1,heatflux_Instrument,1);
        
        %plot Reynolds stress
        figure(i+2*N_AllInstruments); clf;
        plot(tauRe_S1,tauRe_Instrument,'o'); hold on;
        plot([0 0.07],[0 0.07]); %plot 1-1 line
        plot([0 0.07],polyval(tauRe_linearfit,[0 0.07])); %plot linear fit
        xlabel('\tau_{Re,S1} (Pa)','FontSize',14);
        ylabel(['\tau_{Re,',Instrument,'} (Pa)'],'FontSize',14);
        set(gca,'FontSize',14);
        h_legend = legend('data','1-1 line',...
            ['fit: ',num2str(tauRe_linearfit(1)),'\tau_{Re,S1} + ',num2str(tauRe_linearfit(2))],...
            'Location','NorthWest');
        set(h_legend,'FontSize',14);
        title(['\epsilon_{raw}=',num2str(tauRe_nocal_relerr),...
            ', \epsilon_{cal}=',num2str(tauRe_cal_relerr)]);
        print([folder_Plots,'tauRe_Calibration_S1_',Instrument,'.png'],'-dpng');
        
        %plot heatflux
        figure(i+3*N_AllInstruments); clf;
        plot(heatflux_S1,heatflux_Instrument,'o'); hold on;
        plot([0 0.3],[0 0.3]); %plot 1-1 line
        plot([0 0.3],polyval(heatflux_linearfit,[0 0.3])); %plot linear fit
        xlabel('w`T`_{S1} (K m/s)','FontSize',14);
        ylabel(['w`T`_{',Instrument,'} (K m/s)'],'FontSize',14);
        set(gca,'FontSize',14);
        h_legend = legend('data','1-1 line',...
            ['fit: ',num2str(heatflux_linearfit(1)),'w`T`_{S1} + ',num2str(heatflux_linearfit(2))],...
            'Location','NorthWest');
        set(h_legend,'FontSize',14);
        title(['\epsilon_{raw}=',num2str(heatflux_nocal_relerr),...
            ', \epsilon_{cal}=',num2str(heatflux_cal_relerr)]);
        print([folder_Plots,'heatflux_Calibration_S1_',Instrument,'.png'],'-dpng');
    end
    
    %look at dust concentrations
    if strcmp(AllInstruments_Type(i),'Dust')
        
        %get indices for comparison to D1
        ind_comparison = intersect(find(~isnan(n1_matrix(:,i))),ind_onshorewind);
        
        %get dust concentrations for D1 and instrument
        n1_D1 = n1_matrix(ind_comparison,ind_D1);
        n1_Instrument = n1_matrix(ind_comparison,i);
        n2_D1 = n2_matrix(ind_comparison,ind_D1);
        n2_Instrument = n2_matrix(ind_comparison,i);
        n3_D1 = n3_matrix(ind_comparison,ind_D1);
        n3_Instrument = n3_matrix(ind_comparison,i);
        n4_D1 = n4_matrix(ind_comparison,ind_D1);
        n4_Instrument = n4_matrix(ind_comparison,i);
        n5_D1 = n5_matrix(ind_comparison,ind_D1);
        n5_Instrument = n5_matrix(ind_comparison,i);
        n6_D1 = n6_matrix(ind_comparison,ind_D1);
        n6_Instrument = n6_matrix(ind_comparison,i);
        n7_D1 = n7_matrix(ind_comparison,ind_D1);
        n7_Instrument = n7_matrix(ind_comparison,i);
        n8_D1 = n8_matrix(ind_comparison,ind_D1);
        n8_Instrument = n8_matrix(ind_comparison,i);

        %compute simple fit
        cal_n1_simple = mean(n1_Instrument./n1_D1);
        cal_n2_simple = mean(n2_Instrument./n2_D1);
        cal_n3_simple = mean(n3_Instrument./n3_D1);
        cal_n4_simple = mean(n4_Instrument./n4_D1);
        cal_n5_simple = mean(n5_Instrument./n5_D1);
        cal_n6_simple = mean(n6_Instrument./n6_D1);
        cal_n7_simple = mean(n7_Instrument./n7_D1);
        cal_n8_simple = mean(n8_Instrument./n8_D1);
        
        %compute linear fit
        cal_n1_linearfit = polyfit(n1_D1,n1_Instrument,1);
        cal_n2_linearfit = polyfit(n2_D1,n2_Instrument,1);
        cal_n3_linearfit = polyfit(n3_D1,n3_Instrument,1);
        cal_n4_linearfit = polyfit(n4_D1,n4_Instrument,1);
        cal_n5_linearfit = polyfit(n5_D1,n5_Instrument,1);
        cal_n6_linearfit = polyfit(n6_D1,n6_Instrument,1);
        cal_n7_linearfit = polyfit(n7_D1,n7_Instrument,1);
        cal_n8_linearfit = polyfit(n8_D1,n8_Instrument,1);

        %get concentrations for simply calibrated Instrument
        n1_simplecal_Instrument = n1_Instrument/cal_n1_simple;
        n2_simplecal_Instrument = n2_Instrument/cal_n2_simple;
        n3_simplecal_Instrument = n3_Instrument/cal_n3_simple;
        n4_simplecal_Instrument = n4_Instrument/cal_n4_simple;
        n5_simplecal_Instrument = n5_Instrument/cal_n5_simple;
        n6_simplecal_Instrument = n6_Instrument/cal_n6_simple;
        n7_simplecal_Instrument = n7_Instrument/cal_n7_simple;
        n8_simplecal_Instrument = n8_Instrument/cal_n8_simple;
        
        %get concentration values for linear calibrated Instrument
        n1_linearcal_Instrument = (n1_Instrument-cal_n1_linearfit(2))/cal_n1_linearfit(1);
        n2_linearcal_Instrument = (n2_Instrument-cal_n2_linearfit(2))/cal_n2_linearfit(1);
        n3_linearcal_Instrument = (n3_Instrument-cal_n3_linearfit(2))/cal_n3_linearfit(1);
        n4_linearcal_Instrument = (n4_Instrument-cal_n4_linearfit(2))/cal_n4_linearfit(1);
        n5_linearcal_Instrument = (n5_Instrument-cal_n5_linearfit(2))/cal_n5_linearfit(1);
        n6_linearcal_Instrument = (n6_Instrument-cal_n6_linearfit(2))/cal_n6_linearfit(1);
        n7_linearcal_Instrument = (n7_Instrument-cal_n7_linearfit(2))/cal_n7_linearfit(1);
        n8_linearcal_Instrument = (n8_Instrument-cal_n8_linearfit(2))/cal_n8_linearfit(1);
   
        %compute relative error for dust concentration
        n1_nocal_relerr = mean(abs((n1_Instrument-n1_D1)./(n1_D1)));
        n1_simplecal_relerr = mean(abs((n1_simplecal_Instrument-n1_D1)./(n1_D1)));
        n1_linearcal_relerr = mean(abs((n1_linearcal_Instrument-n1_D1)./(n1_D1)));
        
        n2_nocal_relerr = mean(abs((n2_Instrument-n2_D1)./(n2_D1)));
        n2_simplecal_relerr = mean(abs((n2_simplecal_Instrument-n2_D1)./(n2_D1)));
        n2_linearcal_relerr = mean(abs((n2_linearcal_Instrument-n2_D1)./(n2_D1)));
        
        n3_nocal_relerr = mean(abs((n3_Instrument-n3_D1)./(n3_D1)));
        n3_simplecal_relerr = mean(abs((n3_simplecal_Instrument-n3_D1)./(n3_D1)));
        n3_linearcal_relerr = mean(abs((n3_linearcal_Instrument-n3_D1)./(n3_D1)));
        
        n4_nocal_relerr = mean(abs((n4_Instrument-n4_D1)./(n4_D1)));
        n4_simplecal_relerr = mean(abs((n4_simplecal_Instrument-n4_D1)./(n4_D1)));
        n4_linearcal_relerr = mean(abs((n4_linearcal_Instrument-n4_D1)./(n4_D1)));
        
        n5_nocal_relerr = mean(abs((n5_Instrument-n5_D1)./(n5_D1)));
        n5_simplecal_relerr = mean(abs((n5_simplecal_Instrument-n5_D1)./(n5_D1)));
        n5_linearcal_relerr = mean(abs((n5_linearcal_Instrument-n5_D1)./(n5_D1)));
        
        n6_nocal_relerr = mean(abs((n6_Instrument-n6_D1)./(n6_D1)));
        n6_simplecal_relerr = mean(abs((n6_simplecal_Instrument-n6_D1)./(n6_D1)));
        n6_linearcal_relerr = mean(abs((n6_linearcal_Instrument-n6_D1)./(n6_D1)));
        
        n7_nocal_relerr = mean(abs((n7_Instrument-n7_D1)./(n7_D1)));
        n7_simplecal_relerr = mean(abs((n7_simplecal_Instrument-n7_D1)./(n7_D1)));
        n7_linearcal_relerr = mean(abs((n7_linearcal_Instrument-n7_D1)./(n7_D1)));

        n8_nocal_relerr = mean(abs((n8_Instrument-n8_D1)./(n8_D1)));
        n8_simplecal_relerr = mean(abs((n8_simplecal_Instrument-n8_D1)./(n8_D1)));
        n8_linearcal_relerr = mean(abs((n8_linearcal_Instrument-n8_D1)./(n8_D1)));
        
        %plot dust concentrations - n1
        figure(i+1*N_AllInstruments); clf;
        plot(n1_D1,n1_Instrument,'o'); hold on;
        plot([0 max(n1_D1)],cal_n1_simple*[0 max(n1_D1)]); %plot simple fit calibration
        plot([0 max(n1_D1)],polyval(cal_n1_linearfit,[0 max(n1_D1)])); %plot linear fit calibration
        xlabel('n1_{D1}','FontSize',14);
        ylabel(['n1_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n1_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n1_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n1_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n1_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n1_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n2
        figure(i+2*N_AllInstruments); clf;
        plot(n2_D1,n2_Instrument,'o'); hold on;
        plot([0 max(n2_D1)],cal_n2_simple*[0 max(n2_D1)]); %plot simple fit calibration
        plot([0 max(n2_D1)],polyval(cal_n2_linearfit,[0 max(n2_D1)])); %plot linear fit calibration
        xlabel('n2_{D1}','FontSize',14);
        ylabel(['n2_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n2_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n2_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n2_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n2_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n2_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n3
        figure(i+3*N_AllInstruments); clf;
        plot(n3_D1,n3_Instrument,'o'); hold on;
        plot([0 max(n3_D1)],cal_n3_simple*[0 max(n3_D1)]); %plot simple fit calibration
        plot([0 max(n3_D1)],polyval(cal_n3_linearfit,[0 max(n3_D1)])); %plot linear fit calibration
        xlabel('n3_{D1}','FontSize',14);
        ylabel(['n3_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n3_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n3_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n3_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n3_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n3_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n4
        figure(i+4*N_AllInstruments); clf;
        plot(n4_D1,n4_Instrument,'o'); hold on;
        plot([0 max(n4_D1)],cal_n4_simple*[0 max(n4_D1)]); %plot simple fit calibration
        plot([0 max(n4_D1)],polyval(cal_n4_linearfit,[0 max(n4_D1)])); %plot linear fit calibration
        xlabel('n4_{D1}','FontSize',14);
        ylabel(['n4_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n4_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n4_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n4_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n4_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n4_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n5
        figure(i+5*N_AllInstruments); clf;
        plot(n5_D1,n5_Instrument,'o'); hold on;
        plot([0 max(n5_D1)],cal_n5_simple*[0 max(n5_D1)]); %plot simple fit calibration
        plot([0 max(n5_D1)],polyval(cal_n5_linearfit,[0 max(n5_D1)])); %plot linear fit calibration
        xlabel('n5_{D1}','FontSize',14);
        ylabel(['n5_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n5_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n5_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n5_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n5_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n5_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n6
        figure(i+6*N_AllInstruments); clf;
        plot(n6_D1,n6_Instrument,'o'); hold on;
        plot([0 max(n6_D1)],cal_n6_simple*[0 max(n6_D1)]); %plot simple fit calibration
        plot([0 max(n6_D1)],polyval(cal_n6_linearfit,[0 max(n6_D1)])); %plot linear fit calibration
        xlabel('n6_{D1}','FontSize',14);
        ylabel(['n6_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n6_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n6_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n6_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n6_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n6_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n7
        figure(i+7*N_AllInstruments); clf;
        plot(n7_D1,n7_Instrument,'o'); hold on;
        plot([0 max(n7_D1)],cal_n7_simple*[0 max(n7_D1)]); %plot simple fit calibration
        plot([0 max(n7_D1)],polyval(cal_n7_linearfit,[0 max(n7_D1)])); %plot linear fit calibration
        xlabel('n7_{D1}','FontSize',14);
        ylabel(['n7_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n7_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n7_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n7_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n7_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n7_Calibration_D1_',Instrument,'.png'],'-dpng');
        
        %plot dust concentrations - n8
        figure(i+8*N_AllInstruments); clf;
        plot(n8_D1,n8_Instrument,'o'); hold on;
        plot([0 max(n8_D1)],cal_n8_simple*[0 max(n8_D1)]); %plot simple fit calibration
        plot([0 max(n8_D1)],polyval(cal_n8_linearfit,[0 max(n8_D1)])); %plot linear fit calibration
        xlabel('n8_{D1}','FontSize',14);
        ylabel(['n8_{,',Instrument,'} (Pa)'],'FontSize',14);
        title(strcat('Cal factor = ',num2str(cal_n8_simple)),'FontSize',14);
        h_legend = legend(...
            ['data, \epsilon = ',num2str(n8_nocal_relerr)],...
            ['simple fit, \epsilon =',num2str(n8_simplecal_relerr)],...
            ['linear fit, \epsilon =',num2str(n8_linearcal_relerr)],...
            'Location','NorthWest');
        set(gca,'FontSize',14);
        set(h_legend,'FontSize',14);
        print([folder_Plots,'n8_Calibration_D1_',Instrument,'.png'],'-dpng');
    end
end

%% Restore function path to default value
restoredefaultpath;