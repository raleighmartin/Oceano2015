%% PLOT ALL VARIABLES WITH FIVE MINUTE MOVING AVERAGE

%% clear existing data and load processed data and metadata
%clearvars;

%% indicate where to save plots and resulting processed data
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output

%% set time interval for window average
T_wa = duration(0,1,0);

%% information about sites for analysis
Site_Names = {'Jericoacoara','RanchoGuadalupe','Oceano'};
Folder_Names = {'Jericoacoara2014','RanchoGuadalupe2015','Oceano2015'};
N_Sites = length(Site_Names);

%% load processed data and metadata for each site
%Data = cell(N_Sites,1); %initialize cell array containing data
%Metadata = cell(N_Sites,1); %initialize cell array containing metadata
SaveData_Path = cell(N_Sites,1); %initialize list of paths for saving data
folder_Plots = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Site_Names{i});
    Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Site_Names{i});
    %Data{i} = load(ProcessedData_Path); %load processed data
    %Metadata{i} = load(Metadata_Path); %load metadata
    SaveData_Path{i} = ProcessedData_Path; %list of paths for saving data
    folder_Plots{i} = ['../../',Folder_Names{i},'/PlotOutput/DiagnosticPlots/'];
end

%go through each site to produce plots
%for i = 1:N_Sites
for i = 3;
    InstrumentTypes = unique(Metadata{i}.InstrumentMetadata.InstrumentType);
    N_InstrumentTypes = length(InstrumentTypes);
    
    %go through all instrument types
    for j = N_InstrumentTypes
    %for j = 1:N_InstrumentTypes
        Instruments = fieldnames(Data{i}.ProcessedData.(InstrumentTypes{j}));
        N_Instruments = length(Instruments);

        %add additional variables for certain instruments
        AdditionalVariables = cell(0,0);
        AdditionalSubvariables = cell(0,0);
        AdditionalUnits = cell(0,0);
        if strcmp(InstrumentTypes{j},'Wenglor')
            AdditionalVariables{1} = 'flux';
            AdditionalSubvariables{1} = 'qz';
            AdditionalUnits{1} = 'g/m^2/s';
            AdditionalVariables{2} = 'flux';
            AdditionalSubvariables{2} = 'qzPerCount';
            AdditionalUnits{2} = '(g/m^2/s) / count';
        end
        N_AdditionalVariables = length(AdditionalVariables);
        
        %go through each instrument of this type
        %for k = 1:N_Instruments
        for k = 1;
            Variables = Metadata{i}.InstrumentVariables.VarNameGeneral(...
                strcmp(Metadata{i}.InstrumentVariables.Instrument,Instruments{k}));
            N_Variables = length(Variables);
            N_Intervals = length(Data{i}.ProcessedData.(InstrumentTypes{j}).(Instruments{k}));

            %get relative duration of each interval for plotting
            Interval_Durations = zeros(N_Intervals,1);
            for l = 1:N_Intervals
                Interval_Durations(l) = hours(Data{i}.ProcessedData.(InstrumentTypes{j}).(Instruments{k})(l).EndTime-...
                    Data{i}.ProcessedData.(InstrumentTypes{j}).(Instruments{k})(l).StartTime);
            end
            Total_Duration = sum(Interval_Durations);
            Interval_RelativeDurations = round(Interval_Durations/min(Interval_Durations));
            
            %convert relative durations to subplot indices for interval
            Interval_SubplotIndices = cell(N_Intervals,1);
            StartIndex_SubplotIndices = cumsum(Interval_RelativeDurations)-Interval_RelativeDurations;
            for l = 1:N_Intervals
                Interval_SubplotIndices{l} = (1:Interval_RelativeDurations(l))+StartIndex_SubplotIndices(l);
            end
            N_Subplot = sum(Interval_RelativeDurations);
            
            %inialize plots for variables
            close all;
            for m = 1:N_Variables
                figure(m); clf;
                %subplot(1,N_Subplot,Interval_SubplotIndices{1});
                hold on;
                Units = Data{i}.ProcessedData.(InstrumentTypes{j}).(Instruments{k})(1).(Variables{m}).Units; %get units
                ylabel([Site_Names{i},', ',Instruments{k},', ',Variables{m},' (',Units,')'],'FontSize',16);
            end
            
            for m = 1:N_AdditionalVariables
                figure(m+N_Variables); clf;
                %subplot(1,N_Subplot,Interval_SubplotIndices{1});
                hold on;
                Units = AdditionalUnits{m}; %get units
                ylabel([Site_Names{i},', ',Instruments{k},', ',AdditionalVariables{m},' (',Units,')'],'FontSize',16);
            end
            
            %go through each interval for this instrument, generate plots
            %for l = 1:N_Intervals
            for l = 7
                %get data for instrument                
                InstrumentDataInterval = Data{i}.ProcessedData.(InstrumentTypes{j}).(Instruments{k})(l);
                t = InstrumentDataInterval.t.int;
                
                %construct a timeseries of error pulses (0 = no error, 1 = error)
                [~, t_errorind] = intersect(t,InstrumentDataInterval.t.err); %indices of timesteps with error
                errorpulse = zeros(size(t));
                errorpulse(t_errorind) = 1;
                
                %window-averaged error pulse series and corresponding window-averaged times
                [errorpulse_wa, t_wa] = window_average(errorpulse, t, T_wa); 
                errorpulse_wa = ceil(errorpulse_wa); %convert window-averaged values to binary
                
                %go through and plot each variable
                for m = 1:N_Variables
                    VarData = InstrumentDataInterval.(Variables{m}).int; %get raw variable data
                    [VarData_wa, ~] = window_average(VarData, t, T_wa); %get window-averaged data              
                    
                    %plot data
                    figure(m);
                    %subplot(1,N_Subplot,Interval_SubplotIndices{l}); %subplot window size depends on data duration
                    plot(t_wa,VarData_wa); hold on;
                    
                    %plot error data
                    ErrorPlot_t = t_wa(errorpulse_wa==1);
                    ErrorPlot_y = VarData_wa(errorpulse_wa==1).*errorpulse_wa(errorpulse_wa==1);
                    plot(ErrorPlot_t,ErrorPlot_y,'.r','LineWidth',2);
                    
                    %set xlimits and label subplot
                    xlim([datenum(min(t)) datenum(max(t))])
                    title(datestr(InstrumentDataInterval.Date));
                    set(gca,'FontSize',16);
                end
                
                %go through and plot additional variables
                for m = 1:N_AdditionalVariables
                    VarData = InstrumentDataInterval.(AdditionalVariables{m}).(AdditionalSubvariables{m}); %get raw variable data
                    [VarData_wa, ~] = window_average(VarData, t, T_wa); %get window-averaged data              
                    
                    %plot data
                    figure(m+N_Variables);
                    %subplot(1,N_Subplot,Interval_SubplotIndices{l}); %subplot window size depends on data duration
                    plot(t_wa,VarData_wa); hold on;
                    
                    %plot error data
                    ErrorPlot_t = t_wa(errorpulse_wa==1);
                    ErrorPlot_y = VarData_wa(errorpulse_wa==1).*errorpulse_wa(errorpulse_wa==1);
                    plot(ErrorPlot_t,ErrorPlot_y,'.r','LineWidth',2);
                    
                    %set xlimits and label subplot
                    xlim([datenum(min(t)) datenum(max(t))])
                    title(datestr(InstrumentDataInterval.Date));
                    set(gca,'FontSize',16);
                end
            end
            
            %print all plots - basic variables
            for m = 1:N_Variables
                figure(m);
                %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*Total_Duration 5])
                print([folder_Plots{i},Site_Names{i},'_',Instruments{k},'_',Variables{m},'.png'],'-dpng');
            end
            
            %print all plots - additional variables
            for m = 1:N_AdditionalVariables
                figure(N_Variables+m);
                %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*Total_Duration 5])
                print([folder_Plots{i},Site_Names{i},'_',Instruments{k},'_',AdditionalSubvariables{m},'.png'],'-dpng');
            end
        end
    end
end