%% SCRIPT TO GENERATE TIMESERIES OF FLUX AND WIND DATA
% Dependencies: ExtractVariableTimeInterval

%%clear existing data and load processed data and metadata
clear all;
folder_ProcessedData = '../../../../../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/RawPlots/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to processed data
load(Metadata_Path); %load meta data
load(ProcessedData_Path); %load processed data

%set start and end time
StartTime = datetime(2015,6,3,17,20,00);
EndTime = datetime(2015,6,3,17,21,00);

[u_Sonic, t_Sonic, ~, ~] = ExtractVariableTimeInterval(ProcessedData.Sonic.S1,StartTime,EndTime,'u','int','int');
[n_Wenglor, t_Wenglor, ~, ~] = ExtractVariableTimeInterval(ProcessedData.Wenglor.W1,StartTime,EndTime,'n','int','int');
n_Wenglor = n_Wenglor*25;

%generate plot
figure(1); clf;
[AX,H1,H2] = plotyy(t_Sonic,u_Sonic,t_Wenglor,n_Wenglor);
set(H1,'LineWidth',1,'Color','b');
set(H2,'LineWidth',1,'Color','r');
ylabel(AX(1),'u_{0.5m} (m/s)','FontSize',12);
ylabel(AX(2),'n_{W1} (1/s)','FontSize',12);
title(['Oceano ',datestr(StartTime,'YYYY-mm-dd')]);
set(AX(1),'FontSize',12,'YColor','b');
set(AX(2),'FontSize',12,'YColor','r');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 5])
print([folder_Plots,'S1_W1_raw_',datestr(StartTime,'yyyy-mm-dd_HHMM'),'.png'],'-dpng');