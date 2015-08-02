%% SCRIPT TO LOOK AT BSNE PROFILE GRAIN SIZES

%%clear existing data and load processed data and metadata
%clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/GrainSize/Profiles'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
%load(Metadata_Path); %load metadata
%load(ProcessedData_Path); %load processed data

Locations = {'A4','A3','A2','A1','Wenglor','B1','B2_B3','B4'};
linespec = {':',':',':',':','-','--','--','--','--'};
N_Profiles = length(ProcessedData.BSNE);

for i = 1:N_Profiles
    
    %get data for profile
    BSNE_i = ProcessedData.BSNE(i);  
    d_10 = BSNE_i.d_Profile.d_10;
    d_50 = BSNE_i.d_Profile.d_50;
    d_90 = BSNE_i.d_Profile.d_90;
    z = BSNE_i.z;
    name = BSNE_i.name;
    Q = BSNE_i.Q;
    Q_time = mean([BSNE_i.StartTime; BSNE_i.EndTime]);
    
    %restrict to entries with values
    z = z(~isnan(d_50));
    name = name(~isnan(d_50));
    d_10 = d_10(~isnan(d_10));
    d_50 = d_50(~isnan(d_50));
    d_90 = d_90(~isnan(d_90));
    
    %get surface values from mean of all samples that day
    ind_Surface = find([ProcessedData.GrainSize.Surface(:).Date]==BSNE_i.Date);
    d_10_Surface = mean([ProcessedData.GrainSize.Surface(ind_Surface).d_10_mm]);
    d_50_Surface = mean([ProcessedData.GrainSize.Surface(ind_Surface).d_50_mm]);
    d_90_Surface = mean([ProcessedData.GrainSize.Surface(ind_Surface).d_90_mm]);
    
%     %add surface samples to profiles
%     z = [0; z];
%     d_10 = [d_10_Surface; d_10];
%     d_50 = [d_10_Surface; d_50];
%     d_90 = [d_10_Surface; d_90];
    
    %plot profile
    figure(1); clf; hold on;
    title([datestr(Q_time),', Q = ',num2str(Q),' g/m/s'],'FontSize',16);

    plot(d_10,z,'o','MarkerSize',10);
    plot(d_50,z,'x','MarkerSize',10);
    text(d_50,z,name)
    plot(d_90,z,'+','MarkerSize',10);
    plot(d_10_Surface,0,'ko','MarkerSize',20);
    plot(d_50_Surface,0,'kx','MarkerSize',20);
    plot(d_90_Surface,0,'k+','MarkerSize',20);
    xlabel('d (mm)','FontSize',16);
    ylabel('z (m)','FontSize',16);
    h_legend = legend({'d_{10}','d_{50}','d_{90}'},'Location','EastOutside');
    set(h_legend,'FontSize',16);
    print([folder_Plots,'/GrainSizeProfile_',datestr(Q_time,'yyyy-mm-dd_HHMM'),'.png'],'-dpng');
end