%% SCRIPT TO LOOK AT CHANGES IN SURFACE GRAIN SIZE THROUGH TIME

%%clear existing data and load processed data and metadata
%clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/GrainSize/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to saving meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to saving processed data
%load(Metadata_Path); %load metadata
%load(ProcessedData_Path); %load processed data

Locations = {'A4','A3','A2','A1','Wenglor','B1','B2_B3','B4'};
linespec = {':',':',':',':','-','--','--','--','--'};
N_Locations = length(Locations);

%initialize figure
figure(1); clf; hold on;
xlabel('time','FontSize',16);
ylabel('d_{10} (mm)','FontSize',16);

figure(2); clf; hold on;
xlabel('time','FontSize',16);
ylabel('d_{50} (mm)','FontSize',16);

figure(3); clf; hold on;
xlabel('time','FontSize',16);
ylabel('d_{90} (mm)','FontSize',16);

for i = 1:N_Locations
    Locations{i}
    ind_Location = find(strcmp({ProcessedData.GrainSize.Surface(:).Location},Locations{i}));
    
    figure(1);
    plot([ProcessedData.GrainSize.Surface(ind_Location).CollectionTime],...
        [ProcessedData.GrainSize.Surface(ind_Location).d_10_mm],...
        linespec{i});
    mean([ProcessedData.GrainSize.Surface(ind_Location).d_10_mm])
    
    figure(2);
    plot([ProcessedData.GrainSize.Surface(ind_Location).CollectionTime],...
        [ProcessedData.GrainSize.Surface(ind_Location).d_50_mm],...
        linespec{i});
    mean([ProcessedData.GrainSize.Surface(ind_Location).d_50_mm])

    
    figure(3);
    plot([ProcessedData.GrainSize.Surface(ind_Location).CollectionTime],...
        [ProcessedData.GrainSize.Surface(ind_Location).d_90_mm],...
        linespec{i});
    mean([ProcessedData.GrainSize.Surface(ind_Location).d_90_mm])

end

figure(1);
legend(Locations,'Location','NorthWest');
print([folder_Plots,'Surface_d10.png'],'-dpng');

figure(2);
legend(Locations,'Location','NorthWest');
print([folder_Plots,'Surface_d50.png'],'-dpng');

figure(3);
legend(Locations,'Location','NorthWest');
print([folder_Plots,'Surface_d90.png'],'-dpng');