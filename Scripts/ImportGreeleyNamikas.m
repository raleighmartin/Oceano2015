%% clear all
clear all;

%% filepath information
folder_LiteratureData = '../../../Google Drive/Data/Literature/'; %location where Namikas and Greeley data are stored
filename_Greeley96 = [folder_LiteratureData,'Greeley96.dat'];
filename_Namikas03 = [folder_LiteratureData,'Namikas03verMF.dat'];
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots
SaveData_Path = strcat(folder_ProcessedData,'GreeleyNamikasData');

%% manually enter basic info
N_Greeley96 = 10;
N_Namikas03 = 9;
ust_mpers_Greeley96 = [0.47; 0.51; 0.49; 0.31; 0.35; 0.31; 0.54; 0.49; 0.41; 0.42]; %values used in Kok et al. (2008)
ust_mpers_Namikas03 = [0.32; 0.27; 0.32; 0.30; 0.38; 0.37; 0.38; 0.47; 0.63];
d50_mm_Greeley96 = 0.230;
d50_mm_Namikas03 = 0.250;

%% Initialize processing variables
delimiter = '\t';
startRow = 2;

%% Format string for each line of text:
formatSpec_Greeley96 = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
formatSpec_Namikas03 = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text files
fileID_Greeley96 = fopen(filename_Greeley96,'r');
fileID_Namikas03 = fopen(filename_Namikas03,'r');

%% Read columns of data according to format strings
dataArray_Greeley96 = textscan(fileID_Greeley96, formatSpec_Greeley96, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
dataArray_Namikas03 = textscan(fileID_Namikas03, formatSpec_Namikas03, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text files
fclose(fileID_Greeley96);
fclose(fileID_Namikas03);

%% Process Greeley data

% Initialize cell arrays with all values
z_m_Greeley96 = cell(N_Greeley96,1);
q_gperm2s_Greeley96 = cell(N_Greeley96,1);
RunName_Greeley96 = cell(N_Greeley96,1);

% Allocate imported array to column variable names and cell arrays
z_cm_1a = dataArray_Greeley96{:, 1};
q_gpercm2s_1a = dataArray_Greeley96{:, 2};
z_m_Greeley96{1} = 1e-2*z_cm_1a(~isnan(z_cm_1a));
q_gperm2s_Greeley96{1} = 1e4*q_gpercm2s_1a(~isnan(q_gpercm2s_1a));
RunName_Greeley96{1} = '1a';

z_cm_1b = dataArray_Greeley96{:, 3};
q_gpercm2s_1b = dataArray_Greeley96{:, 4};
z_m_Greeley96{2} = 1e-2*z_cm_1b(~isnan(z_cm_1b));
q_gperm2s_Greeley96{2} = 1e4*q_gpercm2s_1b(~isnan(q_gpercm2s_1b));
RunName_Greeley96{2} = '1b';

z_cm_1c = dataArray_Greeley96{:, 5};
q_gpercm2s_1c = dataArray_Greeley96{:, 6};
z_m_Greeley96{3} = 1e-2*z_cm_1c(~isnan(z_cm_1c));
q_gperm2s_Greeley96{3} = 1e4*q_gpercm2s_1c(~isnan(q_gpercm2s_1c));
RunName_Greeley96{3} = '1c';

z_cm_2 = dataArray_Greeley96{:, 7};
q_gpercm2s_2 = dataArray_Greeley96{:, 8};
z_m_Greeley96{4} = 1e-2*z_cm_2(~isnan(z_cm_2));
q_gperm2s_Greeley96{4} = 1e4*q_gpercm2s_2(~isnan(q_gpercm2s_2));
RunName_Greeley96{4} = '2';

z_cm_3 = dataArray_Greeley96{:, 9};
q_gpercm2s_3 = dataArray_Greeley96{:, 10};
z_m_Greeley96{5} = 1e-2*z_cm_3(~isnan(z_cm_3));
q_gperm2s_Greeley96{5} = 1e4*q_gpercm2s_3(~isnan(q_gpercm2s_3));
RunName_Greeley96{5} = '3';

z_cm_4 = dataArray_Greeley96{:, 11};
q_gpercm2s_4 = dataArray_Greeley96{:, 12};
z_m_Greeley96{6} = 1e-2*z_cm_4(~isnan(z_cm_4));
q_gperm2s_Greeley96{6} = 1e4*q_gpercm2s_4(~isnan(q_gpercm2s_4));
RunName_Greeley96{6} = '4';

z_cm_5a = dataArray_Greeley96{:, 13};
q_gpercm2s_5a = dataArray_Greeley96{:, 14};
z_m_Greeley96{7} = 1e-2*z_cm_5a(~isnan(z_cm_5a));
q_gperm2s_Greeley96{7} = 1e4*q_gpercm2s_5a(~isnan(q_gpercm2s_5a));
RunName_Greeley96{7} = '5a';

z_cm_5b = dataArray_Greeley96{:, 15};
q_gpercm2s_5b = dataArray_Greeley96{:, 16};
z_m_Greeley96{8} = 1e-2*z_cm_5b(~isnan(z_cm_5b));
q_gperm2s_Greeley96{8} = 1e4*q_gpercm2s_5b(~isnan(q_gpercm2s_5b));
RunName_Greeley96{8} = '5b';

z_cm_6a = dataArray_Greeley96{:, 17};
q_gpercm2s_6a = dataArray_Greeley96{:, 18};
z_m_Greeley96{9} = 1e-2*z_cm_6a(~isnan(z_cm_6a));
q_gperm2s_Greeley96{9} = 1e4*q_gpercm2s_6a(~isnan(q_gpercm2s_6a));
RunName_Greeley96{9} = '6a';

z_cm_6b = dataArray_Greeley96{:, 19};
q_gpercm2s_6b = dataArray_Greeley96{:, 20};
z_m_Greeley96{10} = 1e-2*z_cm_6b(~isnan(z_cm_6b));
q_gperm2s_Greeley96{10} = 1e4*q_gpercm2s_6b(~isnan(q_gpercm2s_6b));
RunName_Greeley96{10} = '6b';

% Calculate saltation flux and height
zbar_m_Greeley96 = zeros(N_Greeley96,1);
Q_gperm2s_Greeley96 = zeros(N_Greeley96,1);

for i=1:N_Greeley96
    [zbar,q0,Q] = qz_profilefit(z_m_Greeley96{i}, q_gperm2s_Greeley96{i});
    zbar_m_Greeley96(i) = zbar;
    Q_gperm2s_Greeley96(i) = Q;
    
    figure(i); clf; hold on;
    plot(q_gperm2s_Greeley96{i},z_m_Greeley96{i},'x')
    z_range = linspace(0, max(z_m_Greeley96{i}),50);
    q_pred = q0*exp(-z_range./zbar);
    plot(q_pred,z_range);
    set(gca,'xscale','log');
    set(gca,'FontSize',16);
    xlabel('q (g m^{-2} s^{-1})','FontSize',16);
    ylabel('z (m)','FontSize',16);
    title(['Run = ',RunName_Greeley96{i},', u_{*} = ',num2str(ust_mpers_Greeley96(i)),' m/s, Q = ',num2str(Q),' g/m/s, z_{q} = ',num2str(zbar),' m']);
    print([folder_Plots,'Greeley96_qzprofile_',RunName_Greeley96{i},'.png'],'-dpng');
end

%% Process Namikas data

% Initialize cell arrays with all values
z_m_Namikas03 = cell(N_Namikas03,1);
q_gperm2s_Namikas03 = cell(N_Namikas03,1);
RunName_Namikas03 = cell(N_Namikas03,1);

% Allocate imported array to column variable names and cell arrays
z_cm_5 = dataArray_Namikas03{:, 1};
q_kgperm2s_5 = dataArray_Namikas03{:, 2};
z_m_Namikas03{1} = z_cm_5(~isnan(q_kgperm2s_5));
q_gperm2s_Namikas03{1} = 1e3*q_kgperm2s_5(~isnan(q_kgperm2s_5));
RunName_Namikas03{1} = '5';

z_cm_3 = dataArray_Namikas03{:, 1};
q_kgperm2s_3 = dataArray_Namikas03{:, 3};
z_m_Namikas03{2} = z_cm_3(~isnan(q_kgperm2s_3));
q_gperm2s_Namikas03{2} = 1e3*q_kgperm2s_3(~isnan(q_kgperm2s_3));
RunName_Namikas03{2} = '3';

z_cm_4 = dataArray_Namikas03{:, 1};
q_kgperm2s_4 = dataArray_Namikas03{:, 4};
z_m_Namikas03{3} = z_cm_4(~isnan(q_kgperm2s_4));
q_gperm2s_Namikas03{3} = 1e3*q_kgperm2s_4(~isnan(q_kgperm2s_4));
RunName_Namikas03{3} = '4';

z_cm_8 = dataArray_Namikas03{:, 1};
q_kgperm2s_8 = dataArray_Namikas03{:, 5};
z_m_Namikas03{4} = z_cm_8(~isnan(q_kgperm2s_8));
q_gperm2s_Namikas03{4} = 1e3*q_kgperm2s_8(~isnan(q_kgperm2s_8));
RunName_Namikas03{4} = '8';

z_cm_10 = dataArray_Namikas03{:, 1};
q_kgperm2s_10 = dataArray_Namikas03{:, 6};
z_m_Namikas03{5} = z_cm_10(~isnan(q_kgperm2s_10));
q_gperm2s_Namikas03{5} = 1e3*q_kgperm2s_10(~isnan(q_kgperm2s_10));
RunName_Namikas03{5} = '10';

z_cm_6 = dataArray_Namikas03{:, 1};
q_kgperm2s_6 = dataArray_Namikas03{:, 7};
z_m_Namikas03{6} = z_cm_6(~isnan(q_kgperm2s_6));
q_gperm2s_Namikas03{6} = 1e3*q_kgperm2s_6(~isnan(q_kgperm2s_6));
RunName_Namikas03{6} = '6';

z_cm_9 = dataArray_Namikas03{:, 1};
q_kgperm2s_9 = dataArray_Namikas03{:, 8};
z_m_Namikas03{7} = z_cm_9(~isnan(q_kgperm2s_9));
q_gperm2s_Namikas03{7} = 1e3*q_kgperm2s_9(~isnan(q_kgperm2s_9));
RunName_Namikas03{7} = '9';

z_cm_13 = dataArray_Namikas03{:, 1};
q_kgperm2s_13 = dataArray_Namikas03{:, 9};
z_m_Namikas03{8} = z_cm_13(~isnan(q_kgperm2s_13));
q_gperm2s_Namikas03{8} = 1e3*q_kgperm2s_13(~isnan(q_kgperm2s_13));
RunName_Namikas03{8} = '13';

z_cm_14 = dataArray_Namikas03{:, 1};
q_kgperm2s_14 = dataArray_Namikas03{:, 10};
z_m_Namikas03{9} = z_cm_14(~isnan(q_kgperm2s_14));
q_gperm2s_Namikas03{9} = 1e3*q_kgperm2s_14(~isnan(q_kgperm2s_14));
RunName_Namikas03{9} = '14';

% Calculate saltation flux and height
zbar_m_Namikas03 = zeros(N_Namikas03,1);
Q_gperm2s_Namikas03 = zeros(N_Namikas03,1);

for i=1:N_Namikas03
    [zbar,q0,Q] = qz_profilefit(z_m_Namikas03{i}, q_gperm2s_Namikas03{i});
    zbar_m_Namikas03(i) = zbar;
    Q_gperm2s_Namikas03(i) = Q;
    
    figure(i+N_Namikas03); clf; hold on;
    plot(q_gperm2s_Namikas03{i},z_m_Namikas03{i},'x')
    z_range = linspace(0, max(z_m_Namikas03{i}),50);
    q_pred = q0*exp(-z_range./zbar);
    plot(q_pred,z_range);
    set(gca,'xscale','log');
    set(gca,'FontSize',16);
    xlabel('q (g m^{-2} s^{-1})','FontSize',16);
    ylabel('z (m)','FontSize',16);
    title(['Run ',RunName_Namikas03{i},', u_{*} = ',num2str(ust_mpers_Namikas03(i)),' m/s, Q = ',num2str(Q),' g/m/s, z_{q} = ',num2str(zbar),' m']);
    print([folder_Plots,'Namikas03_qzprofile_',RunName_Namikas03{i},'.png'],'-dpng');
end

%plot zbar(e) versus u*
figure(20); clf;
plot(ust_mpers_Greeley96,zbar_m_Greeley96,'r^',...
    ust_mpers_Namikas03,zbar_m_Namikas03,'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{q,e} (m)','FontSize',16);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_zsalt_e.png'],'-dpng');

%plot z50 versus u*
figure(21); clf;
plot(ust_mpers_Greeley96,zbar_m_Greeley96*log(2),'r^',...
    ust_mpers_Namikas03,zbar_m_Namikas03*log(2),'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{q,50} (m)','FontSize',16);
ylim([0.03 0.04]);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_zsalt_50.png'],'-dpng');

%plot Q versus u*
figure(22); clf;
plot(ust_mpers_Greeley96,Q_gperm2s_Greeley96,'r^',...
    ust_mpers_Namikas03,Q_gperm2s_Namikas03,'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('Q (g/m^2/s)','FontSize',16);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_Q.png'],'-dpng');

%% Clear temporary variables
clearvars filename* fileID* dataArray* formatSpec*;

%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03');