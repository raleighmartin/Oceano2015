%% SCRIPT TO LOOK AT CHANGES IN SALTATION HEIGHT WITH SURFACE GRAIN SIZE

clear all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationProfile/'; %folder for plots

%information about sites for analysis
Sites = {'Jericoacoara','RanchoGuadalupe','Oceano'};
Markers = {'x','o','v'};
N_Sites = length(Sites);

for i = 1:N_Sites
    %load metadata and processed data
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    Data{i} = load(ProcessedData_Path); %load processed data
end

rho_a = 1.23; %air density kg/m^3

for i = 1:N_Sites
    
    %initialize lists of daily values for plotting
    d50_Daily = [];
    zbar_Daily = [];
    
    %initialize list of interval values for plotting
    ust_Intervals = []; 
    zbar_Intervals = [];
    znorm_Intervals = [];
    d50_Intervals = [];
    
    %get list of dates
    Dates = unique([Data{i}.ProcessedData.GrainSize.Surface.Date]);
    N_Dates = length(Dates);
    
    for j = 1:N_Dates
        ind_Date = find([Data{i}.ProcessedData.GrainSize.Surface.Date]==Dates(j));
        d50_Date = median([Data{i}.ProcessedData.GrainSize.Surface(ind_Date).d_50_mm]);
        
        ind_Date = find([Data{i}.ProcessedData.BSNE.Date]==Dates(j));
        zbar_Interval = [Data{i}.ProcessedData.BSNE(ind_Date).zbar];
        
        if ~isempty(zbar_Interval)
            zbar_Interval = [zbar_Interval.mid];
            
            zbar_Date = median(zbar_Interval);
            zbar_Daily = [zbar_Daily, zbar_Date];
            d50_Daily = [d50_Daily, d50_Date];
            
            zbar_Intervals = [zbar_Intervals, zbar_Interval];
            znorm_Intervals = [znorm_Intervals, 1000*zbar_Interval/d50_Date];
            d50_Intervals = [d50_Intervals, d50_Date*ones(size(zbar_Interval))];
        end
        
        %compute shear velocity for each interval
        N_Intervals = length(zbar_Interval);
        
        if strcmp(Sites{i},'Oceano')
            WindData = Data{i}.ProcessedData.Sonic.S1;
        elseif strcmp(Sites{i},'RanchoGuadalupe')
            WindData = Data{i}.ProcessedData.Ultrasonic.U1;
        elseif strcmp(Sites{i},'Jericoacoara')
            WindData = Data{i}.ProcessedData.Ultrasonic.U1;
        end
        
        for k = ind_Date
            StartTime = Data{i}.ProcessedData.BSNE(k).StartTime;
            EndTime = Data{i}.ProcessedData.BSNE(k).EndTime;
            u = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'u','int','int');
            v = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'v','int','int');
            w = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'w','int','int');
            [~,ustRe] = computeReynoldsStress(u, v, w, rho_a);
            ust_Intervals = [ust_Intervals, ustRe];
        end
        
    end
    
    %add to cell arrays of values
    zbar_Daily_all{i} = zbar_Daily;
    d50_Daily_all{i} = d50_Daily;
    ust_Intervals_all{i} = ust_Intervals;
    zbar_Intervals_all{i} = zbar_Intervals;
    znorm_Intervals_all{i} = znorm_Intervals;
    d50_Intervals_all{i} = d50_Intervals;
end

%fit all zbar versus d50
P = polyfit([d50_Daily_all{:}],[zbar_Daily_all{:}],1);

%plot z_salt versus d_50
figure(1); clf; hold on;
for i = 1:N_Sites
    plot(d50_Daily_all{i},zbar_Daily_all{i},Markers{i},'MarkerSize',10);
end
plot([0 0.6],polyval(P,[0 0.6]));
h_legend = legend([Sites,'fit']);
set(h_legend,'FontSize',16,'Location','NorthWest');
xlabel('d_{50} (mm)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
set(gca,'FontSize',16);
ylim([0 0.12]);
print([folder_Plots,'d50_zSalt.png'],'-dpng');

%plot z_salt/d_50 versus u*
figure(2); clf; hold on;
for i = 1:N_Sites
    plot(ust_Intervals_all{i},znorm_Intervals_all{i},Markers{i},'MarkerSize',10);
end
h_legend = legend(Sites);
set(h_legend,'FontSize',16,'Location','NorthWest');
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{salt}/d_{50}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ust_zNorm.png'],'-dpng');

%plot z_salt/z_salt_pred versus u*
figure(3); clf; hold on;
for i = 1:N_Sites
    zbar_obs = zbar_Intervals_all{i};
    zbar_fit = polyval(P,d50_Intervals_all{i});
    plot(ust_Intervals_all{i},zbar_obs./zbar_fit,Markers{i},'MarkerSize',10);
end
h_legend = legend(Sites);
set(h_legend,'FontSize',16,'Location','NorthOutside');
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{salt,obs}/z_{salt,fit}','FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'ust_zSaltOverFit.png'],'-dpng');