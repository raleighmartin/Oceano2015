%% SCRIPT TO GENERATE PROFILES OF FLUX DATA AND COMPARE TO SHEAR VELOCITIES
% Dependencies: CreateTimeBlocks, ExtractVariableTimeInterval,
% reorient_anemometer_vanboxel2004

%%clear existing data and load processed data and metadata
% clear all;
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots
Metadata_Path = strcat(folder_ProcessedData,'Metadata_Oceano'); %get path to meta data
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_Oceano'); %get path to processed data
load(Metadata_Path); %load meta data
% load(ProcessedData_Path); %load processed data

%set physical parameters
rho_a = 1.23; %air density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,15,0); %10 minutes

%set anemometer for shear velocity calculations
ShearCalcSonic = 'S1';
ShearCalcType = 'Sonic';

%get start times and end times for sonic and total flux
SonicStartTimes = [ProcessedData.(ShearCalcType).(ShearCalcSonic).StartTime];
SonicEndTimes = [ProcessedData.(ShearCalcType).(ShearCalcSonic).EndTime];
FluxStartTimes = [ProcessedData.TotalFlux.StartTime];
FluxEndTimes = [ProcessedData.TotalFlux.EndTime];

%create time blocks based on time intervals
[BlockStartTimes, BlockEndTimes] = CreateTimeBlocks(SonicStartTimes, SonicEndTimes, ProfileTimeInterval);
N_Blocks = length(BlockStartTimes);

%initialize values for zbar and range
zbar_Block = zeros(N_Blocks,1)*NaN;
zbar_Block_max = zeros(N_Blocks,1)*NaN;
zbar_Block_min = zeros(N_Blocks,1)*NaN;

%initialize values for flux and shear velocity
Q_Block = zeros(N_Blocks,1)*NaN;
ust_Block = zeros(N_Blocks,1)*NaN;

%go through blocks to get total fluxes
for i = 1:N_Blocks
    i
    %get qz values
    [~, qz_t, qz_interval, qz_ind] = ExtractVariableTimeInterval(ProcessedData.TotalFlux,BlockStartTimes(i),BlockEndTimes(i),'qz','calc','calc');
    
    %check to make sure times fill whole interval window
    if (~isempty(qz_t)&&length(qz_interval)==1)&&((min(qz_t)==BlockStartTimes(i))&&(max(qz_t)==BlockEndTimes(i)))
        %get profile data for interval
        z = ProcessedData.TotalFlux(qz_interval).z.calc;
        qz = mean(ProcessedData.TotalFlux(qz_interval).qz.calc(cell2mat(qz_ind),:));
        
        %fit profile to get zbar
        [zbar,q0,Q] = qz_profilefit(z, qz);
        
        %calculate Reynolds stress
        u = ExtractVariableTimeInterval(ProcessedData.(ShearCalcType).(ShearCalcSonic),BlockStartTimes(i),BlockEndTimes(i),'u','int','int');
        v = ExtractVariableTimeInterval(ProcessedData.(ShearCalcType).(ShearCalcSonic),BlockStartTimes(i),BlockEndTimes(i),'v','int','int');
        w = ExtractVariableTimeInterval(ProcessedData.(ShearCalcType).(ShearCalcSonic),BlockStartTimes(i),BlockEndTimes(i),'w','int','int');
        [tauRe,ustRe] = computeReynoldsStress(u, v, w, rho_a);
        
        %add to lists
        ust_Block(i) = ustRe;
        
        %flux values only added if there were sufficient data points for profile calc
        if isstruct(zbar);
            zbar_Block(i) = zbar.mid;
            zbar_Block_max(i) = zbar.max;
            zbar_Block_min(i) = zbar.min;
            Q_Block(i) = Q;
        end
    end
end

%remove empty or inf entries
ind_good = find((~isnan(zbar_Block))&(abs(zbar_Block)~=Inf));
zbar_Block = zbar_Block(ind_good);
zbar_Block_min = zbar_Block_min(ind_good);
zbar_Block_max = zbar_Block_max(ind_good);
ust_Block = ust_Block(ind_good);
Q_Block = Q_Block(ind_good);

%plot results - saltation height versus u*
figure(1); clf;
errorbar(ust_Block,zbar_Block,zbar_Block_min,zbar_Block_max,'o');
xlabel('u_{*,Re} (m/s)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
ylim([0 0.15]);
set(gca,'FontSize',16);
print([folder_Plots,'FluxHeightUst_Wenglor.png'],'-dpng');

%plot results - saltation flux versus u*
figure(2); clf;
plot(ust_Block,Q_Block,'o');
xlabel('u_{*,Re} (m/s)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
corr_matrix = corrcoef(ust_Block,Q_Block);
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
ylim([0 35]);
set(gca,'FontSize',16);
print([folder_Plots,'FluxUst_Wenglor.png'],'-dpng');

%plot results - saltation flux versus tau
figure(3); clf;
plot(rho_a*ust_Block.^2,Q_Block,'o');
xlabel('\tau (Pa)','FontSize',16);
ylabel('Q (g/m/s)','FontSize',16);
corr_matrix = corrcoef(ust_Block.^2,Q_Block);
title(['r = ',num2str(corr_matrix(2))],'FontSize',16)
ylim([0 35]);
set(gca,'FontSize',16);
print([folder_Plots,'FluxTau_Wenglor.png'],'-dpng');