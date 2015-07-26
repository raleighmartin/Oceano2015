%% function to compute cross-correlation between timeseries 'X' and 'Y' (of
% equal length) with timestep dt (seconds) and maximum lag 'maxlag'
% (timesteps)
% Dependencies: NONE
% Used by: WindFluxLagAnalysis

%T_max -- how much does Y lag behind X?
%R_max -- what is correlation associated with this lag?
function [R, lags, T_max, R_max] = xcorrelation(X,Y,dt,maxlag)

%set maxlag if not specified by user
if nargin == 3
  maxlag = 1000;
end

%initialize output variables
R = zeros(2*maxlag+1,1);
lags = (-maxlag:maxlag)*dt;

%compute crosscorrelation for each lag
Xbar = mean(X);
Xstd = std(X);
Ybar = mean(Y);
Ystd = std(Y);
for i = -maxlag:maxlag
    if i<=0
        X0 = X((1-i):end)-Xbar;
        Y1 = Y(1:(end+i))-Ybar;
    else
        X0 = X(1:(end-i))-Xbar;
        Y1 = Y((1+i):end)-Ybar;
    end
    R(i+1+maxlag) = mean(X0.*Y1)/(Xstd*Ystd);
end

%find peak correlation
R_max_ind = find(R==max(R));
R_max = mean(R(R_max_ind));
T_max = mean(lags(R_max_ind));