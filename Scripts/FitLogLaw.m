%% function to fit log law to profile of horizontal wind velocities
% INPUTS:
% 'z_profile' is a list of heights (m)
% 'u_profile' is a list of horizontal wind velocities (m/s)
% 'kappa' is the von Karman parameter (assumed 0.39 if none given)
% 'zL_profile' is the profile of stability parameter (assumed 0 if none given)
% OUTPUTS:
% 'ust_raw' is the value of u* (m/s) without stability correction
% 'z0_raw' is the value of z0 (m) without stability correction
% 'ust_stabcorr' is the value of u* (m/s) with stability correction
% 'z0_stabcorr' is the value of z0 (m) with stability correction

% DEPENDENCIES:
% USED BY:

function [ust_raw, z0_raw, ust_stabcorr, z0_stabcorr] = ...
    FitLogLaw(z_profile,u_profile,zL_profile,kappa)



%set z/L = 0 if no argument given
if nargin<3
    zL_profile = zeros(size(u_profile));
end

%set kappa = 0.39 if no argument given
if nargin<4
    kappa = 0.39;
end

%fit to get parameters for log law
P = polyfit(log(z_profile),u_profile,1);
ust_raw = P(1)*kappa; %raw u* for log law
z0_raw = exp(-P(2)/P(1)); %raw z0 for log law

%fourth-order polynomial coefficients to determine psi for stability correction from z/L, from Kaimal and Finnigan (1994) Table 1.1
P_psi = [-0.2473 -1.2570 -2.3943 -2.4641 0.0312];

%compute psi profile, law of wall stability correction
psi_profile = polyval(P_psi,zL_profile);

%for stability correction, modify values of z based on Kaimal and Finnigan (1994) Eqn. 1.37
%u_profile_stabcorr = u_profile+psi_profile*(kappa/ust_raw); %use raw u* for stability correction
P = polyfit(log(z_profile)-psi_profile,u_profile,1);
ust_stabcorr = P(1)*kappa; %stability corrected u* for log law
z0_stabcorr = exp(-P(2)/P(1)); %stability corrected z0 for log law

%see also "AOS 224 - Lecture 15 - Similarity theory #2.pdf" slide 7 to make
%sense of this