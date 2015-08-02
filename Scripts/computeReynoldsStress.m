%% function to compute Reynolds stress tau and u* from wind velocities
%u,v,w - raw wind velocity components (m/s)
%rho_a - air density (kg/m^3)

function [tauRe,ustRe] = computeReynoldsStress(u, v, w, rho_a)

%rotate instrument
[u, ~, w] = reorient_anemometers_vanboxel2004(u, v, w);

%make computations
u_bar = mean(u); %mean wind
w_bar = mean(w); %mean wind
tauRe = real(-rho_a*mean((u-u_bar).*(w-w_bar))); %Reynolds stress
ustRe = sqrt(tauRe/rho_a); %associated shear velocity