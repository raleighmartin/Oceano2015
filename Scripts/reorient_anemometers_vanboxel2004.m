%% function to perform anemometer reorientation described by Van Boxel et al., 2014
% (u_0, v_0, w_0) are raw vectors
% (u_3, v_3, w_3) are reoriented vectors
% Dependencies: NONE
% Used by: ProcessUltrasonics, GetVelocityProfiles

function [u_3, v_3, w_3] = reorient_anemometers_vanboxel2004(u_0, v_0, w_0)

%% yaw rotation
if mean(u_0)>0
    theta = atan(mean(v_0)/mean(u_0));
else
    theta = atan(mean(v_0)/mean(u_0))+pi;
end
u_1 = u_0*cos(theta)+v_0*sin(theta);
v_1 = -u_0*sin(theta)+v_0*cos(theta);
w_1 = w_0;

%% pitch rotation
phi = atan(mean(w_1)/mean(u_1));
u_2 = u_1*cos(phi)+w_1*sin(phi);
v_2 = v_1;
w_2 = -u_1*sin(phi)+w_1*cos(phi);

%% roll rotation
psi = 0.5*atan(2*mean(v_2.*w_2)/(mean(v_2.^2)-mean(w_2.^2)));
u_3 = u_2;
v_3 = v_2*cos(psi)+w_2*sin(psi);
w_3 = -v_2*sin(psi)+w_2*cos(psi);

end