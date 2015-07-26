%% function, given BSNE mass, height of opening, width of opening, and duration of observation
% to determine the mass flux
% Dependencies: NONE
% Used by: ProcessBSNEs

function [qz_g_m2_s] = qzCalc(MassBSNE_g, HeightBSNE_cm, WidthBSNE_cm, DurationBSNE)

qz_g_m2_s = (MassBSNE_g)/(1e-4*HeightBSNE_cm*WidthBSNE_cm*seconds(DurationBSNE));

end