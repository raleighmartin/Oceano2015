%% function to compute saltation flux profile assuming exponential
%z - observation heights (m)
%qz - vertical fluxes (g/m^2/s)
function [zbar,q0,Q] = qz_profilefit(z, qz)
    
%height-integrated flux from exponential fit (using only nonzero values of qz)
z_fit = z(qz~=0); %list of z values to fit
logqz_fit = log(qz(qz~=0)); %list of ln(qz) values to fit 

if length(z_fit)>2
    %linear fit to these values
    fit_params = polyfit(z_fit,logqz_fit,1);
    m_fit = fit_params(1); %slope of fit
    c_fit = fit_params(2); %intercept of fit

    %from this, get range of flux e-folding heights
    zbar = -1/m_fit; %e-folding saltation height (m)

    %also compute total flux from fit
    q0 = exp(c_fit); %basal flux (g/m^2/s)
    Q = q0*zbar; %g/m/s
else
    q0 = NaN;
    Q = NaN;
    zbar = NaN;
end