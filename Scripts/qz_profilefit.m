%% function to compute saltation flux profile assuming exponential
%z - observation heights (m)
%qz - vertical fluxes (g/m^2/s)
%Error calculation based on Wilkinson (1984)
function [zbar,q0,Q] = qz_profilefit(z, qz, sigmaz, sigmaqz)
    
%height-integrated flux from exponential fit (using only nonzero values of qz)
z_fit = z(qz~=0); %list of z values to fit
logqz_fit = log(qz(qz~=0)); %list of ln(qz) values to fit 

if length(z_fit)>2
    
    %linear fit to these values
    fit_params = polyfit(z_fit,logqz_fit,1);
    m_fit = fit_params(1); %slope of fit
    c_fit = fit_params(2); %intercept of fit

    %correlation coefficient for fit
    corr_params = corrcoef(z_fit,logqz_fit);
    r_fit = corr_params(2); %correlation for fit
    df_fit = length(z_fit); %number of degrees of freedom for fit

    %compute range of slopes based on t-statistic
    t_fit = t_95(df_fit); %t-value for 95% confidence interval
    %dm_fit = (t_fit/sqrt(df_fit-2))*(std(logqz_fit)/std(z_fit))*sqrt(1-r_fit^2); %original method, range of slope values from t-statistic
    dm_fit = (t_fit/sqrt(df_fit-2))*(log(mean(sigmaz))/mean(sigmaz))*sqrt(1-r_fit^2); %alternative method, range of slope values from t-statistic
    m_min = (m_fit-dm_fit); %minimum fit slope
    m_max = (m_fit+dm_fit); %maximum fit slope

    %from this, get range of flux e-folding heights
    zbar_mid = -1/m_fit; %e-folding saltation height (m)
    zbar_min = -1/m_min; %minimum e-folding height from error range
    zbar_max = -1/m_max; %maximum e-folding height from error range

    %also compute total flux from fit
    q0 = exp(c_fit); %basal flux (g/m^2/s)
    Q = q0*zbar_mid; %g/m/s

    %put all saltation height values together into structured array
    zbar = struct('mid',zbar_mid,'min',zbar_min,'max',zbar_max);
else
    q0 = NaN;
    Q = NaN;
    zbar = NaN;
end