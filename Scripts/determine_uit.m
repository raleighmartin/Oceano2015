function [uit_out, uit_se_out] = determine_uit(u_fr, flux, rho_a)

plot_switch = false;
uit_min_rel_err = 0.08; %minimum relative error for the impact threshold shear velocity
fit_margin = 0.20; %sets fit margin for uit. that is, the data points within (1+fit_margin*)uit will be included in the linear fit for uit
max_fit_points = 1000; %maximum data points used in linear uit fit
max_ufr_diff = 0.50; %data points within this margin from uit will be considered for the fit for uit
min_fit_points = 3; %minimum number of data points for linear uit fit to be reasonable
second_order_fit = true; %sets whether the threshold is fit using a 2nd-order Taylor expansion instead of a 1st order (linear)
min_ufr = 0.16; %minimum value of u_fr for which measurements can be considered above threshold in m/s
rho_a_0 = 1.2250; %density at standard atmosphere
use_fit_uit = false; %sets whether to use fit_uit, which uses explicit non-linear squares fitting but gives unrealistic results for some cases

n = 1;
uit(1) = min(u_fr);
convergence = false;
tolerance = 0.0001;

while (convergence == false)
	n=n+1;
    i_end_ufr = size(u_fr,2);
    i_end_fit_margin = size(find(u_fr<(1+fit_margin)*max([min(u_fr),uit(n-1)])),2);
    %i_end_fit_margin = size(find(u_fr<(1+fit_margin)*uit(n-1)),2);
    i_end_ufr_diff = size(find(u_fr<min(u_fr)+max_ufr_diff),2);
    i_end = min([i_end_ufr,i_end_fit_margin,i_end_ufr_diff,max_fit_points]);
    if (i_end < min_fit_points)
        error('fit does not converge')
    end
    x = u_fr(1:i_end);
    y = flux(1:i_end);
    y_err = 0; %unknown error
    if (second_order_fit)
        if (use_fit_uit)
            [uit(n), uit_se(n), C_min] = fit_uit(x, y, 0);
        else
            [intercept, slope, SE_intercept, SE_slope, covariance] = linearfit (x.^2, y, y_err); %intercept = a, and slope = b in the notation of Bevington and Robinson (2003)
            if (slope > 0) %making sure that flux increases with u_Fr
                uit(n) = (sqrt(-intercept/slope)+uit(n-1))/2;
                uit_se(n) = sqrt(SE_intercept^2/slope^2+SE_slope^2*intercept^2/slope^4-2*covariance^2*intercept/slope^3);
            else
                error('Flux vs u_fr does not give good regression')
            end        
        end
        if(imag(uit(n))~=0)
        	uit(n)=uit(n-1)-0.02;
        end
        if (uit(n)-uit(n-1)>0.02)
            uit(n)=uit(n-1)+0.02;
        elseif (uit(n-1)-uit(n)>0.02)
        	uit(n)=uit(n-1)-0.02;
        end           
    else
        [intercept, slope, SE_intercept, SE_slope, covariance] = linearfit (x, y, y_err); %intercept = a, and slope = b in the notation of Bevington and Robinson (2003)
        if (slope > 0) %making sure that flux increases with u_Fr
            uit(n) = max((-intercept/slope+uit(n-1))/2,min_ufr*sqrt(rho_a_0/rho_a));
            if (uit(n) == min_ufr*sqrt(rho_a_0/rho_a))
                uit_se(n) = uit(n)*uit_min_rel_err;
            else
                uit_se(n) = sqrt(SE_intercept^2/slope^2+SE_slope^2*intercept^2/slope^4-2*covariance^2*intercept/slope^3);
            end
        else
            error('Flux vs u_fr does not give good regression')
        end        
    end
    if (n > 20)
        uit(n) = mean(uit(20:n));
        abs((uit(n)-uit(n-1))/uit(n-1));
        1;
    elseif (n > 10)
        uit(n) = mean(uit(10:n));
    end
    if (abs((uit(n)-uit(n-1))/uit(n-1))<tolerance)
        convergence = true;
        if (plot_switch)
            x_plot = uit(n-1):0.001:1.01*max(x);
            if (second_order_fit)
                if(use_fit_uit)
                    figure(1); clf; plot(x,y,'ko','MarkerSize',10,'LineWidth',2); hold on; plot(x_plot,C_min*(x_plot.^2-uit(n)^2),'b--','LineWidth',3);
                else
                    figure(1); clf; plot(x,y,'ko','MarkerSize',10,'LineWidth',2); hold on; plot(x_plot,intercept+slope*x_plot.^2,'b--','LineWidth',3);
                end
                axis([0.02*min(floor(min(50*x)),floor(50*uit(n))) 0.02*ceil(max(50*x)) 0 1.1*max(y)]);
                set(gca,'FontSize',22); 
                xlabel('Soil friction speed, \itu\rm_{*} (m/s)','FontSize',26); 
                ylabel('Saltation flux, \itQ\rm (impacts per second)','FontSize',24); 
            else
                figure(1); clf; semilogy(x,y,'s'); hold on; plot(x_plot,intercept+slope*x_plot,'r');
            end
            
        end
    end   
end
uit_out = max(uit(n),sqrt(rho_a_0/rho_a)*min_ufr); %uit can't be smaller than min_ufr
uit_se_out = uit_min_rel_err*uit_out;