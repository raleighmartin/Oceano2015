%% From Jasper: "For the uncertainty estimation I don't have a very clear
% script, but see lines 122 and 123 below.  These lines do a finite
% difference calculation of the second derivative of xi square towards the
% fitting parameter to obtain the error in the parameter. This is based on
% p. 146-148 in Bevington and Robinson (2003), which you can borrow from me
% again. This method only works when you have either one fitting parameter,
% or your fitting parameters are uncorrelated. Otherwise, things get a bit
% complex...
%%

function [uit_min, uit_SE, C_min] = fit_uit (x, y, a)

global plot_switch Fd_rel_err uit_min_rel_err accurate uit_log_fit

local_plot_switch=false;

y_SD = Fd_rel_err*y;
%y_SD = 0.1*mean(y)*ones(size(y,2),1);
log_y = log10(y);
log_y_SD = log10(y_SD); %assume equal errors in log-space

x_int = 0.01:0.01:50;

if (plot_switch)
    figure(2); clf; errorbar(x,y,y_SD,'s');
end

delta_C = 0.1;
if(accurate)
    delta_C = delta_C/5;
end
C_est = max(y)/(max(x)^2-min(x)^2);
C = 10.^(log10(C_est)-1:delta_C:log10(C_est)+1);

delta_uit = 0.01;
if(accurate)
    delta_uit = delta_uit/5;
end
uit = 0.1*min(x):delta_uit:1.0*min(x);

for i=1:size(C,2)   
    i;
    for j=1:size(uit,2)
        j;
        for k=1:size(x,2) 
            if (x(k)<uit(j))
                calc(i,j,k)=0.01;
            else                
                calc(i,j,k)=C(i)*(x(k)^2-uit(j)^2)*(x(k)/uit(j))^a;                
            end
            if (uit_log_fit)
                diff(i,j,k)=(1/log_y_SD(k))^2*(log10(calc(i,j,k))-log_y(k))^2;                
            else
                diff(i,j,k)=(1/y_SD(k))^2*(calc(i,j,k)-y(k))^2;
            end            
        end
        chi_square(i,j)=sqrt(sum(diff(i,j,:)));   
    end
end

if (plot_switch)
    A=ones(size(C,2),size(uit,2)); for p=1:size(uit,2) A(:,p)=C; end
    B=ones(size(C,2),size(uit,2)); for p=1:size(C,2) B(p,:)=uit; end
    figure(4); clf; mesh(A,B,chi_square); %contour(A,B,chi_square);
end
[r,j_min]=min(min(chi_square));
[r,i_min]=min(chi_square(:,j_min));
if (chi_square(i_min,j_min)-min(min(chi_square))~=0)
    error('something wrong')
end
chi_square(i_min,j_min);
C_min = C(i_min);
uit_min = uit(j_min);

delta_C_old = delta_C;
delta_uit_old = delta_uit;
delta_C=delta_C_old/5;
delta_uit=delta_uit_old/5;

clear C; clear uit; clear chi_square

C = 10.^(log10(C_min)-2*delta_C_old:delta_C:log10(C_min)+2*delta_C_old);
uit = (uit_min-2*delta_uit_old):delta_uit:(uit_min+2*delta_uit_old);

for i=1:size(C,2)   
    i;
    for j=1:size(uit,2)
        j;
        for k=1:size(x,2)
            if (x(k)<uit(j))
                calc(i,j,k)=0.01;
            else
                calc(i,j,k)=C(i)*(x(k)^2-uit(j)^2)*(x(k)/uit(j))^a;
            end
            if (uit_log_fit)
                diff(i,j,k)=(1/log_y_SD(k))^2*(log10(calc(i,j,k))-log_y(k))^2;                
            else
                diff(i,j,k)=(1/y_SD(k))^2*(calc(i,j,k)-y(k))^2;
            end
        end
        chi_square(i,j)=sqrt(sum(diff(i,j,:)));        
    end
end

if (plot_switch)
    A=ones(size(C,2),size(uit,2)); for p=1:size(uit,2) A(:,p)=C; end
    B=ones(size(C,2),size(uit,2)); for p=1:size(C,2) B(p,:)=uit; end
    figure(4); clf; mesh(A,B,chi_square); %contour(A,B,chi_square);
end
[r,j_min]=min(min(chi_square));
[r,i_min]=min(chi_square(:,j_min));
if (chi_square(i_min,j_min)-min(min(chi_square))~=0)
    error('something wrong')
end
chi_square(i_min,j_min);
%C_min = C(i_min)
%uit_min = uit(j_min)
if (j_min==1 || j_min==size(uit,2))
    uit_min = uit(j_min);
    uit_SE = uit_min*uit_min_rel_err;
else       
    uit_min = uit(j_min+1)-(delta_uit)*((chi_square(i_min,j_min+1)-chi_square(i_min,j_min))/(chi_square(i_min,j_min-1)-2*chi_square(i_min,j_min)+chi_square(i_min,j_min+1))+0.5);
    uit_SE = (delta_uit)*sqrt(2/(chi_square(i_min,j_min-1)-2*chi_square(i_min,j_min)+chi_square(i_min,j_min+1)));
end
if (i_min==1 || i_min==size(C,2))
    C_min = C(i_min);
else
    C_min = C(i_min+1)-(delta_C)*((chi_square(i_min+1,j_min)-chi_square(i_min,j_min))/(chi_square(i_min-1,j_min)-2*chi_square(i_min,j_min)+chi_square(i_min+1,j_min))+0.5);
end
%C_SE = (delta_C)*sqrt(2/(chi_square(i_min-1,j_min)-2*chi_square(i_min,j_min)+chi_square(i_min+1,j_min)));

%calc_nonlinear_range
%uit_SE = (uit_upper - uit_lower)/2
%C_SE = (C_upper - C_lower)/2
%log10(C_min)-log10(C_lower);
%C_SE/(C_min*log(10));

x_plot = 0.99*min(x):(max(x)-min(x))/10:1.01*max(x);

if (local_plot_switch)
    if (log_fit)
        figure(2); clf; semilogy(x_plot, C_min*(x_plot.^2-uit_min^2).*(x_plot/uit_min).^a,'r'); hold on; errorbar(x,y,y_SD,'s'); 
    else
        figure(2); clf; semilogy(x_plot, C_min*(x_plot.^2-uit_min^2).*(x_plot/uit_min).^a,'r'); hold on; errorbar(x,y,y_SD,'s'); 
    end
end

%calculating the coefficient of determination, R-squared
for k=1:size(x,2)     
    y_calc(k)=C_min*x(k)^uit_min;
    log_y_calc(k)=log10(C_min*x(k)^uit_min);
end

SST = sum((y-mean(y)).^2);
SSR = sum((y_calc-mean(y)).^2);
SSE = sum((y-y_calc).^2);
SSR/SST;
R_squared = 1-SSE/SST;

log_SST = sum((log_y-mean(log_y)).^2);
log_SSR = sum((log_y_calc-mean(log_y)).^2);
log_SSE = sum((log_y-log_y_calc).^2);
log_SSR/log_SST;
log_R_squared = 1-log_SSE/log_SST;

1;

