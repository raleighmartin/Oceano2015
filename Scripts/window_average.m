%time_out gives timeseries of first time for each window

function [data_out, time_out] = window_average(data_in, t, dt_avg)

L_data = length(data_in);
dt = mode(diff(t));
L_chunk = round(dt_avg/dt);
N_chunk = floor(L_data/L_chunk);

data_out = zeros(N_chunk,1);
for n=1:N_chunk
    data_out(n) = mean(data_in((n-1)*(1:L_chunk)+1));
end

time_out = (t(1)+(0:(N_chunk-1))*dt_avg)';