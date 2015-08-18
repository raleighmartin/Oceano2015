function [ind_init_windows, ind_cess_windows, t_init, t_cess] = fluxevent_detection(q,t,T_min,T_window)

N_t = length(t); %determine duration of timeseries
dt = mode(diff(t)); %determine time increment
N_min = round(T_min/dt); %number of timesteps of continuous transport/rest for initiation/cessation events 
N_window = round(T_window/dt); %total number of timesteps to be included in event

%create timeseries of current flux detection state, =1 if flux detected, =0 if no flux detected
q_current = (q~=0); 

%create timeseries of previous flux states, =proportion of previous N_min timesteps (including current) with flux detected
q_previous = zeros(N_t,1);
q_previous(1:N_min)=cumsum(q_current(1:N_min))./(1:N_min)'; %deal with first few entries
for i=0:(N_min-1)
    q_previous(N_min:N_t) = q_previous(N_min:N_t)+(q_current((N_min-i):(N_t-i)))/N_min;
end

%create timeseries of next flux states, =proportion of next N_min timesteps (including current) with flux detected
q_next = zeros(N_t,1);
for i=0:(N_min-1)
    q_next(1:(N_t-N_min+1)) = q_next(1:(N_t-N_min+1))+(q_current((1+i):(N_t-N_min+1+i)))/N_min;
end
q_next(N_t:-1:(N_t-N_min+1))=cumsum(q_current(N_t:-1:(N_t-N_min+1)))./(1:N_min)'; %deal with last few entries

%determine indices of dead (full) transport, i.e. times when q_next = 0 (1)
q_next = round(q_next*N_min)/N_min; %kluge to deal with rounding issues
ind_full = find(q_next==1); %find indices of times with 'full' transport
ind_dead = find(q_next==0); %find indices of times with 'dead' transport

%locate initiation and cessation events based on alternations between dead and full transport
ind_init = []; %initialize initiation indices list
ind_cess = []; %initialize cessation indices list

%set initial flux state (q_state =0 for no transport, =1 for transport)
q_state = round(mean(q_current(1:N_window)));

%index t_ind will record progress through timeseries.
%first possible i is N_event, since we only want to identify
%initiation/cessation with +/- N_event timesteps before and after
t_ind = N_window;
while t_ind <= (N_t-N_window) %keep going until t_ind reaches N_window
    if q_state == 0 %if current state is no transport, look for next full transport
        ind_next = min(ind_full(ind_full>t_ind)); %index of next full transport
        if isempty(ind_next)||((ind_next+N_window)>N_t) %quit cycle if none is to be found, or if next event will cause window to go past timeseries
            break
        else %otherwise, add this to list and update time and state
            ind_init = [ind_init; ind_next];
            q_state = 1;
            t_ind = ind_init(end);
        end
    elseif q_state == 1 %if current state is transport, look for next dead transport
        ind_next = min(ind_dead(ind_dead>t_ind)); %index of next dead transport
        if isempty(ind_next)||((ind_next+N_window)>N_t) %quit cycle if none is to be found, or if next event will cause window to go past timeseries
            break
        else %otherwise, add this to list and update time and state
            ind_cess = [ind_cess; ind_next];
            q_state = 0;
            t_ind = ind_cess(end);
        end
    end
end

%get associated times
t_init = t(ind_init);
t_cess = t(ind_cess);

%create matrix with indices for windows of initiation and cessation
if ~isempty(ind_init);
    ind_init_windows = (ones(size(ind_init))*(-N_window:N_window))+(ind_init*ones(1,N_window*2+1));
else
    ind_init_windows = [];
end
if ~isempty(ind_cess);
    ind_cess_windows = (ones(size(ind_cess))*(-N_window:N_window))+(ind_cess*ones(1,N_window*2+1));
else
    ind_cess_windows = [];
end

% %determine transport states (combines info about currently detected and previous fluxes)
% q_state = zeros(N_t,1); %initialize transport states
% q_state((q_current==1)&(q_previous>0)) = 1;
% q_state((q_current==0)&(q_previous==1)) = 1;
% 
% %set transport state for first N_min timesteps
% if sum(q_current(1:N_min))==0
%     q_state(1:N_min) = 0; %0 if no transport detected in first few timesteps
% else
%     q_state(1:N_min) = 1; %1 if transport detected in first few timesteps
% end
% 
% 
% 
% N_start = find(diff(q_current)==1)+1; %timesteps for initiation of motion
% N_stop = find(diff(q_current)==-1)+1; %timesteps for cessation of motion
% 
% %find significant initiation and cessation events
% if min(N_start)<min(N_stop) %first event is start
%     if max(N_start)>max(N_stop) %last event is start
%         trans_N = N_stop - N_start(1:end-1);
%         hiatus_N = N_start(2:end) - N_stop;
%     else %last event is stop
%         trans_N = N_stop - N_start;
%         hiatus_N = N_start(2:end) - N_stop(1:end-1);
%     end
%     start_long_ind = intersect(find(trans_N>=N_min),(find(hiatus_N>=N_min)+1));
%     stop_long_ind = intersect(find(trans_N>=N_min),find(hiatus_N>=N_min));    
% else %first event is stop
%     if max(N_start)>max(N_stop) %last event is start
%         trans_N = N_stop(2:end) - N_start(1:end-1);
%         hiatus_N = N_start - N_stop;
%     else %last event is stop
%         trans_N = N_stop(2:end) - N_start;
%         hiatus_N = N_start - N_stop(1:end-1);
%     end
%     start_long_ind = intersect(find(trans_N>=N_min),find(hiatus_N>=N_min));
%     stop_long_ind = intersect((find(trans_N>=N_min)+1),find(hiatus_N>=N_min));
% end
% N_start_long = N_start(start_long_ind);
% N_stop_long = N_stop(stop_long_ind);
% 
% %remove values too close to edge of timeseries
% N_start_keep_ind = intersect(find(N_start_long>N_window),find(N_start_long<(N_segment-N_window)));
% N_stop_keep_ind = intersect(find(N_stop_long>N_window),find(N_stop_long<(N_segment-N_window)));
% N_start_long = N_start_long(N_start_keep_ind);
% N_stop_long = N_stop_long(N_stop_keep_ind);