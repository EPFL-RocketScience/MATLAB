clear all;
close all;

%% Define parameters of rocket and environment

M = 2.079;
rho = 1.2;
S = pi*(3*0.0247)^2/4;

%% Get Data

tmp = load('FlightData/VOL_JULIEN.mat');
data = tmp.data;
clear tmp

t = data(1, :);         % experiment time [ms]
h = data(12, :);        % Altitude (non calibrated) [m]
h_off = 745;

%% Get Flight Data

% flight time is generally no more than 2 min and time to apogee about 10
% sec
t_fly = 120000;     % [ms]
t_apo = 15000;      % [ms]

% flight happens around peak of altitude
[h_max, h_max_i] = max(h);

% get flight data
t_start = t(h_max_i)-t_apo;
t_end   = t(h_max_i)+t_fly-t_apo;
t_start_i = find(t>t_start, 1, 'first');
t_end_i = find(t<t_end, 1, 'last');
f_data  = data(:,t_start_i:t_end_i);

% get apogee
apo = max(h) - h_off;
display(['Apogee: ' num2str(apo)]);

% save only data from flight
t = f_data(1, :)-t_start;
h = f_data(12, :);
az = f_data(5, :);
ax = f_data(6, :);
ay = f_data(7, :);

%% Plot raw data
figure; hold on;
title('Altitude vs. time')
plot(t/1000, h);
xlabel('t [s]');
ylabel('h [m]');

figure; hold on;
title('Acceleration vs. time');
subplot(3, 1, 1);
plot(t/1000, az);
ylabel('a_z');
subplot(3, 1, 2);
plot(t/1000, ax);
ylabel('a_x');
subplot(3, 1, 3);
plot(t/1000, ay);
ylabel('a_y');
xlabel('t');

%% Plot altitude vs. time
% remove all altitudes lower than initial altitude and after apogee 
h_i = h(1);
[~, h_apogee_i] = max(h);
h_filt_i = find(h>=h_i+2);
h_filt_i = h_filt_i(find(h_filt_i<h_apogee_i));

figure; hold on;
title('Altitude vs. time (filtered)')
plot((t(h_filt_i)-t(h_filt_i(1)))/1000, h(h_filt_i)-h_i);
xlabel('t [s]');
ylabel('h [m]');

%% Get deceleration data
figure;
title('Select deceleration time span')
plot(t/1000, az);
ylabel('a_z');

% input from user
display('select start and end of deceleration time');
[t_in, az_in] = ginput(2);
t_in_i(1) = find(t>t_in(1)*1000, 1, 'first');
t_in_i(2) = find(t<t_in(2)*1000, 1, 'last');

% get data
t_decel = t(t_in_i(1):t_in_i(2));
az_decel = az(t_in_i(1):t_in_i(2));
h_decel = h(t_in_i(1):t_in_i(2));

% plot deceleration phase
figure; hold on;
title('Deceleration phase (along z axis)');
plot(t_decel, az_decel, 'b', 'DisplayName', 'data');
xlabel('t');
ylabel('a_z');

% fit second order function to data
p = fit(t_decel', az_decel', 'poly2');
p_az = [p.p1, p.p2, p.p3];
plot(t_decel, polyval(p_az, t_decel), 'r--', 'DisplayName', 'fit poly2');
legend show;

%% Get velocity data
h_decel = h_decel(find(h_decel>100));
t_decel_tmp = t_decel(find(h_decel>100));
figure; hold on;
title('Altitude during deceleration phase');
plot(t_decel_tmp, h_decel, 'b', 'DisplayName', 'data');
p = fit(t_decel_tmp', h_decel', 'poly2');
p_h = [p.p1, p.p2, p.p3];
plot(t_decel_tmp, polyval(p_h, t_decel_tmp), 'r--', 'DisplayName', 'fit poly2');
legend show;

% h_decel = h_decel(find(h_decel>100));
% t_decel_tmp = t_decel(find(h_decel>100));
% 
% v = (h_decel(1:end-2)-2*h_decel(2:end-1)+h_decel(3:end))/(mean(diff(t_decel_tmp))/1000)^2;
% v = v(find(v>0));
% 
% figure;
% plot(t_decel_tmp(find(v>0)), v);
%% Compute Drag coefficient

Cd = -2*M*polyval(p_az, t_decel)./(rho*S*(polyval([2*p_h(1) p_h(2)], t_decel)*1000).^2);

figure;
plot(t_decel, Cd);
%plot(t_decel, polyval([2*p_h(1) p_h(2)], t_decel)*1000);
%plot(t_decel, polyval(p_az, t_decel));
