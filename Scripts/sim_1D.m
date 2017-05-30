% 1D simulation of rocket with thrust curve and mass variation
clear all
close all

% System variables
Cd = 0.55;       % estimated drag coefficient
mi = 27.765;        % initial mass [kg]
mp = 7.183;     % propellant mass
mf = mi-mp;     % final mass
bt = 5.4;        % burn time [s]
S = pi*(0.17)^2/4;   % exposed section [m^2]

% Sim parameters
tSpan = [0, 30];
x0 = [0;0];

data = load('Thrust_Curves/Aerotech_N2220.mat');    % Thrust Curve data Tdat(:,1) = time [s], Tdat(:,2) = Thrust [N]
Tdat = data.data;
%Tdat = data.Tdat;
Tdat = Tdat(find(~isnan(Tdat(:,2))), :);

% Add values to thrust curve at initial time and final time 
if(Tdat(1,2) == 0)
    Tdat = [Tdat; bt+0.1, 0; tSpan(2), 0];
else
    Tdat = [0, Tdat(1,2); Tdat; bt+0.1, 0; tSpan(2), 0];
end

% assume linear mass variation
m = @(t) ((mf-mi)/bt * t + mi)*(t<=bt) + mf*(t>bt);
m_dot = @(t) (mf-mi)/bt*(t<=bt);

% Interpolate Thrust curve from Thrust data
Ft = @(t) interp1(Tdat(:,1)', Tdat(:,2)', t);

% Simulate 

odefun = @(t, x) xdot_1D(t,x,m, m_dot, Ft, Cd, S);
[t, x] = ode45(odefun, tSpan, x0);

[T, a, p, rho] = stdAtmos(x(:,1));

display(['Apogee: ' num2str(max(x(:,1) - x0(1))) ' m'])
display(['Mach max: ' num2str(max(x(:, 2)./a))])
display(['gmax: ' num2str(max(diff(x(:, 2))./diff(t)/9.8)) ' m/s^2']);

figure

plot(t, x(:, 1)-x0(1));
title 'Altitude';
ylabel 'z [m]';
xlabel 't [s]';
box off;

figure

subplot(4,1,1)
plot(t, x(:, 1)-x0(1));
title 'Altitude';
ylabel 'z [m]';

subplot(4,1,2);
plot(t, x(:, 2)./a);
title 'Mach'
ylabel 'Ma = v/a [-]'

subplot(4,1,3);
plot(t(1:end-1), diff(x(:, 2))./diff(t)/9.8);
title 'g'
ylabel 'g = a/g_0 [-]'

subplot(4,1,4)
plot(t, Ft(t));
title 'Thrust'
xlabel 't [s]'
ylabel 'Thrust [N]'



