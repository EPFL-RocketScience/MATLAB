% Test Rocket class

close all;
clear all;

% parameters
rho_carbon = 1550; % kg/m^3
data = load('HyperTEK_Mgrain.mat');
thrustCurve = data.Tdat;
bt = 10; % burn time
K = 1.1; % coefficient correctif de force aerodynamique normale

% Create a rocket
R = Rocket();

% nose(L, D, e, rho)
R.nose(0.5, 0.15, 0.002, rho_carbon);
% tail(D1, D2, L, e, z, rho)
R.tail(0.15, 0.5, 0.2, 0.002, 2.5, rho_carbon);
% stage(id, z, L, Dout, e, rho)
R.stage('tube', 0.5, 2, 0.15, 0.002, rho_carbon);
R.stage('motor', 2, 0.5, 0.15, 0.002, rho_carbon);
% motor(id, z, m, D, L, thrustCurve, bt)
m = @(t) (10 - 0.4*t).*(t<=bt)+6.*(t>bt);
R.motor('grain', 2, m, 0.098, 0.5, thrustCurve, bt);
m = @(t) (5 - 0.5*t).*(t<=bt)+0.*(t>bt);
R.motor('tank', 2, m, 0.098, 0.5, thrustCurve, bt);
% payload(id, z, m, L, D)
R.payload('main', 0.8, 5, 0.3, 0.15);
% parachute(id, z, m, D, Cd)
R.parachute('main', 0.5, 1, 2, 1.5);
R.parachute('secondary', 0.5, 0.5, 0.5, 0.75);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(2, 4, 0.2, 0.3, 0.1, 0.1, 0.075, 0.005, rho_carbon);
% point(id, z, m)
R.point('tuyere', 2.5, 1);

% Check definition
R.check();

% Print specs
R.printSpecs();

% Plot stuff
t = linspace(0, max([R.Motor.bt])+10, 100);
for it = 1:length(t) 
   mass(it) = R.m(t(it));
   cm(it) = R.cm(t(it));
   Iz(it) = R.Iz(t(it));
   Ir(it) = R.Ir(t(it));
end
figure
subplot(4,1,1)
plot(t, mass);
title('Mass properties');
ylabel('m [kg]');
subplot(4,1,2)
plot(t, cm);
ylabel('CM [m]');
subplot(4,1,3)
plot(t, Iz);
ylabel('Iz [kg*m^2]');
subplot(4,1,4)
plot(t, Ir);
ylabel('Ir [kg*m^2]');
xlabel('t [s]');


alpha = linspace(-30, 30, 100)*pi/180;
figure
for alpha_i = 1:length(alpha)
   [CNa(alpha_i), zCP(alpha_i)] = R.aeroCoeff(alpha(alpha_i), 0, 0, K); 
end
subplot(2,1,1)
plot(alpha*180/pi, CNa.*alpha);
title('Aerodynamic properties at M = 0 & \theta = 0.')
ylabel('CN');
grid on;
subplot(2,1,2)
plot(alpha*180/pi, zCP);
ylabel('zCP [m]');
grid on;
xlabel('alpha [^\circ]');

