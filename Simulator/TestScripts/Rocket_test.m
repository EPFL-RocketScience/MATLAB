% Test Rocket class

close all;
clear all;

% parameters
rho_carbon = 1550; % kg/m^3
data = load('HyperTEK_Mgrain.mat');
thrustCurve = data.Tdat;
bt = 10; % burn time

% Create a rocket
R = Rocket();

% nose(L, D, e, rho)
R.nose(0.5, 0.15, 0.002, rho_carbon);
% tail(D1, D2, L, e, z, rho)
R.tail(0.15, 0.5, 0.2, 0.002, 2.5, rho_carbon);
% stage(id, z, L, Dout, e, rho)
R.stage('tube', 0.5, 2, 0.15, 0.002, rho_carbon);
R.stage('motor', 2, 0.5, 0.15, 0.002, rho_carbon);
% motor(z, m, D, L, thrustCurve, bt)
m = @(t) (10 - 0.4*t).*(t<=bt)+6.*(t>bt);
R.motor(2, m, 0.098, 0.5, thrustCurve, bt);
% payload(id, z, m, L, D)
R.payload('main', 1, 5, 0.3, 0.15);
% parachute(id, z, m, D, Cd)
R.parachute('main', 0.5, 1, 2, 1.5);
R.parachute('secondary', 0.5, 0.5, 0.5, 0.75);
% fins(z, N, a, b, gamma, phi, e, rho);
R.fins(2, 4, 0.15, 0.35, 45, 90, 0.005, rho_carbon);
% point(id, z, m)
R.point('tuyere', 2.5, 1);

% Calculate stuff
R.update();

% Plot stuff
t = linspace(0, R.Motor.bt+10, 100);
figure
subplot(3,1,1)
plot(t, R.Mass(t));
ylabel('m [kg]');
subplot(3,1,2)
plot(t, R.CM(t));
ylabel('CM [m]');
subplot(3,1,3)
plot(t, R.Iz(t));
ylabel('Iz [kg*m^2]');
xlabel('t [s]');

