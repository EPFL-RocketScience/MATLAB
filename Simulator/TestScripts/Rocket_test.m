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
R.payload('main', 0.8, 5, 0.3, 0.15);
% parachute(id, z, m, D, Cd)
R.parachute('main', 0.5, 1, 2, 1.5);
R.parachute('secondary', 0.5, 0.5, 0.5, 0.75);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(2, 4, 0.25, 0.35, 0.1, 0.15, 0.075, 0.005, rho_carbon);
% point(id, z, m)
R.point('tuyere', 2.5, 1);

% Calculate stuff
R.update();

% Plot stuff
t = linspace(0, R.Motor.bt+10, 100);
figure
subplot(4,1,1)
plot(t, R.Mass(t));
ylabel('m [kg]');
subplot(4,1,2)
plot(t, R.CM(t));
ylabel('CM [m]');
subplot(4,1,3)
plot(t, R.Iz(t));
ylabel('Iz [kg*m^2]');
subplot(4,1,4)
for ti = 1:length(t)
    Ir(ti) = R.Ir(t(ti)); 
end
plot(t, Ir);
ylabel('Ir [kg*m^2]');
xlabel('t [s]');


