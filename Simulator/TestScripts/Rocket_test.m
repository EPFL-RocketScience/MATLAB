% Test Rocket class

close all;
clear all;

% parameters
rho_carbon = 1550; % kg/m^3
data = load('HyperTEK_Mgrain.mat');
thrustCurve = data.Tdat;

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
R.motor(2, @(t) 10, 0.098, 0.5, thrustCurve, 10);
% payload(id, z, m, L, D)
R.payload('main', 1, 5, 0.3, 0.15);
% parachute(id, z, m, D, Cd)
R.parachute('main', 0.5, 1, 0.15, 1.5);
% fins(z, N, a, b, gamma, phi, e, rho);
R.fins(2, 4, 0.15, 0.35, 45, 90, 0.005, rho_carbon);
% point(id, z, m)
R.point('tuyere', 2.5, 1);