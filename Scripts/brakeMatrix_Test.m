clear all;
close all;

Cd0 = 0.85;     % estimated drag coefficient
Cdb = 1;        % estimated drag coefficient of brakes   
m0 = 2.264;     % initial mass [kg]
mp = 0.345;     % propellant mass [kg]
BT = 10;        % burn time [s]
D = 4*0.024;    % Diameter [m]
w = 0.05;       % Brake width [m]
l = 0.2;       % Brake length [m]
n = 4;          % Number of brake surfaces
Arefb = n*w*l;  % Brake reference area [m^2] 

data = load('Thrust_Curves/Aerotech_H123W.mat');    % Thrust Curve data Tdat(:,1) = time [s], Tdat(:,2) = Thrust [N]
TC = data.data;
TC = TC(find(~isnan(TC(:,2))), :);

Z2 = brakeMatrix( m0, mp, D, Arefb, Cd0, Cdb, TC, BT, 270);