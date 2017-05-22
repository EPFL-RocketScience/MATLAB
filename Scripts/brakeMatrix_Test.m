clear all;
close all;

% supress warnings
warning('off', 'Altitude:max');
warning('off', 'Altitude:min');

Cd0 = 0.55;     % estimated drag coefficient
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

xt = 200; % target

[Z2, X10, V10, qual, err, err_max] = brakeMatrix( m0, mp, D, Arefb, Cd0, Cdb, TC, BT, xt);

display('Quality of deployment matrix:');
display([num2str(100*qual) '% unreachable values']);

figure; hold on;
title(['Z2 : deployment altitudes for target = ' num2str(xt) ' m']);
surf(X10, V10, Z2);
colorbar;
xlabel('X1_0 [m]');
ylabel('V1_0 [m/s]');