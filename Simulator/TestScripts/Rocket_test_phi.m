% Test Rocket class

close all;
clear all;

R = TestRocket();

% time step simulation
tspan = [0 10];
[tsim, Xsim] = ode45(@(t, x) stateEquation_phi(t, x, R, [0.1; 10; 0], 1), tspan, [0, 0]);
figure; hold on;
title('Rocket angle')
plot(tsim, Xsim(:,1));