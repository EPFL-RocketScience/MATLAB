% Main test for simulation of any rocket and mission

close all;
clear all;
clc;

% Creation des deux objets de la simulation
R = RocketJuju();
S = SimulationZurich(R);

% Simulate part
[tsim, Xsim, alpha, calibre, T, M] = Simulate( R, S);

% DrawGraphs
DrawGraph( tsim, Xsim, R, S, alpha, calibre, T, M);
