% Test Rocket class and simulator

close all;
clear all;
clc;

% Creation des deux objets de la simulation
R = RocketJuju();
S = SimulationZurich(R);
G = [ 1 1 1 1 1 1 1 1 1];%Dessin des graphes

% Simulate part
[tsim, Xsim, alpha, calibre, T, M] = Simulate( R, S);

% DrawGraphs
DrawGraph( tsim, Xsim, R, S, alpha, calibre, T, M, G);
