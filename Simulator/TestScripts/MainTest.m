% Main test for simulation of any rocket and mission

close all;
clear all;
clc;

% Creation des deux objets de la simulation
R = RocketIREC();
S = SimulationIREC();

% Simulate part
results = Simulate( R, S);

% DrawGraphs
DrawGraph(R, S, results);

% Abaqus_Analysis(20,10,0.2,max(abs(results.accel.Vz_dot)),R);