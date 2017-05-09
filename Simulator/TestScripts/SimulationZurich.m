function S = SimulationZurich( R)
% Environment de Zurich
% INPUT
%   - R     : Rocket


%   Iinitalisation des variable
    tfin = 11;
    phi0 = 2*pi/360*45;
    v_vent = 10;
    tquer = [1 3 5 7 8];
    xquer = linspace(0, R.Tail.z + R.Tail.L, 10);%Pas JOLI DE REMPRENDRE R ICI
    L_ramp = 2;
    % Graphe
    %   Choix des graphes à afficher :
    %   Figure1     : Mass properties
    %   Figure2     : Aerodynamic properties at M = 0 & \theta = 0
    %   Figure3     : Drawing of the rocket with CM and CP
    %   Figure4     : Altitude
    %   Figure5     : Rocket Trajectory
    %   Figure6     : Rocket angle
    %   Figure7     : Angle of attack
    %   Figure8     : Stabiltie statique lors du vol
    %   Figure9     : Flexion load
    Figure1 = false;  
    Figure2 = false;  
    Figure3 = false;
    Figure4 = true;
    Figure5 = true;
    Figure6 = true;
    Figure7 = true;
    Figure8 = true;
    Figure9 = false;
    
    %create a envirnment
    S = Simulation();
    %initiat environement
    inititSimulation(S, L_ramp, phi0, v_vent, tfin, tquer, xquer);
    % Choix des graphes
    G = [Figure1 Figure2 Figure3 Figure4 Figure5 Figure6 Figure7 Figure8 Figure9];
    S.setGraph(G);
    
    %Afficher environment
    S.printSpecs();
    
end

