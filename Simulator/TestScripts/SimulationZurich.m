function S = SimulationZurich( R)
% Environment de zurich
%   Iinitalisation des variable
    tfin = 11;
    phi0 = 0;
    v_vent = 10;
    tquer = [1 3 5 7 8];
    xquer = linspace(0, R.Tail.z + R.Tail.L, 10);%Pas JOLI DE REMPRENDRE R ICI
    L_ramp = 2;
    
    %create a envirnment
    S = Simulation();
    %initiat environement
    inititSimulation(S, L_ramp, phi0, v_vent, tfin, tquer, xquer);
    %Afficher environment
    S.printSpecs();
    
end

