function S = SimulationZurich()
% Environment de zurich
%   Iinitalisation des variable
    tfin = 11;
    phi0 = 0;
    tquer = [1 3 5 7 8];
    xquer = 0;%linspace(0, R.Tail.z + R.Tail.L, 10);
    L_ramp = 2;
    v_vent = 10;
    
    %create a envirnment
    S = Simulation();
    %initiat environement
    inititSimulation(S, L_ramp, phi0, tfin, tquer, xquer, v_vent);
    %Afficher environment
    S.printSpecs();
    
end

