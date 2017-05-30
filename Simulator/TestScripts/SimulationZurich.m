function S = SimulationZurich()
% Environment de Zurich
% INPUT
%   - R     : Rocket
%   - v_vent: un vent positif va de gauche ? droite
%                     negatif va de droite ? gauche


    % Iinitalisation des variable
    tfin = 11;
    phi0 = 0;
    v_vent = 10;
    tquer = [1 3 5 7 8];
    nquer = 10;
    L_ramp = 2;
   
    %initiate environement
    S = Simulation(L_ramp, phi0, v_vent, tfin, tquer, nquer);
    
    % Choix des graphes
    S.setGraph('mass', 1);
    S.setGraph('aero', 1);
    S.setGraph('rocket', 1);
    S.setGraph('altitude', 1);
    S.setGraph('trajectory', 1);
    S.setGraph('verticalAngle', 1);
    S.setGraph('attackAngle', 1);
    S.setGraph('stability', 1);
    % S.setGraph('flexion', 1);
    
    %Afficher environment
    S.printSpecs();
    
end

