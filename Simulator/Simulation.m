classdef Simulation < matlab.System
    % Simulation
    %  
    %       - V0    : vitesse du vent [m/s]
    %       - tfin  : temps de simulation maximal
    %       - phi0  : angle de d?part de la rampe en [rad]
    %       - l_ramp: longueure de la rampe de lancement [m]
    %       - tquer : times at which special calculations should be done (flexion)
    %       - xquer : positions along rocket where special calculation values
    %                 are requested
    %       - G     : Choix des graphes � afficher :
    %   Figure1     : Mass properties
    %   Figure2     : Aerodynamic properties at M = 0 & \theta = 0
    %   Figure3     : Drawing of the rocket with CM and CP
    %   Figure4     : Altitude
    %   Figure5     : Rocket Trajectory
    %   Figure6     : Rocket angle
    %   Figure7     : Angle of attack
    %   Figure8     : Stabiltie statique lors du vol
    %   Figure9     : Flexion load
    

    properties(SetAccess = public, GetAccess = public)
        L_ramp
        phi0
        v_vent
        tfin
        tquer
        xquer
        G
        
    end

    methods(Access = public)
        function obj =  Simulation()
            % Simulation
            % Constructor of Simulation, empty for now
        end
    end
    methods
        function inititSimulation(obj, L_ramp, phi0,  v_vent, tfin, tquer, xquer)
            % inititSimulation
            % Initialise la simulation
            % INPUTS
            %   - L_ramp    : Longueur de la rampe de lancement
            %   - pho0      : Angle de la rampe de lancement par rapport � la verticale 
            %   - tfin      : Temps de fin de la simulation au maximum
            %   - tquer     : times at which special calculations should be done (flexion)
            %   - xquer     : positions along rocket where special calculation values
            %                 are requested
            %   - v_vent    : Vitesse du vent
            
            % Assignation des proprietes
            obj.L_ramp = L_ramp;
            obj.phi0 = phi0;
            obj.tfin = tfin;
            obj.tquer = tquer;
            obj.xquer = xquer; 
            obj.v_vent = v_vent;
            obj.G = [ 0 0 0 0 0 0 0 0 0];%Dessin des graphes
        end
        
        function setGraph(obj, G)
           
            obj.G = G
                    
        end
        
        function printSpecs(obj)
           % affiche les specs de l environment
           
           display('*********************');
           display('Environment Specifications');
           display('*********************');
           
           % length ramp
           display(['* Longueur rampe   : ' num2str(obj.L_ramp) ' [m]']);
           
           % angle inclinaison
           display(['* Angle rampe   : ' num2str(obj.phi0) ' [rad]']);
           
           % wind
           display(['* Vitesse vent : ' num2str(obj.v_vent) ' [m/s]']);       
                    
        end
        
    end
    
end
