classdef Environment < matlab.System
    % 

    properties(SetAccess = public, GetAccess = public)
        L_ramp
        phi0
        tfin
        tquer
        xquer
        v_vent
        
    end

    methods(Access = public)
        function obj =  Environment()
            % Environment
            % Constructor of Environment, empty for now
        end
    end
    methods
        function inititEnvironment(obj, L_ramp, phi0, tfin, tquer, xquer, v_vent)
            % inititEnvironment
            % Initialise l'environnement
            % INPUTS
            %   - L-ramp    : Longueur de la rampe de lancement
            %   - pho0      : Angle de la rampe de lancement par rapport à la verticale 
            %   - tfin      : Temps de fin de la simulation au maximum
            %   - tquer     :
            %   - xquer     : 
            %   - v_vent    : Vitesse du vent
            
            % Assignation des proprietes
            obj.L_ramp = L_ramp;
            obj.phi0 = phi0;
            obj.tfin = tfin;
            obj.tquer = tquer;
            obj.xquer = xquer; 
            obj.v_vent = v_vent;

        end
        
        function printSpecs(obj)
           % affiche les specs de l'environment
           
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
