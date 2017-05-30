classdef Simulation < matlab.System
    % Simulation
    % D?fintion des parametres de simulation et d'affichage.
    

    properties(SetAccess = public, GetAccess = public)
        L_ramp
        phi0
        v_vent
        tfin
        tquer
        nquer
        
        
        %%%%%%%%%%%%%%%%%%%
        % Graph attributes
        %%%%%%%%%%%%%%%%%%%
        G = []
        
    end

    methods
        function obj = Simulation(L_ramp, phi0,  v_vent, tfin, tquer, nquer)
            % Constructeur
            %
            % INPUTS
            %   - L_ramp    : Longueur de la rampe de lancement
            %   - phi0      : Angle de la rampe de lancement par rapport ? la verticale 
            %   - v_vent    : Vitesse du vent
            %   - tfin      : Temps de fin de la simulation au maximum
            %   - tquer     : times at which special calculations should be done (flexion)
            %   - nquer     : number of query points for flexion moment
            %                 plot
            
            % Assignation des proprietes
            obj.L_ramp = L_ramp;
            obj.phi0 = phi0;
            obj.tfin = tfin;
            obj.tquer = tquer;
            obj.nquer = nquer; 
            obj.v_vent = v_vent;
            
            % Initiation des parametres de graphes
            plot.id = 'mass';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'aero';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'rocket';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'altitude';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'trajectory';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'verticalAngle';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'attackAngle';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'stability';
            plot.value = 0;
            obj.G = [obj.G, plot];
            plot.id = 'flexion';
            plot.value = 0;
            obj.G = [obj.G, plot];
        end
        
        function setGraph(obj, varargin)
           % setGraph(obj, G) for?age manuel de la structure de param?tres
           % d'affichage.
           % 
           % setGraph(obj, id, value) attribution de la valeure 'value' au
           % param?tre avec l'identifiant 'id'. Si 'id' n'est pas un
           % identifiant connu, la m?thode g?n?rera une erreure. 
           
            if length(varargin) == 1
            
                obj.G = varargin{1};
                    
            elseif length(varargin) == 2
               
                id = varargin{1};
                value = varargin{2};
                
                index = find(strcmp(id, {obj.G.id}));
                
                if (~isempty(index) && length(index) == 1)
                    obj.G(index).value = value;   
                else
                    error(['id: ' id ' doesn''t correspond to an existing plot or more than one plot with the given id exists.'])
                end
                
            else
                error('unexpected number of arguments given to setGraph.')
            end
        end
        
        
        function val = getGraphValue(obj, id)
            % val = getGraphValue(obj, id) retourne la valeure du param?tre
            % de graphe avec l'identifiant 'id'.
            % Si 'id' n'est pas un identifiant connu, la m?thode 
            % g?n?rera une erreure.
            
            index = find(strcmp(id, {obj.G.id}));
                
            if (~isempty(index) && length(index) == 1)
                val = obj.G(index).value;  
            else
                error(['id: ' id ' doesn''t correspond to an existing plot or more than one plot with the given id exists.'])
            end
                
        end
        
        function printSpecs(obj)
           % printSpecs(obj) 
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
