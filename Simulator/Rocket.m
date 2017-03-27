classdef Rocket < handle
    
    properties (SetAccess = private, GetAccess = public)
        Nose
        Stage
        Tail
        Fins
        Motor
        Payload
        Parachute
        Points
    end
    
    methods
        function obj =  Rocket()
            % Rocket
            % Constructor of Rocket, empty for now
        end
    end
    
    methods
        function nose(obj, L, D, e, rho, type)
            % nose
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - L     :
            %   - D     :
            %   - z     :
            %   - e     :
            %   - rho   :
            %   - type  :
            
            % Assignation des proprietes
            obj.Nose.L = L;
            obj.Nose.D = D;
            obj.Nose.z = z;
            obj.Nose.rho = rho;
            obj.Nose.type = type;
            
            % Calcule des proprietes de masse
            obj.Nose.m = rho*pi/3*(((D/2)^2*L)-(((D/2)-e)^2*(L-e)));
            obj.Nose.cm = 3*L/4;
            obj.Nose.Iz = 0;
            obj.Nose.Ir = 0;
            
            % Calcule des proprietes aerodynamiques
            obj.Nose.CN = 0; % coefficient aerodynamique normal
            obj.Nose.zCP = 0; % position du centre de pression relatif au haut du cone
        end
        
        function tail(obj, z, L, e, D1, D2, rho)
            % D1 : diametre du haut
            % D2 : diametre du bas
        end
        
        function stage(obj, id, z, L, Din, Dout, rho)
            tmpStage.id = id; % 
            tmpStage.L = L;
            tmpStage.Din = Din;
            tmpStage.Dout= Dout;
            tmpStage.rho = rho;
            
            
            
            obj.Stage = [obj.Stage tmpStage];
        end
        
        function motor(obj, z, m, thrustCurve, bt)
            if(~isa(m, 'function_handle'))
                error('La masse doit ?tre une fonction')
            end
            
            obj.Motor.z = z;
            obj.Motor.m = m;
            obj.Motor.ThrustCurve = ThrustCurve;
            obj.Motor.bt = bt;
            
        end
        
        function payload(obj, z, m, L, D)
            % payload
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   position du haut de la payload
            %   - m     :   masse
            %   - L     :   longueur de la payload
            %   - D     :   diametre de la paylaod
            
            % Assignation des proprietes
            obj.Payload.z = z;
            obj.Payload.m = m;
            obj.Payload.L = L;
            obj.Payload.D = D;
            
            % Calcule des proprietes de masse
            obj.Nose.m = m;
            obj.Nose.cm = (1/2)*L;
            obj.Nose.Iz = 0;
            obj.Nose.Ir = 0;
        end
        
        function parachute(obj, z, m, L, D)
            
        end
        
        function fins(obj, z, N, a, b, gamma, phi, e, rho)
            % fins
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   position du haut de la fins
            %   - N     :   nombre de fins
            %   - a     :   longueur de la petite base
            %   - b     :   longueur de la grande base
            %   - gamma :   angle � l'avant de la fins
            %   - phi   :   angle � l'arri�re de la fins
            %   - e     :   epaisseur de la fins
            %   - rho   :   densite 
           
            % Assignation des proprietes
            obj.Fins.z = z;
            obj.Fins.N = N;
            obj.Fins.a = a;
            obj.Fins.b = b;
            obj.Fins.gamma = gamma;
            obj.Fins.phi = phi;
            obj.Fins.e = e;
            obj.Fins.rho = rho;
            
            % Calcule des proprietes de masse
            obj.Fins.h = (b-a)/(1/tan(gamma)+1/tan(phi)); 
            obj.Fins.cm = (h/3)*(2*a+b)/(a+b);
            obj.Fins.m = rho*N*(1/2)*(a+b)*h*e;
        end
        
        function point(obj, z, m)
        
        end
        
        function getDataFromFile(obj, path)
            CADFILE = importdata(path);
            
            LN  = CADFILE.data(13)  ;    %length of nose  
            D  = CADFILE.data(9)   ; %diameter at base of nose  
            LB = CADFILE.data(10)    ; %length of Body tube 
            DT  =  CADFILE.data(13)   ;  %diameter at rear of transition  
            LT  = CADFILE.data(12)    ;%length of transition  
            XT  =  CADFILE.data(10)+obj.LN  ;   %distance from tip of nose to front of transition (=tube lenght + nose as transition at the end) 
            CR  =  CADFILE.data(1)    ; %fins root chord  
            CT  =  CADFILE.data(2) ;    %fins tip chord  
            S  =  (obj.CR-obj.CT)/(1/tand(CADFILE.data(4))+1/tand(CADFILE.data(3)))   ;%fins semispan  
            LF  =  (obj.S^2+(obj.CR/2-obj.CT/2-obj.S*tand(CADFILE.data(3)))^2)^.5  ;   %length of fin mid-chord line  
            XRT  = obj.S/tand(CADFILE.data(3))     ; %distance between fin root leading edge and fin tip leading edge parallel to body  
            XF  = obj.XT-obj.CR-CADFILE.data(8)     ;%distance from nose tip to fin root chord leading edge  
            NF  = CADFILE.data(6)      ;%number of fins  
            
        end
    end

end