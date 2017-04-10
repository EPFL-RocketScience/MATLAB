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
            %   - L     :   longueur du cone
            %   - D     :   Diamètre de la base du cone
            %   - e     :   epaisseur du cone
            %   - rho   :   densite
            %   - type  :   cone ou ogive
            
            % Assignation des proprietes
            obj.Nose.L = L;
            obj.Nose.D = D;
            obj.Nose.e = e;
            obj.Nose.rho = rho;
            obj.Nose.type = type;
            
            % Calcule intermediaire
            V = e*pi*sqrt(1+D^2/(4*L^2))*L*D/2;
                  
            % Calcule des proprietes de masse
            obj.Nose.m = rho*V;
            obj.Nose.cm = L*2/3;
            obj.Nose.Iz = rho/8*V*D^2;
            obj.Nose.Ir = rho*V/4*(D^2/4+2*L^2);
        
            % Calcule des proprietes aerodynamiques
            obj.Nose.CN = 0; % coefficient aerodynamique normal
            obj.Nose.zCP = 0; % position du centre de pression relatif au haut du cone
        end
        
        function tail(obj, D1, D2, L, e, z, rho)
            % INPUTS
            % D1 :  diametre du haut
            % D2 :  diametre du bas
            % L :   longueur du tail
            % e :   epaisseur de paroie
            % z :   postion par rapport au haut de la fus?e du haut du tail
            % rho : densite
            
            % Assignation des proprietes
            obj.Tail.D1 = D1;
            obj.Tail.D2 = D2;
            obj.Tail.L = L;
            obj.Tail.e = e;
            obj.Tail.z = z;
            obj.Tail.rho = rho;
            obj.Tail.type = type;
            
            %Calcule intermediaire
            hF = L*D2/(D1-D2);
            hTot = L*(D1/(D1-D2));
            vTot = e*pi*sqrt(1+D1^2/(4*L^2))*hTot*D1/2;
            vF = e*pi*sqrt(1+D2^2/(4*L^2))*hF*D2/2;
            Ixprime = rho/4*(vTot*(D1^2/4+2*hTot^2)-vF*(D2^2/4+2*hF^2));
            hPrime = hF+L*((D2+2*D1)/(3*(D2+D1)));
            
            % Calcule des proprietes de masse
            %on utilise un cone - un plus petit cone
            obj.Tail.m = rho*(vTot-vF);
            obj.Tail.cm = L*(2*D2+D1)/(3*(D2+D1));            
            
            R1 = D1/2;
            R2 = D2/2;
            m = (R2-R1)/L;
            obj.Tail.Iz = pi*rho/(10*m)*((R2^5-R1^5)-((R2-e)^5-(R1-e)^5));
            obj.Tail.Ir = Ix/2 + pi*rho*((m^2*L^5/5+1/2*m*L^4*R1+R1^2*L^3/3)-(m^2*L^5/5+1/2*m*L^4*(R1-e)+(R1-e)^2*L^3/3));
            
            % Calcule des proprietes aerodynamiques
            obj.Tail.CN = 0; % coefficient aerodynamique normal
            obj.Tail.zCP = 0; % position du centre de pression relatif
            
        end
        
        function stage(obj, id, z, L, Dout, e, rho)
            % Assignation des proprietes
            tmpStage.id = id; % 
            tmpStage.L = L;
            tmpStage.Dout = Dout;
            tmpStage.e = e;
            tmpStage.rho = rho;
            
            %Calcul intermediaire
            Din = Dout - e;%Si besoin
            
            % Calcule des proprietes de masse
            %probleme au niveau de l atribution aux differents etages
            obj.tmpStage.m = rho*pi*L*((Dout/2)^2-(Din/2)^2);
            obj.tmpStage.cm = L/2;
            obj.tmpStage.Iz = pi*rho*L*1/2*((Dout/2)^4-(Din/2)^4);
            obj.tmpStage.Ir = pi*rho*L*1/12*(3*((Dout/2)^4-(Din/2)^4)+L^2*((Dout/2)^2-(Din/2)^2));
            
            % Calcule des proprietes aerodynamiques ??? pas de propr. non ?
            obj.Nose.CN = 0; % coefficient aerodynamique normal
            obj.Nose.zCP = 0;
            
            obj.Stage = [obj.Stage tmpStage];
        end
        
        function motor(obj, z, m, D, L, thrustCurve, bt)
            if(~isa(m, 'function_handle'))
                error('La masse doit etre une fonction')
            end
            
            obj.Motor.z = z;
            obj.Motor.m = m;
            obj.Motor.ThrustCurve = ThrustCurve;
            obj.Motor.bt = bt;
            
            % Est-ce qu il y a des dimensions d un moteur en parametre de
            % la fonction ?
            % Calcule des proprietes de masse
            obj.Payload.cm = (1/2)*L;
            obj.Payload.Iz = m*D^2/8;
            obj.Payload.Ir = m/12*(3*(D/2)^2+L^2)+m*(L/2)^2;
            
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
            obj.Payload.cm = (1/2)*L;
            obj.Payload.Iz = m*D^2/8;
            obj.Payload.Ir = 0;
        end
        
        function parachute(obj, z, m, D, Cd)
            % parachutee
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   emplacement du parachute
            %   - m     :   masse du parachute
            %   - D     :   diamètre du parachute ouvert
            %   - Cd    :   Coefficient de trainee
            
            % Assignation des proprietes
            obj.Parachute.z = z;
            obj.Parachute.m = m;
            obj.Parachute.L = L;
            obj.Parachute.D = D;
            obj.Parachute.V = V;      
            
            % Calcule des proprietes de masse
            obj.Parachute.cm = 0;%car masse ponctuelle
            
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
            %   - gamma :   angle ? l'avant de la fins
            %   - phi   :   angle ? l'arri?re de la fins
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
            
            %Calcule intermediaire 
            obj.Fins.h = (b-a)/(1/tan(gamma)+1/tan(phi)); 
            
            % Calcule des proprietes de masse
            obj.Fins.cm = (h/3)*(2*a+b)/(a+b);
            obj.Fins.m = rho*N*(1/2)*(a+b)*h*e;
            obj.Fins.Iz = 0;
            obj.Fins.Ir = 0;
                   
            % Calcule des proprietes aerodynamiques
            obj.Fins.CN = 0; % coefficient aerodynamique normal
            obj.Fins.zCP = 0; % position du centre de pression relatif SUR LE FINS
        end
        
        function point(obj, z, m)
            % point materiel
            % Calcule la masse, le centre de masse, les moments d'inertie
            % INPUTS
            %   - z     :   position de la masse
            %   - m     :   masse (i.e. equilibrage, d?placememtn centre de
            %   gravite)
            
            % Assignation des proprietes
            obj.Point.z = z;
            obj.Point.m = m;
            
            % Calcule des proprietes de masse
            obj.Point.cm =0;%en zero car point massique
            
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