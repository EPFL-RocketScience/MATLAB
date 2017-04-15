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
        
        Mass
        CM
        Iz
    end
    
    methods
        function obj =  Rocket()
            % Rocket
            % Constructor of Rocket, empty for now
        end
    end
    
    methods
        function nose(obj, L, D, e, rho)
            % nose
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - L     :   longueur du cone
            %   - D     :   Diam?tre de la base du cone
            %   - e     :   epaisseur du cone
            %   - rho   :   densite
            %   - type  :   cone ou ogive
            
            % Assignation des proprietes
            obj.Nose.L = L;
            obj.Nose.D = D;
            obj.Nose.e = e;
            obj.Nose.rho = rho;
            
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
            % D1    :   diametre du haut
            % D2    :   diametre du bas
            % L     :   longueur du tail
            % e     :   epaisseur de paroie
            % z     :   postion par rapport au haut de la fus?e du haut du tail
            % rho   :   densite
            
            % Assignation des proprietes
            obj.Tail.D1 = D1;
            obj.Tail.D2 = D2;
            obj.Tail.L = L;
            obj.Tail.e = e;
            obj.Tail.z = z;
            obj.Tail.rho = rho;
            
            % Calcule des proprietes de masse
            %on utilise un cone - un plus petit cone
            V = e*pi*sqrt(1+(D2-D1)^2/4/L^2)*L*(D1+D2)/2;
            obj.Tail.m = rho*V;
            obj.Tail.cm = L*(2*D2+D1)/(3*(D2+D1));            
            
            R1 = D1/2;
            R2 = D2/2;
            m = (R2-R1)/L;
            obj.Tail.Iz = pi*rho/(10*m)*((R2^5-R1^5)-((R2-e)^5-(R1-e)^5));
            obj.Tail.Ir = obj.Tail.Iz/2 + pi*rho*((m^2*L^5/5+1/2*m*L^4*R1+R1^2*L^3/3)-(m^2*L^5/5+1/2*m*L^4*(R1-e)+(R1-e)^2*L^3/3));
            
            % Calcule des proprietes aerodynamiques
            obj.Tail.CN = 0; % coefficient aerodynamique normal
            obj.Tail.zCP = 0; % position du centre de pression relatif
            
        end
        
        function stage(obj, id, z, L, Dout, e, rho)
            % stage
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du payload (string)
            %   - z     :   position du haut de la payload
            %   - m     :   masse
            %   - L     :   longueur de l'etage
            %   - Dout  :   diametre externe
            %   - e     :   epaisseur de la coque
            %   - rho   :   densite 
            
            % Assignation des proprietes
            stage.id = id; % 
            stage.z = z;
            stage.L = L;
            stage.Dout = Dout;
            stage.e = e;
            stage.rho = rho;
            
            %Calcul intermediaire
            Din = Dout - e;%Si besoin
            
            % Calcule des proprietes de masse
            %probleme au niveau de l atribution aux differents etages
            stage.m = rho*pi*L*((Dout/2)^2-(Din/2)^2);
            stage.cm = L/2;
            stage.Iz = pi*rho*L*1/2*((Dout/2)^4-(Din/2)^4);
            stage.Ir = pi*rho*L*1/12*(3*((Dout/2)^4-(Din/2)^4)+L^2*((Dout/2)^2-(Din/2)^2));
            
            % Calcule des proprietes aerodynamiques ??? pas de propr. non ?
            stage.CN = 0; % coefficient aerodynamique normal
            stage.zCP = 0;
            
            obj.Stage = [obj.Stage stage];
        end
        
        function motor(obj, z, m, D, L, thrustCurve, bt)
            % payload
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   position du haut du moteur
            %   - m     :   masse
            %   - D     :   Diametre externe du moteur
            %   - L     :   longueur du moteur
            %   - thrustCurve
            %           :   Courbe de poussee (N vs. s)
            %   - bt    :   temps de combustion (s)
            
            if(~isa(m, 'function_handle'))
                error('La masse doit etre une fonction')
            end
            
            obj.Motor.z = z;
            obj.Motor.m = m;
            obj.Motor.ThrustCurve = thrustCurve;
            obj.Motor.bt = bt;
            
            % Est-ce qu il y a des dimensions d un moteur en parametre de
            % la fonction ?
            % Calcule des proprietes de masse
            
            obj.Motor.cm = (1/2)*L;
            obj.Motor.Iz = @(t) m(t)*D^2/8;
            obj.Motor.Ir = @(t) m(t)/12*(3*(D/2)^2+L^2)+m(t)*(L/2)^2;
            
        end
        
        function payload(obj, id, z, m, L, D)
            % payload
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du payload (string)
            %   - z     :   position du haut de la payload
            %   - m     :   masse
            %   - L     :   longueur de la payload
            %   - D     :   diametre de la paylaod
            
            % Assignation des proprietes
            payload.id = id;
            payload.z = z;
            payload.m = m;
            payload.L = L;
            payload.D = D;
            
            % Calcule des proprietes de masse
            payload.cm = (1/2)*L;
            payload.Iz = m*D^2/8;
            payload.Ir = m/12*(3*(D/2)^2+L^2)+m*(L/2)^2;
            
            obj.Payload = [obj.Payload payload];
        end
        
        function parachute(obj, id, z, m, D, Cd)
            % parachutee
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du parachute
            %   - z     :   emplacement du parachute
            %   - m     :   masse du parachute
            %   - D     :   diametre du parachute ouvert
            %   - Cd    :   Coefficient de trainee
            
            % Assignation des proprietes
            parachute.id = id;
            parachute.z = z;
            parachute.m = m;
            parachute.D = D;
            parachute.Cd = Cd;
            
            % Calcule des proprietes de masse
            parachute.cm = 0;%car masse ponctuelle
            
            obj.Parachute = [obj.Parachute parachute];
            
        end
        
        function fins(obj, z, N, Ct, Cr, gamma, phi, r, e, rho)
            % fins
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   position du haut de la fins
            %   - N     :   nombre de fins
            %   - Ct    :   longueur de la pointe de l'aileron
            %   - Cr    :   longueur de la base de l'aileron
            %   - gamma :   angle ? l'avant de la fins
            %   - phi   :   angle ? l'arriere de la fins
            %   - r     :   distance de l'aileron par rapport ? l'axe
            %   vertical
            %   - e     :   epaisseur de la fins
            %   - rho   :   densite 
           
            % Assignation des proprietes
            obj.Fins.z = z;
            obj.Fins.N = N;
            obj.Fins.Ct = Ct;
            obj.Fins.Cr = Cr;
            obj.Fins.gamma = gamma;
            obj.Fins.phi = phi;
            obj.Fins.r = r;
            obj.Fins.e = e;
            obj.Fins.rho = rho;
            
            %Calcule intermediaire 
            obj.Fins.h = (Cr-Ct)/(1/tan(gamma)+1/tan(phi)); 
            
            % Calcule des proprietes de masse
            obj.Fins.cmz = 0; % coordonn?e du centre de masse dans l'axe de la fus?e
            obj.Fins.cmr = 0; % coordonn?e du centre de masse dans le plan transverse ? l'axe de la fusee
            obj.Fins.m = 0;
            obj.Fins.Iz = 0;
            obj.Fins.Ir = 0;
                   
            % Calcule des proprietes aerodynamiques
            obj.Fins.CN = 0; % coefficient aerodynamique normal
            obj.Fins.zCP = 0; % position du centre de pression relatif SUR LE FINS
        end
        
        function point(obj, id, z, m)
            % point materiel
            % Calcule la masse, le centre de masse, les moments d'inertie
            % INPUTS
            %   - id    :   id de la masse (string)
            %   - z     :   position de la masse
            %   - m     :   masse (i.e. equilibrage, d?placememtn centre de
            %   gravite)
            
            % Assignation des proprietes
            point.id = id;
            point.z = z;
            point.m = m;
            
            % Calcule des proprietes de masse
            point.cm =0;%en zero car point massique
            
            obj.Points = [obj.Points point];
            
        end
    
        function update(obj)
            obj.Mass    = obj.calc_m;
            obj.CM      = obj.calc_cm; 
            obj.Iz      = obj.calc_Iz;
        end
        
    end
    
    methods (Access = private)
        
        function m = calc_m(obj)
            % calcule de la masse totale

            % masse statique
            m_stat = obj.Nose.m;
            m_stat = m_stat + sum([obj.Stage.m]);
            m_stat = m_stat + obj.Tail.m;
            m_stat = m_stat + obj.Fins.m;
            m_stat = m_stat + sum([obj.Payload.m]);
            m_stat = m_stat + sum([obj.Parachute.m]);
            m_stat = m_stat + sum([obj.Points.m]);

            % masse variable (masse du moteur)
            m = @(t) (m_stat + obj.Motor.m(t));

        end
        
        function cm =calc_cm(obj)
            % calcule du centre de mass
            
            cm_stat = obj.Nose.m*obj.Nose.cm;
            cm_stat = cm_stat + sum([obj.Stage.m].*([obj.Stage.z] + [obj.Stage.cm]));
            cm_stat = cm_stat + obj.Tail.m*(obj.Tail.z+obj.Tail.cm);
            cm_stat = cm_stat + obj.Fins.m*(obj.Fins.z+obj.Fins.cm);
            cm_stat = cm_stat + sum([obj.Payload.m].*([obj.Payload.z] + [obj.Payload.cm]));
            cm_stat = cm_stat + sum([obj.Parachute.m].*([obj.Parachute.z] + [obj.Parachute.cm]));
            cm_stat = cm_stat + sum([obj.Points.m].*([obj.Points.z] + [obj.Points.cm]));
            
            cm = @(t) (cm_stat + obj.Motor.m(t).*(obj.Motor.z+obj.Motor.cm))./obj.Mass(t);
        end
        
        function Iz = calc_Iz(obj)
            % calcule du moment d'inertie axial au centre de masse
            Iz_stat = obj.Nose.Iz;
            Iz_stat = Iz_stat + sum([obj.Stage.Iz]);
            Iz_stat = Iz_stat + obj.Tail.Iz;
            
            Iz_stat = Iz_stat + obj.Fins.N*(obj.Fins.Iz +...
                obj.Fins.m*(2*obj.Fins.cmr*obj.Fins.r+ obj.Fins.r^2));
            % moment calculation
            
            Iz_stat = Iz_stat + sum([obj.Payload.Iz]);
            
            Iz = @(t) (Iz_stat + obj.Motor.Iz(t));
        end
        
        function Ir = calc_Ir(obj)
            % calcule du moment d'inertie perpendiculaire ? l'axe de la
            % fus?e au centre de masse
            
            Ir_stat = obj.Nose.Ir + obj.Nose.m*(obj.CM^2-obj.CM*obj.Nose.cm);
            Ir_stat = Ir_stat + sum([obj.Stage.Ir] + [obj.Stage.m].*...
                (obj.CM-[obj.Stage.z]).^2-([obj.Stage.cm]-[obj.Stage.z]).*...
                [obj.Stage.cm])
            Ir_stat = Ir_stat + obj.Tail.Ir + obj.Tail.m*(...
                (obj.CM-obj.Tail.z)^2-(obj.CM-obj.Tail.z)*obj.Tail.cm);
            
            % calcule du moment d'inertie pour les ailerons
            phi = 360/obj.Fins.N; % angle entre les ailerons
            
            for i = 1:obj.Fins.N
                phi_i = phi*(i-1); % angle de rotation du sys de coordonnees
                Ir = cosd(phi_i)^2*obj.Fins.Itheta+sind(phi_i)^2*obj.Fins.Ir;
                Ir_stat = Ir_stat + Ir + obj.Fins.m*(...
                (obj.CM-obj.Fins.z)^2-(obj.CM-obj.Fins.z)*obj.Fins.cmz);
            end
            
            Ir_stat = Ir_stat + obj.
        end
    end
end