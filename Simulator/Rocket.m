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
        
        d = 0; % dim?tre de r?f?rence ? la base du c?ne
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
            obj.d = D;
            obj.Nose.e = e;
            obj.Nose.rho = rho;
            
            % Calcule intermediaire
            V = e*pi*sqrt(1+D^2/(4*L^2))*L*D/2;%volume
                  
            % Calcule des proprietes de masse
            obj.Nose.m = rho*V;
            obj.Nose.cm = L*2/3;
            obj.Nose.Iz = rho/8*V*D^2;
            obj.Nose.Ir = rho*V/4*(D^2/4+2*L^2);
        
            % Calcule des proprietes aerodynamiques
            obj.Nose.CNa = @(alpha) cna_noseCone(alpha); % d?riv?e du coefficient aerodynamique normal
            obj.Nose.zCP = zCP_noseCone(L); % position du centre de pression relatif au haut du cone
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
            obj.Tail.CNa = @(alpha) cna_transition(obj.d, D1, D2, alpha); % d?riv?e du coefficient aerodynamique normal
            obj.Tail.zCP = zCP_transition(L, D1, D2); % position du centre de pression relatif
            
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
            
            % Calcule des proprietes aerodynamiques
            stage.CNa = @(alpha) cna_stage(Dout, L, obj.d, alpha); % coefficient aerodynamique normal
            stage.zCP = zCP_stage(L);
            
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
            
            motor.z = z;
            motor.m = m;
            motor.ThrustCurve = thrustCurve;
            motor.bt = bt;
            
            % Est-ce qu il y a des dimensions d un moteur en parametre de
            % la fonction ?
            % Calcule des proprietes de masse
            
            motor.cm = (1/2)*L;
            motor.Iz = @(t) m(t)*D^2/8;
            motor.Ir = @(t) m(t)/12*(3*(D/2)^2+L^2)+m(t)*(L/2)^2;
            
            obj.Motor = [obj.Motor motor];
            
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
        
        function fins(obj, z, N, Ct, Cr, xt, S, r, e, rho)
            % fins
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - z     :   position du haut de la fins
            %   - N     :   nombre de fins
            %   - Ct    :   longueur de la pointe de l'aileron
            %   - Cr    :   longueur de la base de l'aileron
            %   - xt    :   port?e axiale du bord d'attaque
            %   - S     :   envergure de l'aileron
            %   - r     :   distance de l'aileron par rapport ? l'axe
            %   vertical
            %   - e     :   epaisseur de la fins
            %   - rho   :   densite 
           
            % Assignation des proprietes
            obj.Fins.z = z;
            obj.Fins.N = N;
            obj.Fins.Ct = Ct;
            obj.Fins.Cr = Cr;
            obj.Fins.xt = xt;
            obj.Fins.S = S;
            obj.Fins.r = r;
            obj.Fins.e = e;
            obj.Fins.rho = rho;
            
            %Calcules intermediaires 
            % densit? surfacique
            lamb = rho * e;            
            % Calcule du premier triangle (Surface 1)
            Iz1 = lamb*xt*S^3/36;
            Ir1 = lamb*xt^3*S/36;
            A1  = xt*S/2;        
            cmz1 = 2*xt/3;
            cmr1 = S/3;           
            % Calcule du rectangle (Surface 2)
            Iz2 = lamb*Ct*S^3/3;
            Ir2 = lamb*Ct^3*S/3;
            A2  = Ct*S;       
            cmz2 = Ct/2;
            cmr2 = S/2;           
            % Calcule du deuxi?me triangle (Surface 3)
            xe  = Cr-Ct-xt;
            Iz3 = lamb*xe*S^3/36;
            Ir3 = lamb*xe^3*S/36;
            A3  = xe*S/2;
            cmz3 = xe/3;
            cmr3 = S/3;
            
            % Calcule des proprietes de masse
            Atot    = (Cr + Ct)/2*S;
            obj.Fins.m       = lamb*Atot;
            obj.Fins.cmz     = (A1*cmz1+A2*(xt+cmz2)+A3*(xt+Ct+cmz3))/Atot;
            obj.Fins.cmr     = (A1*cmr1+A2*cmr2+A3*cmr3)/Atot;
            obj.Fins.Iz      = Iz1 + Iz2 + Iz3 + lamb*(A1*cmr1^2+A2*cmr2^2+A3*cmr3^2);
            obj.Fins.Ir      = Ir1 + Ir2 + Ir3 + lamb*(A1*cmz1^2+A2*(xt+cmz2)^2+A3*(xt+Ct+cmz3)^2);
            obj.Fins.Itheta  = obj.Fins.Iz+obj.Fins.Ir;
                   
            % Calcule des proprietes aerodynamiques
            obj.Fins.CNa = @(M, theta) cna_fins(N, r, S, obj.d, Cr, Ct, xt, M, pi*obj.d^2/4, theta); % coefficient aerodynamique normal
            obj.Fins.zCP = zCP_fins(xt, Cr, Ct); % position du centre de pression relatif SUR LE FINS
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
         
    end
    
    methods (Access = public)
        % ?ocket property methods
        
        function m = m(obj, t)
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
            m = (m_stat + sum(cellfun(@(c) c(t), {obj.Motor.m})));

        end
        
        function cm = cm(obj, t)
            % calcule du centre de mass
            
            cm_stat = obj.Nose.m*obj.Nose.cm;
            cm_stat = cm_stat + sum([obj.Stage.m].*([obj.Stage.z] + [obj.Stage.cm]));
            cm_stat = cm_stat + obj.Tail.m*(obj.Tail.z+obj.Tail.cm);
            cm_stat = cm_stat + obj.Fins.m*obj.Fins.N*(obj.Fins.z+obj.Fins.cmz);
            cm_stat = cm_stat + sum([obj.Payload.m].*([obj.Payload.z] + [obj.Payload.cm]));
            cm_stat = cm_stat + sum([obj.Parachute.m].*([obj.Parachute.z] + [obj.Parachute.cm]));
            cm_stat = cm_stat + sum([obj.Points.m].*([obj.Points.z] + [obj.Points.cm]));
            
            Motor_m = cellfun(@(c) c(t), {obj.Motor.m});
            cm = (cm_stat + sum(Motor_m.*([obj.Motor.z]+[obj.Motor.cm])))/obj.m(t);
        end
        
        function Iz = Iz(obj, t)
            % calcule du moment d'inertie axial au centre de masse
            Iz_stat = obj.Nose.Iz;
            Iz_stat = Iz_stat + sum([obj.Stage.Iz]);
            Iz_stat = Iz_stat + obj.Tail.Iz;
            
            Iz_stat = Iz_stat + obj.Fins.N*(obj.Fins.Iz +...
                obj.Fins.m*(2*obj.Fins.cmr*obj.Fins.r+ obj.Fins.r^2));
            % moment calculation
            
            Iz_stat = Iz_stat + sum([obj.Payload.Iz]);
            
            Iz = Iz_stat + sum(cellfun(@(c) c(t), {obj.Motor.Iz}));
        end
        
        function Ir = Ir(obj, t)
            % calcule du moment d'inertie perpendiculaire ? l'axe de la
            % fus?e au centre de masse
            
            CMt = obj.cm(t);
            
            Ir = obj.Nose.Ir + obj.Nose.m*(CMt.^2-CMt*obj.Nose.cm);
            Ir = Ir + sum([obj.Stage.Ir] + [obj.Stage.m].*(...
                (CMt-[obj.Stage.z]).^2-(CMt-[obj.Stage.z]).*...
                [obj.Stage.cm]));
            Ir = Ir + obj.Tail.Ir + obj.Tail.m*(...
                (CMt-obj.Tail.z)^2-(CMt-obj.Tail.z)*obj.Tail.cm);
            
            % calcule du moment d'inertie pour les ailerons
            phi = 360/obj.Fins.N; % angle entre les ailerons
            
            for i = 1:obj.Fins.N
                phi_i = phi*(i-1); % angle de rotation du sys de coordonnees
                Ir_fin = cosd(phi_i)^2*obj.Fins.Itheta+sind(phi_i)^2*obj.Fins.Ir;
                Ir = Ir + Ir_fin + obj.Fins.m*(...
                (CMt-obj.Fins.z)^2-(CMt-obj.Fins.z)*obj.Fins.cmz);
            end
            
            Ir = Ir + sum([obj.Payload.Ir] + [obj.Payload.m].*(...
                (CMt-[obj.Payload.z]).^2-(CMt-[obj.Payload.z]).*...
                [obj.Payload.cm]));
            
            Motor_m = cellfun(@(c) c(t), {obj.Motor.m});
            Motor_Ir = cellfun(@(c) c(t), {obj.Motor.Ir});
            
            Ir = Ir + sum([Motor_Ir] + Motor_m.*(...
                (CMt-[obj.Motor.z]).^2-(CMt-[obj.Motor.z]).*[obj.Motor.cm]));
        end
        
        function [CNa_tot, zCP] = aeroCoeff(obj, alpha, M, theta)
            
            CNa_Nose = obj.Nose.CNa(alpha);
            zCP = CNa_Nose*obj.Nose.zCP;
            CNa_Stage = cellfun(@(c) c(alpha), {obj.Stage.CNa});
            zCP = zCP + sum(CNa_Stage.*([obj.Stage.zCP]+[obj.Stage.z]));
            CNa_Tail = obj.Tail.CNa(alpha);
            zCP = zCP + CNa_Tail*(obj.Tail.zCP+obj.Tail.z);
            CNa_Fins = obj.Fins.CNa(M, theta);
            zCP = zCP + CNa_Fins*(obj.Fins.zCP+obj.Fins.z);
            
            CNa_tot = (CNa_Nose + sum(CNa_Stage) + CNa_Tail + CNa_Fins);
            
            zCP = zCP / CNa_tot;
            
        end
    end
end

% calculate aerodynamic properties
           
function zcp = zCP_transition(L, d1, d2)
    % zCP_transition
    % calculate position of cp on a body element (body of revolution)
    % from tip of element 
    % INPUT:
    %   - L     : length of part [m]
    %   - d1    : diameter before transition
    %   - d2    : diameter after transition
    % RETURN:
    %   - xcp   : center of pressure position along x axis [m]
    zcp = L/3*(1+(1-d1/d2)/(1-(d1/d2)^2));
end

function  zcp = zCP_noseCone(L)
    % zCP_noseCone
    % calculate cp position of a conical nose from tip of nose
    % INPUT:
    %   - L     : nose length [m]
    % RETURN:
    %   - xcp   : center of pressure position along x axis [m]
    zcp = 2/3*L;
end

function zcp = zCP_noseOgive(L)
    % zCP_noseOgive
    % calculate cp position of an ogive shaped nose from tip
    % INPUT:
    %   -L      : nose length [m]
    % RETURN:
    %   - xcp   : center of pressure position along x axis [m]
    zcp = 0.466*L;
end

function zcp = zCP_fins(Xt, Cr, Ct)
    % zCP_fins
    % calculate cp position of any type of trapezoidal fins
    % The position is the intersection of the mean aerodynamic
    % chord and the quarter cord line
    % INPUT: 
    %   - Xt    : longitudinal distance between leading edge of
    %             root chord and leading edge of tip chord. [m]
    %   - Cr    : root chord length [m]
    %   - Ct    : tip chord length [m]
    % RETURN:
    %   - zcp   : center of pressure position along z axis [m]
    zcp = Xt/3*(Cr+2*Ct)/(Cr+Ct) + 1/6*((Cr+Ct)-(Cr*Ct)/(Cr+Ct));
end    

function zcp = zCP_stage(L)
    % zCP_stage
    % calculate cp position of a stage element
    % INPUT:
    % - L   : length of stage [m]
    % RETURN:
    % - zcp : center of pressure position along z axis [m]
    %zcp = L/2;
    zcp = 0;
end

% NORMAL AERODYNAMIC COEFFICIENT DERIVATIVES (NACD)

function cna = cna_noseCone(alpha)
    % cna_noseCone
    % calculate NACD of nose with a conical or ogive shape
    % INPUTS:
    %   - alpha : rocket angle of attack
    % RETURN:
    %   - cna   : NACD [1/rad]
    cna = 2*sind(alpha)/alpha;
end

function cna = cna_stage(D, L, d, alpha)
    % cna_noseCone
    % calculate NACD of nose with a conical or ogive shape
    % INPUTS:
    %   - D     : stage diameter
    %   - L     : stage Length
    %   - d     : reference diameter (base of cone)
    %   - alpha : rocket angle of attack
    % RETURN:
    %   - cna   : NACD [1/rad]
    
    % check d
    if isnan(d) || d <=0
       error(['Reference diameter d is not defined yet or has a',...
       'value of 0. Define a nose cone with non-zero base diameter.']);
    end
    
    %cna = 1.1*(D*L)/(pi*d^2/4)*sind(alpha)^2/alpha;
    cna = 0;
end

function cna = cna_transition(d, d1, d2, alpha)
    % cna_transition
    % calculate NACD of a body element (shoulder, boattail, tube)
    % INPUTS:
    %   - d     : reference diameter (diameter at base of cone)
    %   - d1    : diameter before transition
    %   - d2    : diameter after transition
    %   - alpha : rocket angle of attack
    % RETURN:
    %   - cna   : NACD [1/rad]    

    % check d
    if isnan(d) || d <=0
       error(['Reference diameter d is not defined yet or has a',...
       'value of 0. Define a nose cone with non-zero base diameter.']);
    end

    cna = 2/d^2*(d2^2-d1^2)*sind(alpha)/alpha;    
end

function CNa = cna_fins(N, rt, S, d, Cr, Ct, xt, M, Aref, theta)
    % cna_fins
    % calculate NACD of N (3 or 4) trapezoidal fins.
    % INPUTS:
    %   - N     : number of fins
    %   - rt    : radius of body at fin attachement
    %   - S     : span length of fin
    %   - d     : reference diameter (diameter at base of cone)
    %   - Cr    : fin root chord
    %   - Ct    : fin tip chord
    %   - xt    : longitudinal distance between leading edge of
    %             root chord and leading edge of tip chord. [m]
    %   - M     : Mach number
    %   - Aref  : Reference area (Generaly nose base area)
    %   - theta : angle of attack of fin along r axis. 
    % RETURN:
    %   - cna   : NACD [1/rad] 
    
    % check d
    if isnan(d) || d <=0
       error(['Reference diameter d is not defined yet or has a',...
       'value of 0. Define a nose cone with non-zero base diameter.']);
    end

    % check for M<1
    if (M>1)
       error('Mach number must be < 1. Code cannot guarantee validity for supersonic speeds.'); 
    end

    % check for 1<=N<=8
    if (N<1 || N>8)
       error('The number of fins must be between 1 and 8.'); 
    end

   beta     = sqrt(1-M^2);
   gammac   = atand((xt+Ct/2-Cr/2)/S); % midchord angle
   Afin     = (Cr+Ct)/2*S; % fin Area
   CNa1     = 2*pi*S^2/Aref/(1+sqrt(1+(beta*S^2/(Afin*cosd(gammac))^2))); % CNa for one fin
   CNa1     = (1+rt/(rt+S))*CNa1; %corrected CNa for fin-body interference

   % Multiple fins
   if (N == 1)
       CNa = CNa1*sind(theta)^2;
   elseif (N == 2)
       CNa = CNa1*(sind(theta)^2+sind(theta+180)^2);
   elseif (N == 3)
       CNa = CNa1*1.5*(1-0.15*cosd(3*theta/2));
   elseif (N == 4)
       CNa = CNa1*2*(1-0.06*cosd(2*theta));
   elseif (N == 5)
       CNa = 2.37*CNa1;
   elseif (N == 6)
       CNa = 2.74*CNa1;
   elseif (N == 7)
       CNa = 2.99*CNa1;
   elseif (N == 8)
       CNa = 3.24*CNa1;
   end

end 