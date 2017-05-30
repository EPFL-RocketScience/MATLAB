classdef Rocket < handle
    
    properties (SetAccess = private, GetAccess = public)
        Nose
        Stage
        Tail
        Fins
        Motor
        Cylinder
        Coupler
        Parachute
        Point
        
        d = 0; % dimetre de reference a la base du cone
    end
    
    properties (SetAccess = public, GetAccess = public)
        k = 1.1; % coefficient de correction des corps de r?volution portants
        CD = 0.8; % coefficient de frottement aerodynamique
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
            obj.Nose.CNa = @(alpha) 2*(sin(alpha)/(alpha*(alpha~=0)+(alpha==0))*(alpha~=0)+(alpha==0)); % d?riv?e du coefficient aerodynamique normal
            obj.Nose.zCP = 2/3*L; % position du centre de pression relatif au haut du cone
            Aplan = L*D/2;
            Aref = pi*D^2/4;
            obj.Nose.CNaBL = @(alpha, K) K*Aplan/Aref*(abs(sin(alpha)^2/(alpha*(alpha~=0)+(alpha==0)))*(alpha~=0));
        end
        
        function tail(obj, D1, D2, L, e, z, rho)
            % INPUTS
            % D1    :   diametre du haut
            % D2    :   diametre du bas
            % L     :   longueur du tail
            % e     :   epaisseur de paroie
            % z     :   postion par rapport au haut de la fusee du haut du tail
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
            %obj.Tail.Ir = obj.Tail.Iz/2 + pi*rho*((m^2*L^5/5+1/2*m*L^4*R1+R1^2*L^3/3)-(m^2*L^5/5+1/2*m*L^4*(R1-e)+(R1-e)^2*L^3/3));
            obj.Tail.Ir = obj.Tail.Iz/2 + pi*rho*((1/2*m*L^4*e+(R1^2-(R1-e)^2)*L^3/3));            
            
            
            % Calcule des proprietes aerodynamiques
            obj.Tail.CNa = @(alpha) 2/obj.d^2*(D2^2-D1^2)*(sin(alpha)/(alpha*(alpha~=0)+(alpha==0))*(alpha~=0)+(alpha==0)); % d?riv?e du coefficient aerodynamique normal
            obj.Tail.zCP = L/3*(1+(1-D1/D2)/(1-(D1/D2)^2)); % position du centre de pression relatif
            Aplan = (D1+D2)/2*L;
            Aref = pi*obj.d^2/4;
            obj.Tail.CNaBL = @(alpha, K) K*Aplan/Aref*(abs(sin(alpha)^2/(alpha*(alpha~=0)+(alpha==0)))*(alpha~=0));
        end
        
        function stage(obj, id, z, L, Dout, e, rho)
            % stage
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du stage (string)
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
            Din = Dout - 2*e;%Si besoin
            
            % Calcule des proprietes de masse
            %probleme au niveau de l atribution aux differents etages
            stage.m = rho*pi*L*((Dout/2)^2-(Din/2)^2);
            stage.cm = L/2;
            stage.Iz = pi*rho*L*1/2*((Dout/2)^4-(Din/2)^4);
            stage.Ir = pi*rho*L*1/12*(3*((Dout/2)^4-(Din/2)^4)+L^2*((Dout/2)^2-(Din/2)^2));
            
            % Calcule des proprietes aerodynamiques
            Aplan = Dout*L;
            Aref = obj.d^2*pi/4;
            stage.CNaBL = @(alpha, K)  K*Aplan/Aref*(abs(sin(alpha)^2/(alpha*(alpha~=0)+(alpha==0)))*(alpha~=0)); % coefficient aerodynamique normal
            stage.zCPBL = L/2;
            
            obj.Stage = [obj.Stage stage];
        end
        
         function coupler(obj, id, z, L, Dout, e, rho)
            % Coupleur (piece d'assemblage sans proprietes aerodynamiques
            % Calcule la masse, le centre de masse, le centre de pression,
            % les moments d'inertie
            % INPUTS
            %   - id    :   id du coupleur (string)
            %   - z     :   position du haut de la payload
            %   - m     :   masse
            %   - L     :   longueur de l'etage
            %   - Dout  :   diametre externe
            %   - e     :   epaisseur de la coque
            %   - rho   :   densite 
            
            % Assignation des proprietes
            coupler.id = id; % 
            coupler.z = z;
            coupler.L = L;
            coupler.Dout = Dout;
            coupler.e = e;
            coupler.rho = rho;
            
            %Calcul intermediaire
            Din = Dout - 2*e;
            
            % Calcule des proprietes de masse
            coupler.m = rho*pi*L*((Dout/2)^2-(Din/2)^2);
            coupler.cm = L/2;
            coupler.Iz = pi*rho*L*1/2*((Dout/2)^4-(Din/2)^4);
            coupler.Ir = pi*rho*L*1/12*(3*((Dout/2)^4-(Din/2)^4)+L^2*((Dout/2)^2-(Din/2)^2));
                      
            obj.Coupler = [obj.Coupler coupler];
        end
        
        function motor(obj, id, z, D, L, e, m, mp, rho, thrustCurve, bt)
            % payload
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du moteur
            %   - z     :   position du haut du moteur
            %   - D     :   Diametre externe du moteur
            %   - L     :   longueur du moteur
            %   - e     :   epaisseur de la paroi 
            %   - m     :   masse initiale
            %   - mp    :   masse du combustible
            %   - rho   :   densite du combustible
            %   - thrustCurve
            %           :   Courbe de poussee (N vs. s)
            %   - bt    :   temps de combustion (s)
            
            motor.id = id;
            motor.z = z;
            motor.D = D;
            motor.L = L;
            motor.e = e;
            motor.m = @(t)(m - mp/bt*t)*(t<=bt)+(m-mp)*(t>bt);
            motor.mp = mp;
            motor.ThrustCurve = thrustCurve;
            motor.bt = bt;
            
            % Calcule des proprietes de masse
            
            rho_shell = (m-mp)/(L*pi/4*(D^2-(D-2*e)^2));
            rho_prop = rho;
            
            if strcmp(id,'motor')
                % le moteur se consume radialement
                
                motor.cm = @(t) L/2;
                Din = @(t) sqrt((D-e)^2-4*motor.m(t)/L/rho/pi);
                % Check inner diameter 
                if ~isreal(Din(0))
                    error(['Le volume du combustible depasse la place disponible. Solution: reduire la masse de combustible ou augmenter sa densite.',...
                            'Le Diametre interne depasse de ', num2str(imag(Din(0))) ' m']);
                end
                motor.Iz = @(t) pi*rho_shell*L/2*((D/2)^4-(D/2-e)^4)...
                                + pi*rho_prop*L/2*((D/2-e)^4-(Din(t))^4);
                motor.Ir = @(t) pi*rho_shell*L/12*(3*((D/2)^4-(D/2-e)^4)...
                                + L^2*((D/2)^2-(D/2-e)^2))...
                                + pi*rho_prop*L/12*(3*((D/2-e)^4-(Din(t)/2)^4)...
                                + L^2*((D/2-e)^2-(Din(t)/2)^2));
                            
            elseif strcmp(id,'tank')
                % le reservoir se vide de haut en bas. 
                
                h = @(t) 4*(motor.m(t)-m+mp)/rho/pi/(D-2*e)^2;
                cm_shell = L/2;
                cm_prop  = @(t) L - h(t)/2*(motor.m(t)-m+mp)/mp;
                motor.cm = @(t) (cm_shell*(m-mp) + cm_prop*(motor.m(t)-m+mp))/motor.m(t);
                motor.Iz = @(t) pi*rho_shell*L/2*((D/2)^4-(D/2-e)^4)...
                                + pi*rho_prop*h(t)/2*((D/2-e)^4);
                Ir_shell = pi*rho_shell*L/12*(3*((D/2)^4-(D/2-e)^4)...
                           + L^2*((D/2)^2-(D/2-e)^2));
                Ir_prop  = @(t) (motor.m(t)-m+mp)/12*(3*(D/2-e)^2+h(t)^2);
                motor.Ir = @(t) Ir_prop(t) + (motor.cm(t)-cm_prop(t))^2*(motor.m(t)-m+mp)...
                                + Ir_shell + (motor.cm(t)-cm_shell)^2*(m-mp);
                            
            else
                
                error('Motor element must either be a motor or a tank.')
                
            end
           
            
            obj.Motor = [obj.Motor motor];
            
        end
        
        function cylinder(obj, id, z, m, L, D)
            % cylinder
            % Calcule la masse, le centre de masse, les moments d'inertie
            % INPUTS
            %   - id    :   id du payload (string)
            %   - z     :   position du haut de la payload
            %   - m     :   masse
            %   - L     :   longueur de la payload
            %   - D     :   diametre de la paylaod
            
            % Assignation des proprietes
            cylinder.id = id;
            cylinder.z = z;
            cylinder.m = m;
            cylinder.L = L;
            cylinder.D = D;
            
            % Calcule des proprietes de masse
            cylinder.cm = (1/2)*L;
            cylinder.Iz = m*D^2/8;
            cylinder.Ir = m/12*(3*(D/2)^2+L^2)+m*(L/2)^2;
            
            obj.Cylinder = [obj.Cylinder cylinder];
        end
        
        function parachute(obj, id, z, L, m, D, Cd)
            % parachutee
            % Calcule la masse, le centre de masse, le centre de pression,
            % le coefficient aerodynamique, les moments d'inertie
            % INPUTS
            %   - id    :   id du parachute
            %   - z     :   emplacement du parachute
            %   - L     :   longueure du parachute pli?
            %   - m     :   masse du parachute
            %   - D     :   diametre du parachute ouvert
            %   - Cd    :   Coefficient de trainee
            
            % Assignation des proprietes
            parachute.id = id;
            parachute.z = z;
            parachute.L = L;
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
            obj.Fins.CNa = @(M, theta) cna_fins(N, r, S, obj.d, Cr, Ct, xt, M, theta); % coefficient aerodynamique normal
            obj.Fins.zCP = xt/3*(Cr+2*Ct)/(Cr+Ct) + 1/6*((Cr+Ct)-(Cr*Ct)/(Cr+Ct)); % position du centre de pression relatif SUR LE FINS
        end
        
        function point(obj, id, z, L, m)
            % point materiel
            % Calcule la masse, le centre de masse, les moments d'inertie
            % INPUTS
            %   - id    :   id de la masse (string)
            %   - z     :   position de la masse
            %   - L     :   longueure r?elle de l'?l?ment de masse
            %   - m     :   masse (i.e. equilibrage, d?placememtn centre de
            %   gravite)
            
            % Assignation des proprietes
            point.id = id;
            point.z = z;
            point.L = L;
            point.m = m;
            
            % Calcule des proprietes de masse
            point.cm =0;%en zero car point massique
            
            obj.Point = [obj.Point point];
            
        end
         
    end
    
    methods (Access = public)
        % rocket property methods
        
        function setCoeffAero(obj, k)
           % changer la valeur du coefficient aerodynamique
           obj.k = k;
        end
        
        function m = m(obj, t)
            % calcule de la masse totale

            % masse statique
            m_stat = obj.Nose.m;
            m_stat = m_stat + sum([obj.Stage.m]);
            m_stat = m_stat + obj.Tail.m;
            m_stat = m_stat + obj.Fins.N*obj.Fins.m;
            m_stat = m_stat + sum([obj.Cylinder.m]);
            m_stat = m_stat + sum([obj.Coupler.m]);
            m_stat = m_stat + sum([obj.Parachute.m]);
            m_stat = m_stat + sum([obj.Point.m]);

            % masse variable (masse du moteur)
            m = (m_stat + sum(cellfun(@(c) c(t), {obj.Motor.m})));

        end
        
        function cm = cm(obj, t)
            % calcule du centre de mass
            % INPUT:
            %   - t     : temps
            % OUTPUT:
            %   - cm    : position du centre de mass
            
            cm_stat = obj.Nose.m*obj.Nose.cm;
            cm_stat = cm_stat + sum([obj.Stage.m].*([obj.Stage.z] + [obj.Stage.cm]));
            cm_stat = cm_stat + obj.Tail.m*(obj.Tail.z+obj.Tail.cm);
            cm_stat = cm_stat + obj.Fins.m*obj.Fins.N*(obj.Fins.z+obj.Fins.cmz);
            cm_stat = cm_stat + sum([obj.Cylinder.m].*([obj.Cylinder.z] + [obj.Cylinder.cm]));
            cm_stat = cm_stat + sum([obj.Coupler.m].*([obj.Coupler.z] + [obj.Coupler.cm]));
            cm_stat = cm_stat + sum([obj.Parachute.m].*([obj.Parachute.z] + [obj.Parachute.cm]));
            cm_stat = cm_stat + sum([obj.Point.m].*([obj.Point.z] + [obj.Point.cm]));
            
            Motor_m = cellfun(@(c) c(t), {obj.Motor.m});
            Motor_cm = cellfun(@(c) c(t), {obj.Motor.cm});
            cm = (cm_stat + sum(Motor_m.*([obj.Motor.z]+Motor_cm)))/obj.m(t);
        end
        
        function Iz = Iz(obj, t)
            % calcule du moment d'inertie axial au centre de masse
            % INPUT:
            %   - t     : temps
            % OUTPUT:
            %   - Iz    : moment d'inertie axial au centre de masse
            
            Iz_stat = obj.Nose.Iz;
            Iz_stat = Iz_stat + sum([obj.Stage.Iz]);
            Iz_stat = Iz_stat + obj.Tail.Iz;
            
            Iz_stat = Iz_stat + obj.Fins.N*(obj.Fins.Iz +...
                obj.Fins.m*(2*obj.Fins.cmr*obj.Fins.r+ obj.Fins.r^2));
            % moment calculation
            Iz_stat = Iz_stat + sum([obj.Cylinder.Iz]);
            Iz_stat = Iz_stat + sum([obj.Coupler.Iz]);
            Iz = Iz_stat + sum(cellfun(@(c) c(t), {obj.Motor.Iz}));
        end
        
        function Ir = Ir(obj, t)
            % calcule du moment d'inertie perpendiculaire a l'axe de la
            % fusee au centre de masse
            % INPUT:
            %   - t     : temps
            % OUTPUT:
            %   - Ir    : moment d'inertie perpendiculaire a l'axe, au centre de masse
            
            CMt = obj.cm(t);
            
            Ir = obj.Nose.Ir + obj.Nose.m*(CMt.^2-CMt*obj.Nose.cm);
            Ir = Ir + sum([obj.Stage.Ir] + [obj.Stage.m].*...
                    (CMt-[obj.Stage.z]-[obj.Stage.cm]).^2);
            Ir = Ir + obj.Tail.Ir + obj.Tail.m*(...
                (CMt-obj.Tail.z)^2-(CMt-obj.Tail.z)*obj.Tail.cm);
            
            % calcule du moment d'inertie pour les ailerons
            phi = 2*pi/obj.Fins.N; % angle entre les ailerons   
            for i = 1:obj.Fins.N
                phi_i = phi*(i-1); % angle de rotation du sys de coordonnees
                Ir_fin = cos(phi_i)^2*obj.Fins.Itheta+sin(phi_i)^2*obj.Fins.Ir;
                Ir = Ir + Ir_fin + obj.Fins.m*(...
                (CMt-obj.Fins.z)^2-(CMt-obj.Fins.z)*obj.Fins.cmz);
            end
            
            Ir = Ir + sum([obj.Cylinder.Ir] + [obj.Cylinder.m].*(...
                (CMt-[obj.Cylinder.z]).^2-(CMt-[obj.Cylinder.z]).*...
                [obj.Cylinder.cm]));
            
            Ir = Ir + sum([obj.Coupler.Ir] + [obj.Coupler.m].*(...
                (CMt-[obj.Coupler.z]).^2-(CMt-[obj.Coupler.z]).*...
                [obj.Coupler.cm]));
            
            Ir = Ir + sum([obj.Parachute.m].*(CMt-[obj.Parachute.z]).^2);
            
            Ir = Ir + sum([obj.Point.m].*(CMt-[obj.Point.z]).^2);
            
            Motor_m = cellfun(@(c) c(t), {obj.Motor.m});
            Motor_Ir = cellfun(@(c) c(t), {obj.Motor.Ir});
            Motor_cm = cellfun(@(c) c(t), {obj.Motor.cm});
            
            Ir = Ir + sum(Motor_Ir + Motor_m.*...
                (CMt-[obj.Motor.z]-Motor_cm).^2);
        end
        
        function [CNa_tot, zCP] = aeroCoeff(obj, alpha, M, theta)
            % calcule des proprietes aerodynamiques de la fusee
            % INPUT:
            %   - alpha : angle d'attaque (rad)
            %   - M     : nombre de Mach
            %   - theta : angle de rotation des ailerons
            
            CNa_Nose = obj.Nose.CNa(alpha);
            zCP_Nose = obj.Nose.zCP;
            CNa_NoseBL = obj.Nose.CNaBL(alpha, obj.k);
            
            CNa_Tail = obj.Tail.CNa(alpha);
            zCP_Tail = obj.Tail.zCP+obj.Tail.z;
            CNa_TailBL = obj.Tail.CNaBL(alpha, obj.k);
            
            CNa_StageBL = cellfun(@(c) c(alpha, obj.k), {obj.Stage.CNaBL});
            zCP_StageBL = [obj.Stage.zCPBL] + [obj.Stage.z];
            
            CNa_Fins = obj.Fins.CNa(M, theta);
            zCP_Fins = obj.Fins.zCP+obj.Fins.z;
            
            CNa_tot =   (CNa_Nose + CNa_NoseBL + CNa_Tail + CNa_TailBL + ...
                         sum(CNa_StageBL) + CNa_Fins);
            
            zCP =   (zCP_Nose*(CNa_Nose+CNa_NoseBL)+zCP_Tail*(CNa_Tail+CNa_TailBL)+...
                     sum(zCP_StageBL.*CNa_StageBL)+zCP_Fins*CNa_Fins)/CNa_tot;          
        end
        
        function check(obj)
            % Controle de la definition de la fusee
            
            % check : toutes les parties sont definies?
        
            if(isempty(obj.Nose))
                error('Must define a nose cone.');
            elseif(isempty(obj.Stage))
                warning('No stage defined.');
            elseif(isempty(obj.Tail))
                warning('No tail defined.');
            elseif(isempty(obj.Fins))
                warning('No fins defined.');
            elseif(isempty(obj.Motor))
                error('Must define a motor.');            
            elseif(isempty(obj.Cylinder))
                warning('No Cylinder defined.');
            elseif(isempty(obj.Coupler))
                warning('No Coupler defined.');
            elseif(isempty(obj.Parachute))
                warning('No parachute defined.');
            elseif(isempty(obj.Point))
                warning('No point mass defined.');            
            end
            
            % check : diametre de reference
            if isnan(obj.d) || obj.d <=0
                 error(['Reference diameter d is not defined yet or has a',...
                        'value of 0. Define a nose cone with non-zero base',...
                        'diameter.']);
            end
        end
        
        function printSpecs(obj)
           % affiche les specs de la fusee
           
           display('*********************');
           display('Rocket Specifications');
           display('*********************');
                                 
           display(['* L    : ' num2str(obj.getLength()) ' [m]']);
           
           % reference diameter
           D = obj.Nose.D;
           display(['* Dref : ' num2str(D) ' [m]']);
           
           % liftoff weight
           m0 = obj.m(0);
           display(['* m0   : ' num2str(m0) ' [kg]']);
           
           % propellant weights
           display('* propellant masses :');
           for motor = obj.Motor
               display(['** ' motor.id ' : ' num2str(motor.mp) ' [kg]']);
           end
           display('* total propellant mass');
           display(['** mp   : ' num2str(sum([obj.Motor.mp])) ' [kg]']);
           
           % static margin
           display('* stabililty margin');
           [~, zcp] = obj.aeroCoeff(0, 0, 0);
           display(['* ' num2str((zcp-obj.cm(0))/obj.d) ' calibers'])
        end
        
        function m = getAll_m(obj, t)
            % retourne toutes les masses dans un ordre arbitraire, mais le
            % meme que getAll_z et getAll_l
            Motor_m = cellfun(@(c) c(t), {obj.Motor.m}); 
            m = [obj.Nose.m, [obj.Stage.m], obj.Tail.m, obj.Fins.m,...
                Motor_m, [obj.Cylinder.m], [obj.Coupler.m], [obj.Parachute.m], [obj.Point.m]];
        end
        
        function z = getAll_z(obj)
            % retourne toutes les positions dans un ordre arbitraire, mais le
            % meme que getAll_m et getAll_l
            z = [0, [obj.Stage.z], obj.Tail.z, obj.Fins.z,...
                [obj.Motor.z], [obj.Cylinder.z], [obj.Coupler.z], [obj.Parachute.z], [obj.Point.z]];
        end
        
        function l = getAll_l(obj)
            % retourne toutes les positions dans un ordre arbitraire, mais le
            % meme que getAll_m et getAll_z
            l = [obj.Nose.L, [obj.Stage.L], obj.Tail.L, obj.Fins.Cr,...
                [obj.Motor.L], [obj.Cylinder.L], [obj.Coupler.L], [obj.Parachute.L], [obj.Point.L]];
        end
        
        function L = getLength(obj)
            % retourne la longueure de la fusee
           L = max(obj.getAll_z()+obj.getAll_l()); 
        end
        
    end
end

function CNa = cna_fins(N, rt, S, d, Cr, Ct, xt, M, theta)
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
       warning('Mach number must be < 1. Code cannot guarantee validity for supersonic speeds.'); 
    end

    % check for 1<=N<=8
    if (N<1 || N>8)
       error('The number of fins must be between 1 and 8.'); 
    end
   Aref     = pi*d^2/4;
   beta     = sqrt(1-M^2);
   gammac   = atan((xt+Ct/2-Cr/2)/S); % midchord angle
   Afin     = (Cr+Ct)/2*S; % fin Area
   CNa1     = 2*pi*S^2/Aref/(1+sqrt(1+(beta*S^2/(Afin*cos(gammac)))^2)); % CNa for one fin
   CNa1     = (1+rt/(rt+S))*CNa1; %corrected CNa for fin-body interference

   % Multiple fins
   if (N == 1)
       CNa = CNa1*sin(theta)^2;
   elseif (N == 2)
       CNa = CNa1*(sin(theta)^2+sin(theta+pi)^2);
   elseif (N == 3)
       CNa = CNa1*1.5*(1-0.15*cos(3*theta/2));
   elseif (N == 4)
       CNa = CNa1*2*(1-0.06*cos(2*theta));
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