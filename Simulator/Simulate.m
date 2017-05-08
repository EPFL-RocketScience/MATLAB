function [tsim, Xsim, alpha, calibre, T, M] = Simulate( R, V0, tfin, phi0, lramp, tquer, xquer)
%SIMULATE effectue la simulation de la fusee
%   INPUTS :
%       - R     : objet 'Rocket'
%       - V0    : vitesse du vent [m/s]
%       - K     : coefficient de correction pour la portance des corps
%       - tfin  : temps de simulation maximal
%       - phi0  : angle de d?part de la rampe en [rad]
%       - lramp : longueure de la rampe de lancement [m]
%       - tquer : times at which special calculations should be done (flexion)
%       - xquer : positions along rocket where special calculation values
%                 are requested.
%   OUTPUTS :
%       - tsim  : etapes de temps d'integration
%       - Xsim  : etat du system a chaque etape
%       - alpha : angle d'attaque
%       - calibre 
%               : calibres entre le CM et le CP 
%       - T     : shear values at query points and query times
%       - M     : flexion values at query points and query times
K=R.k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global simulation variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gravite
    g = 9.8;
    % acceleration
    Vx_dot = 0;
    Vz_dot = 0;
    % angle of attack
    alpha_tmp = 0;
    alpha = [];
    % calibres
    calibre_temp = 0;
    calibre = [];
    % remember last query time
    tquer_i = 1;
    % initialize shear and flexion matrices
    T = zeros(length(tquer), length(xquer));
    M = zeros(length(tquer), length(xquer));
    
    %%%%%%%%%%%%%%%%%%%%%%%%         
    % find start time
    %%%%%%%%%%%%%%%%%%%%%%%%
    % when motor thrust compensates gravity
    FT_i = min(find(R.Motor.ThrustCurve(:,2)>cos(phi0)*R.m(0)*g));
    t_start = R.Motor.ThrustCurve(FT_i,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%         
    % integration options
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % integrator options launch
    options_launch = odeset('Events',@events_launch,'OutputFcn',@odeplot, 'OutputFcn', @output,...
                     'Refine', 1);
    tspan_launch = [t_start, 10];
    
    % integrator options flight
    options_flight = odeset('Events', @events_flight, 'OutputFcn', @output,...
                     'Refine', 1);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%         
    % integration
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    display('Started Sim...')
    
    % launch
    x0_launch = [0, 0, 0, 0, phi0, 0];
    [tsim_launch, Xsim_launch] = ode45(@(t, x) stateEquation(t, x, R, 0, R.k, 0), tspan_launch, x0_launch, options_launch);
    display('Cleared launch pad:');
    display(['t = ' num2str(tsim_launch(end)) ' sec']);
    display(['Z = ' num2str(Xsim_launch(end, 2)) ' m']);
    display(['V = ' num2str(Xsim_launch(end, 4)) ' m/s']);
    
    % propulsed flight
    tspan_flight_withFt = [tsim_launch(end), R.Motor.ThrustCurve(end,1)];
    [tsim_flight_withFt, Xsim_flight_withFt] = ode45(@(t, x) stateEquation(t, x, R, V0, K, 1), tspan_flight_withFt, Xsim_launch(end, :), options_flight);
    display('Motor burnout:');
    display(['t = ' num2str(tsim_flight_withFt(end)) ' sec']);
    display(['Z = ' num2str(Xsim_flight_withFt(end, 2)) ' m']);
    display(['V = ' num2str(Xsim_flight_withFt(end, 4)) ' m/s']);
    
    % coasting
    tspan_flight_withOutFt = [R.Motor.ThrustCurve(end,1), tfin];
    [tsim_flight_withOutFt, Xsim_flight_withOutFt] = ode45(@(t, x) stateEquation(t, x, R, V0, K, 2), tspan_flight_withOutFt, Xsim_flight_withFt(end, :), options_flight);
    display('Apogee:');
    display(['t = ' num2str(tsim_flight_withOutFt(end)) ' sec']);
    display(['Z = ' num2str(Xsim_flight_withOutFt(end, 2)) ' m']);
    display(['V = ' num2str(Xsim_flight_withOutFt(end, 4)) ' m/s']);
    
    
    %%%%%%%%%%%%%%%
    % output
    %%%%%%%%%%%%%%%
    
    tsim = [tsim_launch; tsim_flight_withFt; tsim_flight_withOutFt];
    Xsim = [Xsim_launch; Xsim_flight_withFt; Xsim_flight_withOutFt];
    
    
    %%%%%%%%%%%%%%%%%%
    % equation d'etat
    %%%%%%%%%%%%%%%%%%
    
    function deriv = stateEquation(t, x, R, V0, K, phase)
    % deriv = stateEquation(t, x, R, V0, K) Equation d'etat pour
    % l'integration


        % definition des variables a int?grer
        X = x(1);       % position horizontale dans le repere terrestre
        Z = x(2);       % position verticale dans le rep?re terrestre
        Vx = x(3);
        Vz = x(4);
        phi = x(5);     % angle de rotation par rapport a la verticale
        phi_dot = x(6); % derive de l'angle de rotation

        % definition des constantes

        % matrice de rotation: fusee -> terrestre
        Q = [cos(phi), sin(phi); -sin(phi), cos(phi)];

        % derivation discrete
        % etape de temps de derivation
        dt = 0.1;

        % environnement 
        % etat de l'air
        %   - a     : vitesse du son
        %   - rho   : densit?
        [~, a, ~, rho] = stdAtmos(Z);
        % vent relatif (incident) 
        Vi = [V0 - Vx; -Vz];

        % Nombre de Mach
        Mach = sqrt(Vx^2+Vz^2)/a;

        % Proprietes geometriques
        % diametre de reference
        d = R.d;
        % Aire de reference
        Aref = pi*d^2/4;
        L = R.Tail.z+R.Tail.L;

        % Proprietes de masse
        bt = max([R.Motor.bt]);
        m = R.m(t);
        if(t==0)
            m_dot = (R.m(dt)-R.m(0))/dt;
        elseif(t<bt)
            m_dot = (R.m(t+dt)-R.m(t-dt))/(2*dt);
        else
            m_dot = 0;
        end
        Ir = R.Ir(t);
        if(t==0)
            Ir_dot = (R.Ir(dt)-R.Ir(0))/dt;
        elseif(t<bt)
            Ir_dot = (R.Ir(t+dt)-R.Ir(t-dt))/(2*dt);
        else
            Ir_dot = 0;
        end
        CM = R.cm(t);

        % Proprietes aerodynamiques
        % angle de rotation
        theta = 0;
        % angle d'attaque
        a = [sin(phi); cos(phi)];
        if(norm(Vi) == 0)
            alpha_tmp = 0;
        else
            cr = cross([a;0],[-Vi;0]);
            alpha_tmp = atan2(cr(3),dot(a,-Vi));
        end
        % coefficient de force normale et position du centre de pouss?e
        [CNa, CP] = R.aeroCoeff(alpha_tmp, Mach, theta);
        % force normale
        N = 0.5*rho*norm(Vi)^2*Aref*CNa*alpha_tmp;
        % force de train?e
        CD = 0.85; % TODO: define CD
        D = 0.5*rho*norm(Vi)^2*Aref*CD;

        % Proprietes du moteur
        if(t<R.Motor.ThrustCurve(1,1))
            Ft = 0;
        elseif(t<R.Motor.ThrustCurve(end,1));
            Ft = interp1(R.Motor.ThrustCurve(:,1)', R.Motor.ThrustCurve(:,2)', t);
        else
            Ft = 0;
        end
        epsilon = 0;
        
        %Calcul du calibre
        calibre_temp = (CP-CM)/d;

        %Calcul du calibre
        calibre_temp = (CP-CM)/d;
        
        % Forces dans le repere (n, a)
        % Poussee
        Force_T = [sin(epsilon); cos(epsilon)]*Ft;
        Force_N = [N;0];
        Force_D = [sin(alpha_tmp); -cos(alpha_tmp)]*D;
        Force_G = [sin(phi); -cos(phi)]*m*g;

        % Equation d'etat
        
        % launch
        if phase == 0
            % Force du rail sur la fusee
            Force_R = [sin(phi); 0]*m*g;
            
            X_dot   = Vx;
            Z_dot   = Vz;
            V_dot   = (Q*(Force_T+Force_N+Force_D+Force_G-Force_R) - m_dot*[Vx; Vz])/m;
            Vx_dot  = V_dot(1);
            Vz_dot  = V_dot(2);
            phi_dot = 0;
            phi_ddot= 0;
        
        % propulsed flight    
        elseif phase == 1
        
            X_dot   = Vx;
            Z_dot   = Vz;
            V_dot   = (Q*(Force_T+Force_N+Force_D+Force_G) - m_dot*[Vx; Vz])/m;
            Vx_dot  = V_dot(1);
            Vz_dot  = V_dot(2);
            phi_dot = phi_dot;
            phi_ddot= -((N+D*sin(alpha_tmp))*(CP-CM)+sin(epsilon)*Ft*(L-CM)+Ir_dot*phi_dot)/Ir;

              
        % coast    
        elseif phase == 2
        
            X_dot   = Vx;
            Z_dot   = Vz;
            V_dot   = (Q*(Force_N+Force_D+Force_G) - m_dot*[Vx; Vz])/m;
            Vx_dot  = V_dot(1);
            Vz_dot  = V_dot(2);
            phi_dot = phi_dot;
            phi_ddot= -((N+D*sin(alpha_tmp))*(CP-CM)+sin(epsilon)*Ft*(L-CM)+Ir_dot*phi_dot)/Ir;

        end

        deriv = [X_dot Z_dot Vx_dot Vz_dot phi_dot phi_ddot]';
    end
    
    function status = output(t,x,flag)
    % status = output(t, x, flag) function is called each time integrator
    % does a valid step.
    
        % as long as status = 0, ode45 will keep on runing... can be used
        % to create a stop condition.
        status = 0;

        % at initialization, and after each succesful step
        if isempty(flag) || strcmp(flag, 'init')
            
            % calculate angle of attack
            alpha(end+1) = alpha_tmp;
            calibre(end+1) = calibre_temp;
            
            % calculate flexion
            if(tquer_i<=length(tquer) & t>=tquer(tquer_i))
                Q = [cos(x(5)), sin(x(5)); -sin(x(5)), cos(x(5))];
                a = Q*[Vx_dot; Vz_dot+9.8];
                m = R.getAll_m(t(end));
                z = R.getAll_z();
                l = R.getAll_l();
                [T(tquer_i,:), M(tquer_i,:)] = Flexion(m', l', z', xquer, a(1), xquer(end));
                
                % increment query time index
                tquer_i = tquer_i + 1;
            end
            
        end

    end 

    function [value,isterminal,direction] = events_launch(tsim, Xsim)
        %fonction qui gere l'evenement d'arrivee en fin de rampe
        
        verif_x = Xsim(1)/sin(phi0) - lramp;
        verif_z = Xsim(2)/cos(phi0) - lramp;
        value = [verif_x, verif_z, 0, 0, 0, 0]; % On verifie pour x et z
        isterminal = [1, 1, 0, 0, 0, 0]; % stop the integration
        direction = [0, 0, 0, 0, 0, 0]; % negative direction
        
    end

    function [value, isterminal, direction] = events_flight(tsim, Xsim)
        % fonction qui detecte l'apogee
        
        Vz = Xsim(4);
        value = [0, 0, 0, Vz, 0, 0];
        isterminal = [0, 0, 0, 1, 0, 0];
        direction = [0, 0, 0, -1, 0, 0];
        
    end
end

