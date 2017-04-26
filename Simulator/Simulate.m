function [tsim, Xsim, alpha, T, M] = Simulate( R, V0, K, tspan, tquer, xquer)
%SIMULATE effectue la simulation de la fusee
%   INPUTS :
%       - R     : objet 'Rocket'
%       - V0    : vitesse du vent [m/s]
%       - K     : coefficient de correction pour la portance des corps
%       - tspan : interval de temps de la simulation [t_start t_end]
%       - tquer : times at which special calculations should be done (flexion)
%       - xquer : positions along rocket where special calculation values
%                 are requested.
%   OUTPUTS :
%       - tsim  : etapes de temps d'integration
%       - Xsim  : etat du system a chaque etape
%       - alpha : angle d'attaque
%       - T     : shear values at query points and query times
%       - M     : flexion values at query points and query times

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global simulation variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % acceleration
    Vx_dot = 0;
    Vz_dot = 0;
    % angle of attack
    alpha_tmp = 0;
    alpha = [];
    % remember last query time
    tquer_i = 1;
    % initialize shear and flexion matrices
    T = zeros(length(tquer), length(xquer));
    M = zeros(length(tquer), length(xquer));
    
    % integrator options
    options = odeset('OutputFcn', @output,...
                     'Refine', 1);
    
    % integration
    [tsim, Xsim] = ode45(@(t, x) stateEquation(t, x, R, V0, K), tspan, [0, 0, 0, 0, 0, 0], options);
    
    
    function deriv = stateEquation(t, x, R, V0, K)
    % deriv = stateEquation(t, x, R, V0, K) Equation d'etat pour
    % l'integration


        % d?finition des variables ? int?grer
        X = x(1);       % position horizontale dans le rep?re terrestre
        Z = x(2);       % position verticale dans le rep?re terrestre
        Vx = x(3);
        Vz = x(4);
        phi = x(5);     % angle de rotation par rapport ? la verticale
        phi_dot = x(6); % d?riv?e de l'angle de rotation

        % d?finition des constantes

        % debug?
        debug = 0;

        % matrice de rotation: fusee -> terrestre
        Q = [cos(phi), sin(phi); -sin(phi), cos(phi)];

        % derivation discrete
        % etape de temps de derivation
        dt = 0.1;

        % environnement 
        % gravite
        g = 9.8;
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
        [CNa, CP] = R.aeroCoeff(alpha_tmp, Mach, theta, K);
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

        % Forces dans le repere (n, a)
        % Poussee
        Force_T = [sin(epsilon); cos(epsilon)]*Ft;
        Force_N = [N;0];
        Force_D = [sin(alpha_tmp); -cos(alpha_tmp)]*D;
        Force_G = [sin(phi); -cos(phi)]*m*g;

        % Equation d'?tat

        X_dot   = Vx;
        Z_dot   = Vz;
        V_dot   = (Q*(Force_T+Force_N+Force_D+Force_G) - m_dot*[Vx; Vz])/m;
        Vx_dot  = V_dot(1);
        Vz_dot  = V_dot(2);
        phi_dot = phi_dot;
        phi_ddot= -((N+D*sin(alpha_tmp))*(CP-CM)+sin(epsilon)*Ft*(L-CM)+Ir_dot*phi_dot)/Ir;

        deriv = [X_dot Z_dot Vx_dot Vz_dot phi_dot phi_ddot]';

        % debug
        if debug == 1
            feather(X_dot, Z_dot, 'b');
            v = norm([X_dot, Z_dot]);
            hold on;
            feather(a(1)*v, a(2)*v, 'r');
            hold off;
        end
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
            
            % calculate flexion
            if(t>=tquer(tquer_i))
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

end
