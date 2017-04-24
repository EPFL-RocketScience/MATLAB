function deriv = stateEquation(t, x, R, V0, K)
    % equation d'?tat pour l'int?gration ode45
    % INPUT :
    %   - t     : temps (s)
    %   - x     : vecteur d'int?gration
    %   - R     : objet Rocket
    %   - V0    : vent horizontal
    %   - K     : constante de correction des corps portants
    % OUTPUT :
    %   - deriv : d?riv?e de x
    
    
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
    %   - T     : temp?rature
    %   - a     : vitesse du son
    %   - p     : pression statique
    %   - rho   : densit?
    [~, a, ~, rho] = stdAtmos(Z);
    % vent relatif (incident) 
    Vi = [V0 - Vx; -Vz];
    
    % Nombre de Mach
    M = sqrt(Vx^2+Vz^2)/a;
    
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
        alpha = 0;
    else
        cr = cross([a;0],[-Vi;0]);
        alpha = atan2(cr(3),dot(a,-Vi));
    end
    % coefficient de force normale et position du centre de pouss?e
    [CNa, CP] = R.aeroCoeff(alpha, M, theta, K);
    % force normale
    N = 0.5*rho*norm(Vi)^2*Aref*CNa*alpha;
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
    Force_D = [sin(alpha); -cos(alpha)]*D;
    Force_G = [sin(phi); -cos(phi)]*m*g;
    
    % Equation d'?tat
    
    X_dot   = Vx;
    Z_dot   = Vz;
    V_dot   = (Q*(Force_T+Force_N+Force_D+Force_G) - m_dot*[Vx; Vz])/m;
    Vx_dot  = V_dot(1);
    Vz_dot  = V_dot(2);
    phi_dot = phi_dot;
    phi_ddot= -((N+D*sin(alpha))*(CP-CM)+sin(epsilon)*Ft*(L-CM)+Ir_dot*phi_dot)/Ir;
    
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

