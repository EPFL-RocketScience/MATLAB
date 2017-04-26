function deriv = stateEquation(t, x, R, V0)
    % equation d'?tat pour l'int?gration ode45
    % INPUT :
    %   - t     : temps (s)
    %   - x     : vecteur d'int?gration
    %   - R     : objet Rocket
    %   - V0    : vent horizontal
    % OUTPUT :
    %   - deriv : d?riv?e de x
    
    
    % d?finition des variables ? int?grer
    X = x(1);       % position horizontale dans le rep?re terrestre
    Z = x(2);       % position verticale dans le rep?re terrestre
    Vr = x(3);      % vitesse du CM le long de l'axe de la fus?e  
    Vz = x(4);      % vitesse du CM perpendiculaire ? l'axe de la fus?e
    phi = x(5);     % angle de rotation par rapport ? la verticale
    phi_dot = x(6); % d?riv?e de l'angle de rotation
    
    % d?finition des constantes
    
    % d?rivation discr?te
    % ?tape de temps de d?rivation
    dt = 0.1;
        
    % environnement 
    % gravite
    g = 9.8;
    % ?tat de l'air
    %   - T     : temp?rature
    %   - a     : vitesse du son
    %   - p     : pression statique
    %   - rho   : densit?
    [T, a, p, rho] = stdAtmos(Z);
    % vent relatif (incident) 
    Vi = [V0 - cosd(phi)*Vr+sind(phi)*Vz;-sind(phi)*Vr-cosd(phi)*Vz];
    
    % Nombre de Mach
    M = sqrt(Vr^2+Vz^2)/a;
    
    % Propri?t?s g?om?triques
    % diam?tre de r?f?rence
    d = R.d;
    % Aire de r?f?rence
    Aref = pi*d^2/4;
    L = R.Tail.z+R.Tail.L;
    
    % Propri?t?s de masse
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
    
    % Propri?t?s aerodynamiques
    % angle de rotation
    theta = 0;
    % constante de correction pour les corps portants
    K = 1.1;
    % angle d'attaque
    ez = [sind(phi); cosd(phi)];
    if(norm(Vi) == 0)
        alpha = 0;
    else
        alpha = acosd(ez'*-Vi/(norm(ez)*norm(Vi)));
        %alpha = atan2(norm(cross(ez,Vi)),dot(ez,Vi));
    end
    % coefficient de force normale et position du centre de pouss?e
    [CNa, CP] = R.aeroCoeff(alpha*pi/180, M, theta, K);
    % force normale
    N = 0.5*rho*norm(Vi)^2*Aref*CNa*alpha*pi/180;
    % force de train?e
    CD = 0.8; % TODO: define CD
    D = -0.5*rho*norm(Vi)^2*Aref*CD;
    
    % Propri?t?s du moteur
    if(t<R.Motor.ThrustCurve(1,1))
        Ft = 0;
    elseif(t<R.Motor.ThrustCurve(end,1));
        Ft = interp1(R.Motor.ThrustCurve(:,1)', R.Motor.ThrustCurve(:,2)', t);
    else
        Ft = 0;
    end
    epsilon = 0;
    
    
    % Equation d'?tat
    %phi_dot = phi_dot;
    if(Z<5)
        X_dot = 0;
        Z_dot = Vz;
        Vr_dot = 0;
        Vz_dot = -m_dot/m*Vz-g+1/m*(Ft*cosd(epsilon)+D);
        phi_dot = 0;
        phi_ddot = 0;
    else
        phi_ddot= 1/Ir*((CP-CM)*(N-sind(alpha)*D)+(L-CM)*sind(epsilon)*Ft-Ir_dot*phi_dot);
        V_dot   =   [cosd(phi), -sind(phi); sind(phi), cosd(phi)]*(...
                    phi_dot*[sind(phi), cosd(phi); -cosd(phi), sind(phi)]*[Vr;Vz]...
                    -m_dot/m*[cosd(phi), -sind(phi); sind(phi), cosd(phi)]*[Vr;Vz]...
                    +1/m*[N*cosd(phi)-D*(cosd(phi)*sind(alpha)+sind(phi)*cosd(alpha))+Ft*(cosd(phi)*sind(epsilon)-sind(phi)*cosd(epsilon));...
                          -m*g+N*sind(phi)+D*(-sind(phi)*sind(alpha)+cosd(phi)*cosd(alpha))+Ft*(sind(phi)*sind(epsilon)+cosd(phi)*cosd(epsilon))]);
        Vr_dot  = V_dot(1);
        Vz_dot  = V_dot(2);
        X_dot   = cosd(phi)*Vr-sind(phi)*Vz;
        Z_dot   = sind(phi)*Vr+cosd(phi)*Vz;
    end
    
    deriv = [X_dot, Z_dot, Vr_dot, Vz_dot, phi_dot, phi_ddot]';
end

