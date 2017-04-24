function deriv = stateEquation_phi(t, x, R, V0, K)
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
    phi = x(1);     % angle de rotation par rapport ? la verticale
    phi_dot = x(2); % d?riv?e de l'angle de rotation
    
    % d?finition des constantes
    
    % d?rivation discr?te
    % ?tape de temps de d?rivation
    dt = 0.1;
        
    % environnement 
    
     M = norm(V0)/340;
     rho = 1.2;
    
    % Propri?t?s g?om?triques
    % diam?tre de r?f?rence
    d = R.d;
    % Aire de r?f?rence
    Aref = pi*d^2/4;
    
    % Propri?t?s de masse
    bt = max([R.Motor.bt]);
    Ir = R.Ir(bt);
    CM = R.cm(bt);
    
    % Propri?t?s aerodynamiques
    % angle de rotation
    theta = 0;
    % angle d'attaque
    a = [sin(phi); cos(phi); 0]; % vecteur axial unitaire
    if(norm(V0) == 0)
        alpha = 0;
    else
        alpha = atan2(cross(a,V0),dot(a,V0));
        alpha = alpha(3);
    end
    % coefficient de force normale et position du centre de pouss?e
    [CNa, CP] = R.aeroCoeff(alpha, M, theta, K);
    % force normale
    N = 0.5*rho*norm(V0)^2*Aref*CNa*alpha;
    % force de train?e
    CD = 0.85; % TODO: define CD
    D = 0.5*rho*norm(V0)^2*Aref*CD;
    
    % Equation d'?tat
    phi_ddot= -((N+D*sin(alpha))*(CP-CM))/Ir;
    
    deriv = [phi_dot phi_ddot]';
end


