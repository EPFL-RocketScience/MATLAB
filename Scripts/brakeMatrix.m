function Z2 = brakeMatrix( m0, mp, D, Arefb, Cd0, Cdb, TC, BT, xt)
%BRAKEMATRIX calculate matrix containing data to select braking altitude.
%   brakeMatrix( m, Cd, motor )
%   INPUTS : 
%       - m0    :   total mass of rocket at liftoff (kg)
%       - mp    :   total propellant mass at liftoff (kg)
%       - D     :   Reference diameter, generally at base of cone (m)
%       - Arefb :   Brake reference area (m^2)
%       - Cd0   :   estimated constant coefficient of drag (no brakes)
%       - Cdb   :   estimated constant coefficient of drag (for brakes only)
%       - TC    :   Thrust curve vector (N)
%       - BT    :   Motor burn time (sec)
%       - xt    :   target altitude (m)
%   OUTPUT :
%       - Z2    :   matrix of query points for microcontroller

    %% 0. constants

    % 0.1 reference Area
    Aref0 = pi*D^2/4;
    
    % 0.2 motor uncertainty factor
    eta = 0.2;
    
    % 0.3 grid refinement, number of grid points
    refine = 100; 

    %% 1. determine plausible burnout conditions
    
    % 1.1 Mass: assume linear mass variation
    m       = @(t) (m0 - mp/BT * t) * (t<=BT) + (m0-mp) * (t>BT);
    m_dot   = @(t) - mp/BT * (t<=BT);
    
    % 1.2 Thrust curve : interpolate points
    Ft      = @(t) interp1(TC(:,1)', TC(:,2)', t);
    
    % 1.3 Simulation time : from first motor data point to last motor data
    % point
    tSpan   = [min(TC(:,1)), max(TC(:,1))]; 
    
    % 1.4 Initial conditions : no height / no speed
    x0 = [0,0];

    % 1.3 Simulate
    effic = [1-eta, 1+eta];
    
    for effic_i = 1:2
        odefun  = @(t, x) xdot_1D(t,x,m, m_dot, @(t) Ft(t)*effic(effic_i), Cd0, Aref0);
        [t1, x1]  = ode45(odefun, tSpan, x0);
        x10(effic_i) = x1(end,1);
        v10(effic_i) = x1(end,2);
    end
    
    % 1.4 Display extreme values
    display('*Extreme cases*')
    display(['burn out altitudes: ' num2str(x10) ' m']);
    display(['burn out velocitys: ' num2str(v10) ' m/s']);
    display(' ');
    
    % 1.5 Compose final conditions matrices
    X10 = linspace(x10(1), x10(2), refine);
    V10 = linspace(v10(1), v10(2), refine);
    [X10, V10] = meshgrid(X10, V10);
    
    %% 2. determine t2 numerically
    
    % test numerical method
    [x02, res, max_flag] = fixpoint(@(x02) fn_x02(x02, X10(end, end), V10(end, end), xt),...
                                    (X10(end,end)+xt)/2, 0.1, 100);
    
    Z2 = 0;

    
    %% 3. nested functions
    
    function x02_out = fn_x02(x02, x01, v01, x3)
        
        % constantes
        g       =   9.6;
        rho     =   1.2;
        
        % variables de normalisation
        Vinf1   =   sqrt(2*(m0-mp)*g/(rho*Cd0*Aref0));
        tau1    =   Vinf1/g; 
        Vinf2   =   sqrt(2*(m0-mp)*g/(rho*(Cdb*Arefb+Cd0*Aref0)));
        
        % normalisation des variables
        v01n    =   v01/Vinf1;
        
        % valeurs de la fonction
        V02     = Vinf2*sqrt((exp(2*g*(x3-x02)/Vinf2^2)-1));
        t02     = tau1*(atan(v01/Vinf1)-atan(V02/Vinf1));
        x02_out = x01 + Vinf1^2/g*log((cos(atan(v01/Vinf1)-t02/tau1))/(cos(atan(v01/Vinf1))));
        
    end

    function [x,res, max_flag] = fixpoint(phi, x0, resmax, Nmax)
        x = x0;
        xold = x0;
        max_flag = 0;
        for n = 1:Nmax
            x = phi(x);
            res = abs(xold-x);
            xold = x;
            if res < resmax
               break; 
            elseif n >=Nmax
               max_flag = 1;
               break;
            end           
        end
    end
end

