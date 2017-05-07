function [Z2, Z10, V10, qual, err, err_max] = brakeMatrix( m0, mp, D, Arefb, Cd0, Cdb, TC, BT, xt)
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
%       - Z10   :   meshgrid of altitudes after burnout
%       - V10   :   meshgrid of velocities after burnout
%       - qual  :   quality of brake matrix, percentage of unusable data
%                   points
%       - err   :   matrix of possible errors on deployment altitude
%       - max_err
%               :   max error on deployment altitudes.

    %% 0. constants

    % 0.1 reference Area
    Aref0 = pi*D^2/4;
    
    % 0.2 motor uncertainty factor
    eta = 0.2;
    
    % 0.3 grid refinement, number of grid points
    refine = 50; 
    
    % 0.4 environment
    g       =   9.6;
    rho     =   1.2;

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
    Z10 = linspace(x10(1), x10(2), refine);
    V10 = linspace(v10(1), v10(2), refine);
    [Z10, V10] = meshgrid(Z10, V10);
    
    %% 2. populate deployment matrix Z2
    
    Z2 = zeros(length(Z10), length(V10));
    err = zeros(length(Z10), length(V10));
    
    for i = 1:length(Z10)
        for j = 1:length(V10)
            [Z2(i,j), err(i,j)] = deployAltitude(Z10(i,j), V10(i,j), xt);     
        end
    end

    %% 3. analyse results
    
    % 3.1 count invalid matrix points
    n_invalid   =   length(find(Z2 == -1));
    n_tot       =   prod(size(Z2));
    qual        =   n_invalid/n_tot;
    
    % 3.2 max possible error
    err_max = max(max(err));
    
    %% 4. nested functions
    
    function [x_deploy, err]  = deployAltitude(X10, V10, xt)
        % determine min and max altitude and time to apogee
        [x_max, t_max]  =   apogee(X10, V10, 1);
        [x_min, t_min]  =   apogee(X10, V10, 2);

        % check if target altitude is in the range
        if (xt>x_max)
            warning('Altitude:max', ['target altitude is too high. The max altitude is : ' num2str(x_max)]);
            x_deploy = -1;
            err      = -1;
            return
        elseif(xt<x_min)
            warning('Altitude:min', ['target altitude is too low. The min altitude is : ' num2str(x_min)]);
            x_deploy = -1;
            err      = -1;
            return
        end  
        
        % test all possible values of deployment times
        t2 = linspace(0, t_max,100);

        % find apogee for deployment time values
        [X02, V02]  = pos_speed(X10, V10, t2, 1);
        [X3, ~]     = apogee(X02, V02, 2);

        % find deployment time that gets the rocket nearest to the target
        t2_i = find(X3>=xt, 1, 'first');
        t_deploy = t2(t2_i); 

        % compute corresponding deployment altitude
        [x_deploy, ~] = pos_speed(X10, V10, t_deploy, 1);
        
        % estimate error
        if(t2_i == 0)
            [x_err, ~] = pos_speed(X10, V10, t2(t2_i+1), 1);
            err = abs(x_deploy - x_err);
        elseif(t2_i == length(t2))
            [x_err, ~] = pos_speed(X10, V10, t2(t2_i-1), 1);
            err = abs(x_deploy - x_err);
        else
            [x_errmin, ~] = pos_speed(X10, V10, t2(t2_i-1), 1);
            [x_errmax, ~] = pos_speed(X10, V10, t2(t2_i+1), 1);
            err = max([abs(x_deploy-x_errmin), abs(x_deploy-x_errmax)]);
        end
        
    end
    
    function [x, v] = pos_speed(x0, v0, t, phase)
        if phase == 1
            Vinf    =   sqrt(2*(m0-mp)*g/(rho*Cd0*Aref0));
            tau     =   Vinf/g; 
        elseif phase == 2
            Vinf    =   sqrt(2*(m0-mp)*g/(rho*(Cdb*Arefb+Cd0*Aref0)));
            tau     =   Vinf/g;
        else
            error('phase must either be 1 or 2')
        end
        
        V0  =   v0/Vinf;
        T   =   t/tau;
        
        v   =   Vinf*tan(atan(V0)-T);
        x   =   Vinf*tau*log(cos(atan(V0)-T)/cos(atan(V0)))+x0; 
    end

    function [x, t] = apogee(x0, v0, phase)
        
        if phase == 1
            Vinf    =   sqrt(2*(m0-mp)*g/(rho*Cd0*Aref0));
            tau     =   Vinf/g; 
        elseif phase == 2
            Vinf    =   sqrt(2*(m0-mp)*g/(rho*(Cdb*Arefb+Cd0*Aref0)));
            tau     =   Vinf/g;
        else
            error('phase must either be 1 or 2')
        end
        
        V0  =   v0/Vinf;
        
        t   =   atan(V0)*tau;
        x   =   -Vinf*tau*log(cos(atan(V0)))+x0;
    end
end

