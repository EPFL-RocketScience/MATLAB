function xdot = xdot_2D(t, x, m, m_dot, Ft, D, L, I, I_dot, d, Vp)

% function xdot
% State equation of rocket in 2D radial plane from it's departure point.  
%
% VARIABLES:
% m(t)      : mass
% m_dot(t)  : derivative of mass
% Ft(t)     : motor axial thrust force
% D(alpha)  : Drag
% L(alpha)  : Lift
% I(t)      : moment of inertia
% I_dot(t)  : derivative of moment of inertia
% d(t)      : distance between CG and CP
% Vp(Z)     : velocity of perturbed wind at given altitude

    % State variables
    R = x(1);
    Z = x(2)
    Vr = x(3);
    Vz = x(4);
    omega = x(5);
    theta = x(6);
    
    % Speed of CG in rocket coordinates
    Va = cosd(theta)*Vr + sind(theta)*Vz;
    Vn = -sind(theta)*Vr + cosd(theta)*Vz;
    
    % environment wind velocity in rocket coordinates
    [Vpr, Vpz] = Vp(Z);
    Vpa = cosd(theta)*Vpr + sind(theta)*Vpz;
    Vpn = -sind(theta)*Vpr + cosd(theta)*Vpz;
    
    % angle of attack
    alpha = atand((Vpn-Vn)/(Vpa-Va));
    
    % State equations
    Rdot = Vr;
    Zdot = Vz;
    Vrdot = 1/m(t)*(Ft(t)*cosd(theta) - D(alpha)*cosd(theta)...
                    -L(alpha)*sind(theta) + m_dot(t)*Vr);
    Vzdot = 1/m(t)*(Ft(t)*cosd(theta) - D(alpha)*sind(theta)... 
                    + L(alpha)*cosd(theta) - m(t)*g - m_dot(t)*Vz);
    omegadot = (-d*L(alpha) - I_dot(t)*omega)/I(t);
    thetadot = omega;
    
    % diff variable
    xdot = [Rdot, Zdot, Vrdot, Vzdot, omegadot, thetadot];
    
end