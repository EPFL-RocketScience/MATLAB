function xdot = xdot_1D(t, x, m, m_dot, Ft, Cd, S)

    % state variables
    z = x(1);
    v = x(2);

    [T, a, p, rho] = stdAtmos(z);
    
    % motion equations
    zdot = v;
    vdot = 1/m(t)*(Ft(t) - 0.5*S*rho*Cd*v^2 - m(t)*9.8 - m_dot(t)*v);
    
    xdot = [zdot; vdot];
end