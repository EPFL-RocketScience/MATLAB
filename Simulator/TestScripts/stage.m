function Stage = stage(L, Dout, e, rho)

    Din = Dout-2*e;
    r = Dout/2;
    Stage.m = rho*2*pi*Dout/2*e*L;
    Stage.cm = L/2;
    Stage.Iz = Stage.m*(Dout/2)^2;
    Stage.Ir = pi*rho*L*e/12*(12*r^3-18*e*r^2+L^2*(2*r-e));
    
    % Calcule des proprietes de masse
    %probleme au niveau de l atribution aux differents etages
%     Stage.m = rho*pi*L*((Dout/2)^2-(Din/2)^2);
%     Stage.cm = L/2;
%     Stage.Iz = pi*rho*L*1/2*((Dout/2)^4-(Din/2)^4);
%     Stage.Ir = pi*rho*L*1/12*(3*((Dout/2)^4-(Din/2)^4)+L^2*((Dout/2)^2-(Din/2)^2));
end