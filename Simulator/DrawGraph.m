function A = DrawGraph( tsim, Xsim, R, S, alpha, calibre, T, M)
% DrawGraph
%   Affiche la serie de graphe selon l'envie de l'utilisateur. Le choix est
%   dans l'INPUT S.G(i)
%   Notice : Color star in red is the end of the thrust
%            Color star in blue is the maximum hight
%   INPUTS :
%       - tsim  : etapes de temps d'integration
%       - Xsim  : etat du system a chaque etape
%       - R     : objet 'Rocket'
%       - S     : objet 'Simulation'
%       - alpha : angle d'attaque
%       - calibre 
%               : calibres entre le CM et le CP 
%       - T     : shear values at query points and query times
%       - M     : flexion values at query points and query times
%   OUTPUTS :
%       - A     : Ne retourne rien

    %Valeur des certain point
    i_burnout = min(find(tsim>R.Motor.bt));%Permet d'avoir l'indice de fin du moteur
    i_finSim = max(find(Xsim(:,2)));
    
    
    if S.G(1) == true
    t = linspace(0, max([R.Motor.bt])+15, 100);
    for it = 1:length(t) 
       mass(it) = R.m(t(it));
       cm(it) = R.cm(t(it));
       Iz(it) = R.Iz(t(it));
       Ir(it) = R.Ir(t(it));
    end
    figure
    subplot(5,1,1)
    plot(t, mass);
    title('Mass properties');
    ylabel('m [kg]');
    subplot(5,1,2)
    plot(t, cm);
    ylabel('CM [m]');
    subplot(5,1,3)
    plot(t, Iz);
    ylabel('Iz [kg*m^2]');
    subplot(5,1,4)
    plot(t, Ir);
    ylabel('Ir [kg*m^2]');
    subplot(5, 1, 5)
    plot(R.Motor.ThrustCurve(1:end-1,1), R.Motor.ThrustCurve(1:end-1,2));  
    ylabel('Ft [N]');
    xlabel('t [s]');
    end % fin du if 1
    
    if S.G(2) == true
    alphaAvantSimu = linspace(-30, 30, 100)*pi/180;
    figure
    for alpha_i = 1:length(alphaAvantSimu)
       [CNa(alpha_i), zCP(alpha_i)] = R.aeroCoeff(alphaAvantSimu(alpha_i), 0, 0); 
    end
    subplot(2,1,1)
    plot(alphaAvantSimu*180/pi, CNa.*alphaAvantSimu);
    title('Aerodynamic properties at M = 0 & \theta = 0.')
    ylabel('CN');
    grid on;
    subplot(2,1,2)
    plot(alphaAvantSimu*180/pi, zCP);
    ylabel('zCP [m]');
    grid on;
    xlabel('alpha [^\circ]');
    end % fin du if2
    
    if S.G(3) == true
    drawRocket(R);
    end % fin du if 3
    
    if S.G(4) == true
    figure;hold on;
    title('Altitude')
    plot(tsim(i_burnout,1), Xsim(i_burnout,2), '*r');
    plot(tsim(i_finSim,1), Xsim(i_finSim,2), '*b');
    plot(tsim, Xsim(:,2));
    xlabel('t [s]');
    ylabel('z [m]');
    end % end if 4
    
    if S.G(5) == true
    figure; hold on;
    title('Rocket Trajectory')
    plot(Xsim(:,1), Xsim(:,2));
    plot(Xsim(i_burnout,1), Xsim(i_burnout,2), '*r');
    plot(Xsim(i_finSim,1), Xsim(i_finSim,2), '*b');
    daspect([1 1 1]);
    quiver(Xsim(1:100:end,1), Xsim(1:100:end,2), 0.001*sin(Xsim(1:100:end, 5)), 0.001*cos(Xsim(1:100:end, 5)));
    xlabel('horizontal distance from launchpad');
    ylabel('vertical altitude');
    end % end if 5
    
    if S.G(6) == true
    figure; hold on;
    title('Rocket angle');
    plot(tsim(i_burnout,1), Xsim(i_burnout,5), '*r');
    plot(tsim(i_finSim,1), Xsim(i_finSim,5), '*b');
    plot(tsim, Xsim(:,5));
    xlabel('t [s]');
    ylabel('\phi [rad]');
    end % end if 6
    
    if S.G(7) == true
    figure; hold on;
    title('Angle of attack');
    plot(tsim(i_burnout,1), alpha(i_burnout), '*r');
    plot(tsim(i_finSim,1), alpha(i_finSim), '*b');
    plot(tsim, alpha);
    xlabel('t [s]');
    ylabel('\alpha [rad]');
    end % end if 7
    
    if S.G(8) == true
    figure; hold on;
    title('Stabiltie statique lors du vol');
    plot(tsim(i_burnout,1), calibre(i_burnout), '*r');
    plot(tsim(i_finSim,1), calibre(i_finSim), '*b');
    plot(tsim, calibre);
    xlabel('t [s]');
    ylabel('Calibre [-]');
    end % end if 8

    if S.G(9) == true
    figure; hold on;
    title('Flexion load');
    plot(S.xquer, M);
    legend(['t = ' num2str(S.tquer(1))], ['t = ' num2str(S.tquer(2))], ['t = ' num2str(S.tquer(3))], ['t = ' num2str(S.tquer(4))], ['t = ' num2str(S.tquer(5))])
    xlabel('x [m]');
    ylabel('M(x) [Nm]');
    end % end if 9

end

