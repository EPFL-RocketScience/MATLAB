function DrawGraph(R, S, results)
% DrawGraph
%   Affiche la serie de graphe selon l'envie de l'utilisateur. Le choix est
%   dans l'INPUT S.G(i)
%   Notice : Color star in red is the end of the thrust
%            Color star in blue is the maximum height
%   INPUTS :
%       - R     : objet 'Rocket'
%       - S     : objet 'Simulation'
%       - results : structure de resultats de la simulation. Retourne par
%                 la fonction simulate().
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition des variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tsim = results.tsim;
    Xsim = results.Xsim;
    alpha = results.alpha;
    calibre = results.calibre;
    M = results.flexion.M;


    %%%%%%%%%%%%%%%%%%%%%
    % Points importants
    %%%%%%%%%%%%%%%%%%%%%
    
    % fin de combustion
    i_burnout = min(find(results.tsim>R.Motor.bt));
    
    % apogee, ou conclusion de la simulation
    i_finSim = max(find(results.Xsim(:,2)));
    
    
    %%%%%%%%%%%%%%%%%%%%%
    % Plots
    %%%%%%%%%%%%%%%%%%%%%
    
    % Proprietes de masse et de poussee 
    if S.getGraphValue('mass') == true
        t = linspace(0, max([R.Motor.bt])+15, 100);
        for it = 1:length(t) 
           mass(it) = R.m(t(it));
           cm(it) = R.cm(t(it));
           Iz(it) = R.Iz(t(it));
           Ir(it) = R.Ir(t(it));
        end
        figure
        subplot(5,1,1)
        plot(t, mass, 'Linewidth', 1);
        %title('Mass properties');
        ylabel('m [kg]');
        set(gca,'FontSize',20);
        subplot(5,1,2)
        plot(t, cm, 'Linewidth', 1);
        ylabel('CM [m]');
        set(gca,'FontSize',20);
        subplot(5,1,3)
        plot(t, Iz, 'Linewidth', 1);
        ylabel('Iz [kg*m^2]');
        set(gca,'FontSize',20);
        subplot(5,1,4)
        plot(t, Ir, 'Linewidth', 1);
        ylabel('Ir [kg*m^2]');
        set(gca,'FontSize',20);
        subplot(5, 1, 5)
        plot(R.Motor.ThrustCurve(1:end-1,1), R.Motor.ThrustCurve(1:end-1,2), 'Linewidth', 1);  
        ylabel('Ft [N]');
        xlabel('t [s]');
        set(gca,'FontSize',20);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'mass','eps');
    end
    
    % Proprietes aerodynamiques
    if S.getGraphValue('aero') == true
        alphaAvantSimu = linspace(-30, 30, 100)*pi/180;
        figure
        for alpha_i = 1:length(alphaAvantSimu)
           [CNa(alpha_i), zCP(alpha_i)] = R.aeroCoeff(alphaAvantSimu(alpha_i), 0, 0); 
        end
        subplot(2,1,1)
        plot(alphaAvantSimu*180/pi, CNa.*alphaAvantSimu, 'Linewidth', 1);
        %title('Aerodynamic properties at M = 0 & \theta = 0.')
        ylabel('CN');
        set(gca,'FontSize',30);
        grid on;
        subplot(2,1,2)
        plot(alphaAvantSimu*180/pi, zCP, 'Linewidth', 1);
        ylabel('zCP [m]');
        grid on;
        xlabel('Angle d''attaque \alpha [^\circ]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'aero','eps')
    end
    
    % Dessin de la fusee
    if S.getGraphValue('rocket') == true
        drawRocket(R);
    end
    
    % Altitude vs. temps
    if S.getGraphValue('altitude') == true
        figure;hold on;
        %title('Altitude')
        plot(tsim(i_burnout,1), Xsim(i_burnout,2), '*r');
        plot(tsim(i_finSim,1), Xsim(i_finSim,2), '*b');
        plot(tsim, Xsim(:,2), 'Linewidth', 1);
        xlabel('temps t [s]');
        ylabel('Altitude z [m]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'altitude','eps')
    end
    
    % Trajectoire
    if S.getGraphValue('trajectory') == true
        figure; hold on;
        %title('Rocket Trajectory')
        plot(Xsim(:,1), Xsim(:,2), 'Linewidth', 1);
        plot(Xsim(i_burnout,1), Xsim(i_burnout,2), '*r');
        plot(Xsim(i_finSim,1), Xsim(i_finSim,2), '*b');
        daspect([1 1 1]);
        quiver(Xsim(1:100:end,1), Xsim(1:100:end,2), 0.001*sin(Xsim(1:100:end, 5)), 0.001*cos(Xsim(1:100:end, 5)));
        xlabel('Distance horizontal depuis la rampe de lancement');
        ylabel('Altitude z [m]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'trajectory','eps')
    end
    
    % Angle par rapport ? la verticale
    if S.getGraphValue('verticalAngle') == true
        figure; hold on;
        %title('Rocket angle');
        plot(tsim(i_burnout,1), Xsim(i_burnout,5), '*r');
        plot(tsim(i_finSim,1), Xsim(i_finSim,5), '*b');
        plot(tsim, Xsim(:,5), 'Linewidth', 1);
        xlabel('t [s]');
        ylabel('Angle de la fusée par rapport à la verticale \phi [rad]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'verticalAngle','eps')
    end
    
    % Angle d'attaque
    if S.getGraphValue('attackAngle') == true
        figure; hold on;
        %title('Angle of attack');
        plot(tsim(i_burnout,1), alpha(i_burnout), '*r');
        plot(tsim(i_finSim,1), alpha(i_finSim), '*b');
        plot(tsim, alpha, 'Linewidth', 1);
        xlabel('temps t [s]');
        ylabel('Angle d''attaque \alpha [rad]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'attackAngle','eps')
    end
    
    % Marge de stabilite
    if S.getGraphValue('stability') == true
        figure; hold on;
        %title('Stabiltie statique lors du vol');
        plot(tsim(i_burnout,1), calibre(i_burnout), '*r');
        plot(tsim(i_finSim,1), calibre(i_finSim), '*b');
        plot(tsim, calibre, 'Linewidth', 1);
        xlabel('temps t [s]');
        ylabel('Calibre [-]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'stability','eps')
    end

    % Moment de flexion
    if S.getGraphValue('flexion') == true
        figure; hold on;
        %title('Flexion load');
        plot(results.xquer, M, 'Linewidth', 1);
        legend(['t = ' num2str(S.tquer(1))], ['t = ' num2str(S.tquer(2))], ['t = ' num2str(S.tquer(3))], ['t = ' num2str(S.tquer(4))], ['t = ' num2str(S.tquer(5))])
        xlabel('x [m]');
        ylabel('M(x) [Nm]');
        set(gca,'FontSize',30);
        set(gcf, 'Position', [0, 0, 3840, 2160]);
        saveas(gcf,'flexion','eps')
    end

end

