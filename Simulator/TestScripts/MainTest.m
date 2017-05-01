% Test Rocket class and simulator

close all;
clear all;

R = RocketJuju();
K = 1; % coefficient correctif de force aerodynamique normale

% Plot stuff
t = linspace(0, max([R.Motor.bt])+15, 100);
for it = 1:length(t) 
   mass(it) = R.m(t(it));
   cm(it) = R.cm(t(it));
   Iz(it) = R.Iz(t(it));
   Ir(it) = R.Ir(t(it));
end
figure
subplot(4,1,1)
plot(t, mass);
title('Mass properties');
ylabel('m [kg]');
subplot(4,1,2)
plot(t, cm);
ylabel('CM [m]');
subplot(4,1,3)
plot(t, Iz);
ylabel('Iz [kg*m^2]');
subplot(4,1,4)
plot(t, Ir);
ylabel('Ir [kg*m^2]');
xlabel('t [s]');


alpha = linspace(-30, 30, 100)*pi/180;
figure
for alpha_i = 1:length(alpha)
   [CNa(alpha_i), zCP(alpha_i)] = R.aeroCoeff(alpha(alpha_i), 0, 0, K); 
end
subplot(2,1,1)
plot(alpha*180/pi, CNa.*alpha);
title('Aerodynamic properties at M = 0 & \theta = 0.')
ylabel('CN');
grid on;
subplot(2,1,2)
plot(alpha*180/pi, zCP);
ylabel('zCP [m]');
grid on;
xlabel('alpha [^\circ]');

drawRocket(R);

% time step simulation
tspan = [R.Motor.ThrustCurve(1,1) max(t)];
tquer = [1 3 5 7 tspan(2)];
xquer = linspace(0, R.Tail.z + R.Tail.L, 10);
[tsim, Xsim, alpha, T, M] = Simulate( R, 1, K, tspan, tquer, xquer);

figure;hold on;
title('Altitude')
plot(tsim, Xsim(:,2));
xlabel('t [s]');
ylabel('z [m]');

figure; hold on;
title('Rocket Trajectory')
plot(Xsim(:,1), Xsim(:,2));
i_burnout = min(find(tsim>R.Motor.bt));
plot(Xsim(i_burnout,1), Xsim(i_burnout,2), '*r');
daspect([1 1 1]);
quiver(Xsim(1:100:end,1), Xsim(1:100:end,2), 0.001*sin(Xsim(1:100:end, 5)), 0.001*cos(Xsim(1:100:end, 5)));
xlabel('horizontal distance from launchpad');
ylabel('vertical altitude');

figure; hold on;
title('Rocket angle');
plot(tsim, Xsim(:,5));
xlabel('t [s]');
ylabel('\phi [rad]');

figure; hold on;
title('Angle of attack');
plot(tsim, alpha);
xlabel('t [s]');
ylabel('\alpha [rad]');

figure; hold on;
title('Flexion load');
plot(xquer, M);
legend(['t = ' num2str(tquer(1))], ['t = ' num2str(tquer(2))], ['t = ' num2str(tquer(3))], ['t = ' num2str(tquer(4))], ['t = ' num2str(tquer(5))])
xlabel('x [m]');
ylabel('M(x) [Nm]');