clear all
close all

alt = 0:100:10000;
[T, a, p, rho] = stdAtmos(alt);

figure
title 'self vs. Aerotoolbox standard atmosphere'

subplot(2,2,1)
title 'Temperature'
plot(alt, T, 'DisplayName', 'self');
hold on
xlabel 'h [m]'
ylabel 'T [K]'

subplot(2,2,2)
title 'Speed of sound'
plot(alt,a, 'DisplayName', 'self');
hold on
xlabel 'h [m]'
ylabel 'a [m/s]'

subplot(2,2,3)
title 'Pressure'
plot(alt, p, 'DisplayName', 'self');
hold on
xlabel 'h [m]'
ylabel 'p [Pa]'

subplot(2,2,4)
title 'density'
plot(alt, rho, 'DisplayName', 'self');
hold on
xlabel 'h [m]'
ylabel '\rho [kg/m^3]'

[T, a, p, rho] = atmosisa(alt);

subplot(2,2,1)
plot(alt, T, 'DisplayName', 'aerotool');
legend show

subplot(2,2,2)
plot(alt, a, 'DisplayName', 'aerotool');
legend show

subplot(2,2,3)
plot(alt, p, 'DisplayName', 'aerotool');
legend show

subplot(2,2,4)
plot(alt, rho, 'DisplayName', 'aerotool');
legend show