% Test_tail.m
% Test the tail function in the Rocket.m class
% A copy of the function is in tail.m in the same directory.
% Eric Brunner 03/04/2017

clear all;
close all;

% Values to test
D1 = 0.15; % [m] diameter fore of the tail.
D2 = 0.1:0.01:0.149999; % [m] diameter aft of tail.
L = 0.05:0.01:0.3; % [m] length of tail.
e = 0.02;   % [m] thikness of tail.
rho = 2700 ; % [kg/m^3]

% make data and plot it
figure; hold on;
title('Iz vs. D2 and L for tail');
for D2i = 1:length(D2)
   for Li = 1:length(L)
   Tail =  tail(D1, D2(D2i), L(Li), e, rho);
   Iz(D2i, Li) = Tail.Iz;
   end
   plot(L, Iz(D2i,:), 'DisplayName', ['D2 = ' num2str(D2(D2i))]);
end
xlabel('L');
ylabel('m');
lgd = legend('show');
set(lgd, 'Location', 'northwest');

figure; hold on;
title('Ir vs. D2 and L for stage');
for D2i = 1:length(D2)
   for Li = 1:length(L)
   Stage =  stage(L(Li), D1, e, rho);
   Ir(D2i, Li) = Stage.Ir;
   end
   plot(L, Ir(D2i,:), 'DisplayName', ['D = ' num2str(D2(D2i))]);
end
xlabel('L');
ylabel('Ir');
lgd = legend('show');
set(lgd, 'Location', 'northwest');