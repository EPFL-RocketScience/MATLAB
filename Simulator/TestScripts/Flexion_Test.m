clear all
close all

L_tot = 3.5;
g=9.81;
m = [3 7 10 10 4];
l = [0.8 0.7 1 0.8 2.7];
xs = [0 0.8 1.2 2 0.8];

xx = 0:0.01:L_tot;

[T, M] = Flexion(m', l', xs', xx, 5*g, L_tot);

figure(1)
plot(xx,T,'b','DisplayName','T(x)')
hold on
plot(xx,M,'r','DisplayName','M(x)')
hold off
legend show

%contrainte max apparait le plus loin de l'axe de flexion y = 0.075m
D_ext = 0.15;
e=[0.001 0.002 0.003];
I = pi/64*(D_ext^4-(D_ext-e).^4);
S = pi*D_ext*e;
ymax = D_ext/2;
ld=length(M);
le=length(e);
sigma=zeros(le,ld);
tau = zeros(le,ld);
for i=1:length(e)
    sigma(i,:) = ymax*M/I(i);
    tau(i,:) = 2*T/S(i);
end
figure(2)
plot(xx,sigma(1,:),'b','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=1')
hold on
grid on
plot(xx,sigma(2,:),'k','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=2')
plot(xx,sigma(3,:),'g','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=3')
plot(xx,tau(1,:),'r','LineWidth',1,'DisplayName','\tau_{max} (x)')
hold off
legend show