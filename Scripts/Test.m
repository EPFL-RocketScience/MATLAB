%% Exemple
clear all;
close all;
clc;

L_tot = 3.5;
g=9.81;
m = [3 7 10 10 4];
L = [0.8 0.7 1 0.8 2.7];
d = [0 0.8 1.2 2 0.8];
P = m*5*g./L;

xx = 0:0.01:L_tot;
lm = length(m);
lx = length(xx);
M_tot = zeros(lm,lx);
T_tot = zeros(lm,lx);

for i=1:lm
Ay = L(i)*P(i)*(1-(d(i)+L(i)/2)/L_tot);
By = L(i)*P(i)*(d(i)+L(i)/2)/L_tot;
T1 = []; T2 = []; T3 = [];
Mo1 =[]; Mo2 = []; Mo3 = [];

    for x=xx
        % x<d
        if x<=d(i)
        T1 = [T1 Ay];
        Mo1 = [Mo1 Ay*x];
        
        elseif x>d(i) && x<=d(i)+L(i)
        %pour d<x<d+L
        T2 = [T2 Ay-P(i)*(x-d(i))];
        Mo2 = [Mo2 Ay*x-P(i)*(x-d(i))*(x-d(i))/2];
        
        else
        %pour x>d+L
        T3 = [T3 Ay-P(i)*L(i)];
        Mo3 = [Mo3 Ay*x-P(i)*L(i)*(x-d(i)-L(i)/2)];
        end
    end

    T_tot(i,:) = [T1 T2 T3];
    M_tot(i,:) = [Mo1 Mo2 Mo3];      %valeur des moments pour chaque force

end

%calcul du moment total
Mfin = sum(M_tot);
Tfin = sum(T_tot);

figure(1)
plot(xx,Tfin,'b','DisplayName','T(x)')
hold on
plot(xx,Mfin,'r','DisplayName','M(x)')
hold off
legend show

%contrainte max apparait le plus loin de l'axe de flexion y = 0.075m
D_ext = 0.15;
e=[0.001 0.002 0.003];
I = pi/64*(D_ext^4-(D_ext-e).^4);
S = pi*D_ext*e;
ymax = D_ext/2;
ld=length(Mfin);
le=length(e);
sigma=zeros(le,ld);
tau = zeros(le,ld);
for i=1:length(e)
    sigma(i,:) = ymax*Mfin/I(i);
    tau(i,:) = 2*Tfin/S(i);
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